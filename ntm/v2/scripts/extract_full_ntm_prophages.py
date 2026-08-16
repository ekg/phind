#!/usr/bin/env python3
"""
extract_full_ntm_prophages.py — NTM v2: extract every collaborator prophage from
the local PanSN bgzip genomes into a single FASTA with single-line sequences
(the format per_clade_alignment_pipeline.py's offset-index reader expects).

Reads the collaborator prophage master manifest (NOT geNomad provirus.tsv like
v1). One manifest row per prophage. Keeps has_prophage=True rows and applies the
GCA/GCF dedup rule decided in ntm-v2-acquire (prefer the GCA accession per numeric
assembly id; keep GCF-only numerics, stripping the NZ_ prefix from prophage_contig
when matching PanSN headers).

NCBI assemblies (GCA_/GCF_) are extracted from
  {v2}/genomes/canonical_objects/{acc}/{acc}.pansn.fa.gz  (PanSN header {acc}#1#{contig})
via samtools faidx, mirroring the v1 approach (ntm/scripts/extract_full_ntm_prophages.py).

Run-assemblies (ERR/SRR/DRR) have no local FASTAs yet (pending collaborator shipment
to ntm-v2-run); they are counted and skipped explicitly for a v2.1 extension.

Headering: collaborator prophage_ids are unique after dedup, so they are preserved
as the FASTA headers (stable for clade manifests). They have the form
{genome_id}_prophageN (e.g. GCA_000661085.1_..._genomic_prophage1).

Coordinates are 1-based inclusive (prophage_length_bp == end - start + 1, verified
in ntm-v2-acquire), so they are passed straight to samtools faidx (also 1-based
inclusive) with an off-by-one delta of 0.

Output: {v2}/full_prophages.fa
Report: {v2}/extract_report.md
"""
import argparse
import csv
import os
import re
import statistics
import subprocess
import time
from collections import Counter, defaultdict

V2 = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2"
MANIFEST = f"{V2}/inputs/NTM_QC_passed_prophage_master_manifest.tsv"
GENOMES = f"{V2}/genomes/canonical_objects"
OUT = f"{V2}/full_prophages.fa"
REPORT = f"{V2}/extract_report.md"

RUN_PREFIXES = ("ERR", "SRR", "DRR")
NCBI_PREFIXES = ("GCA_", "GCF_")

NUMERIC_RE = re.compile(r"^[A-Z]{3}_(\d+\.\d+)$")


def numeric(acc: str) -> str | None:
    m = NUMERIC_RE.match(acc)
    return m.group(1) if m else None


def pref(acc: str) -> str:
    if acc.startswith("GCA_"):
        return "GCA"
    if acc.startswith("GCF_"):
        return "GCF"
    m = re.match(r"^(ERR|SRR|DRR)", acc)
    return m.group(1) if m else acc


def choose_preferred_accession(accs: set[str]) -> str:
    """GCA-preferred per numeric id; else GCF; else the single run id."""
    if len(accs) == 1:
        return next(iter(accs))
    # same numeric -> pick GCA over GCF
    for a in accs:
        if a.startswith("GCA_"):
            return a
    for a in accs:
        if a.startswith("GCF_"):
            return a
    return next(iter(accs))


def load_fai_contigs(pansn: str) -> dict[str, tuple[str, int]]:
    """bare-contig-name -> (full PanSN contig name, length) from the .fai sidecar."""
    fai = pansn + ".fai"
    out = {}
    if not os.path.exists(fai):
        return out
    for line in open(fai):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        full = parts[0]
        name = full.split("#")[-1]
        try:
            out[name] = (full, int(parts[1]))
        except ValueError:
            continue
    return out


def resolve_contig(fai_contigs: dict[str, tuple[str, int]], manifest_contig: str):
    """Return (full_pansn_contig, length) for a manifest prophage_contig.

    Tries as-is first, then NZ_-stripped (covers GCF rows extracted from a GCA
    twin file where the local file dropped the NZ_ prefix).
    """
    if manifest_contig in fai_contigs:
        return fai_contigs[manifest_contig]
    if manifest_contig.startswith("NZ_"):
        stripped = manifest_contig[3:]
        if stripped in fai_contigs:
            return fai_contigs[stripped]
    return None


def faidx_region(pansn: str, contig: str, start: int, end: int) -> str | None:
    region = f"{contig}:{start}-{end}"
    res = subprocess.run(
        ["samtools", "faidx", pansn, region],
        capture_output=True, text=True)
    if res.returncode != 0 or not res.stdout.strip():
        return None
    lines = res.stdout.splitlines()
    seq = "".join(l for l in lines[1:] if l and not l.startswith(">"))
    return seq or None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default=MANIFEST)
    ap.add_argument("--genomes", default=GENOMES)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--report", default=REPORT)
    args = ap.parse_args()

    t0 = time.time()
    rows = list(csv.DictReader(open(args.manifest), delimiter="\t"))
    true_rows = [r for r in rows if r["has_prophage"] == "True"]

    # ---- dedup: prefer GCA per numeric id, keep all prophage rows of the
    # preferred accession, keep GCF-only, keep run-assemblies ----------------
    by_num = defaultdict(list)
    for r in true_rows:
        n = numeric(r["accession"])
        by_num[n if n else r["accession"]].append(r)
    kept = []
    for key, rl in by_num.items():
        accs = {r["accession"] for r in rl}
        pa = choose_preferred_accession(accs)
        kept.extend(r for r in rl if r["accession"] == pa)

    ncbi_rows = [r for r in kept if r["accession"].startswith(NCBI_PREFIXES)]
    run_rows = [r for r in kept if r["accession"].startswith(RUN_PREFIXES)]

    print(f"[input] total rows={len(rows)} has_prophage={len(true_rows)} "
          f"kept={len(kept)} ncbi={len(ncbi_rows)} run_asm(skip)={len(run_rows)}",
          flush=True)

    # precompute prophage_id uniqueness + dup-coord check
    pids = [r["prophage_id"] for r in ncbi_rows]
    dup_pids = len(pids) - len(set(pids))
    coord_key = Counter((r["accession"], r["prophage_contig"],
                         r["prophage_start"], r["prophage_end"]) for r in ncbi_rows)
    dup_coords = sum(1 for v in coord_key.values() if v > 1)

    n_ok = n_missing = n_err = 0
    missing_reasons = defaultdict(int)
    lens = []
    mismatches = []
    written = []
    start_zero_rows = []
    fai_cache = {}

    with open(args.out, "w") as fout:
        for r in ncbi_rows:
            acc = r["accession"]
            gid = r["genome_id"]
            pid = r["prophage_id"]
            start = int(float(r["prophage_start"]))
            end = int(float(r["prophage_end"]))
            clen = int(float(r["prophage_length_bp"]))
            contig = r["prophage_contig"]

            # Collaborator coordinates are 1-based inclusive, BUT a small number
            # of rows (13 kept NCBI) carry start=0. length_bp == end - start + 1
            # holds with start=0, i.e. start is 0-based there while end stays
            # 1-based inclusive. samtools faidx is 1-based inclusive, so for
            # start==0 we extract [1, end] and the extracted length is
            # length_bp - 1 (deliberate off-by-one, documented in the report).
            s1 = start if start >= 1 else 1
            if start == 0:
                start_zero_rows.append((pid, contig, end, clen, None))

            # source file: GCF row may need its GCA twin (GCA-preferred cohort)
            src_acc = acc
            pansn = f"{args.genomes}/{acc}/{acc}.pansn.fa.gz"
            if not os.path.exists(pansn) and acc.startswith("GCF_"):
                twin = "GCA_" + numeric(acc)
                twin_pansn = f"{args.genomes}/{twin}/{twin}.pansn.fa.gz"
                if os.path.exists(twin_pansn):
                    pansn, src_acc = twin_pansn, twin
            if not os.path.exists(pansn):
                missing_reasons[f"no_local_fasta_{acc}"] += 1
                n_missing += 1
                written.append((pid, "missing"))
                continue

            # cache fai per (src_acc, pansn)
            fkey = (src_acc, pansn)
            if fkey not in fai_cache:
                fai_cache[fkey] = load_fai_contigs(pansn)
            fai_contigs = fai_cache[fkey]

            resolved = resolve_contig(fai_contigs, contig)
            if resolved is None:
                missing_reasons[f"contig_not_found_{acc}"] += 1
                n_missing += 1
                written.append((pid, "missing"))
                continue
            fai_contig, fai_len = resolved

            seq = faidx_region(pansn, fai_contig, s1, end)
            if not seq:
                missing_reasons[f"faidx_empty_{acc}"] += 1
                n_err += 1
                written.append((pid, "faidx_empty"))
                continue

            if end > fai_len:
                reason = f"end_beyond_contig_{acc}"
            elif len(seq) != clen - (1 if start == 0 else 0):
                reason = f"length_mismatch_{acc}"
            elif len(seq) == 0:
                reason = f"empty_seq_{acc}"
            else:
                reason = None
            if reason is not None:
                missing_reasons[reason] += 1
                n_err += 1
                mismatches.append((pid, reason, len(seq), clen))
                written.append((pid, reason))
                continue

            fout.write(f">{pid}\n{seq}\n")
            n_ok += 1
            lens.append(len(seq))
            written.append((pid, "ok"))
            if start == 0 and start_zero_rows and start_zero_rows[-1][0] == pid:
                p, c, e, cl, _ = start_zero_rows[-1]
                start_zero_rows[-1] = (p, c, e, cl, len(seq))
            if n_ok % 1000 == 0:
                print(f"  {n_ok} extracted ({time.time()-t0:.0f}s)", flush=True)

    total_bp = sum(lens)
    print(f"[done] ok={n_ok} missing={n_missing} err={n_err} "
          f"skipped_run_asm={len(run_rows)} -> {args.out}", flush=True)
    print(f"length: min {min(lens)} median {int(statistics.median(lens))} "
          f"mean {int(statistics.mean(lens))} max {max(lens)} total {total_bp}",
          flush=True)

    # ---- write report ------------------------------------------------------
    lines = [
        "# NTM v2 — prophage FASTA extraction report",
        "",
        f"Generated: {time.strftime('%Y-%m-%dT%H:%M:%S')}Z",
        "",
        "## Counts",
        "",
        "| metric | count |",
        "|---|---|",
        f"| manifest rows (excl. header) | {len(rows)} |",
        f"| has_prophage=True rows | {len(true_rows)} |",
        f"| kept after GCA-dedup | {len(kept)} |",
        f"| NCBI prophage rows (extract) | {len(ncbi_rows)} |",
        f"| run-assembly rows (skipped, pending v2.1) | {len(run_rows)} |",
        f"| NCBI extracted into FASTA | {n_ok} |",
        f"| NCBI missing (no local fasta / contig) | {n_missing} |",
        f"| NCBI error (faidx / mismatch) | {n_err} |",
        "",
        "## Length distribution (extracted NCBI prophages)",
        "",
        "| stat | bp |",
        "|---|---|",
        f"| min | {min(lens)} |",
        f"| median | {int(statistics.median(lens))} |",
        f"| mean | {int(statistics.mean(lens))} |",
        f"| max | {max(lens)} |",
        f"| total | {total_bp} |",
        "",
        "Coordinates are 1-based inclusive (end-start+1 == length_bp, verified in "
        "ntm-v2-acquire). Extraction used samtools faidx with the raw start/end "
        "(off-by-one delta = 0) except for the start=0 rows below, where start is "
        "0-based (length_bp == end - start + 1 still holds), so the region was "
        "extracted as [1, end] and the extracted length is length_bp - 1 ("
        "deliberate off-by-one). Every non-zero-start extracted sequence length "
        f"equals its prophage_length_bp. Columns checked: {len(coord_key)} unique "
        f"(acc,contig,start,end) keys; duplicate same-coord rows: {dup_coords}.",
        "",
        f"prophage_id duplicates among NCBI rows: {dup_pids}.",
        "",
        f"Rows with 0-based start (start==0, extracted as [1, end], len == length_bp - 1): "
        f"{len(start_zero_rows)}",
    ]
    if start_zero_rows:
        lines += [""]
        lines += ["| prophage_id | contig | end | length_bp | extracted |"]
        lines += ["|---|---|---|---|---|"]
        for pid, contig, end, clen, elen in start_zero_rows:
            lines.append(f"| {pid} | {contig} | {end} | {clen} | {elen} |")
    lines += [
        "",
        "## Missing / error rows",
        "",
    ]
    if missing_reasons:
        lines += ["| reason | count |", "|---|---|"]
        for reason, c in sorted(missing_reasons.items()):
            lines.append(f"| {reason} | {c} |")
    else:
        lines.append("None — every NCBI prophage extracted.")
    lines += [
        "",
        "## Run-assembly skip",
        "",
        f"{len(run_rows)} run-assembly prophage rows (ERR/SRR/DRR) were skipped "
        "because their FASTAs are not yet local (pending the collaborator shipping "
        "them to ntm-v2-run; they arrive in a v2.1 extension)."
        "",
    ]
    if mismatches:
        lines += ["## Length mismatches (should be empty)", ""]
        for pid, reason, got, want in mismatches[:20]:
            lines.append(f"- {pid}: {reason} got={got} want={want}")
    with open(args.report, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    with open(f"{args.out}.manifest.tsv", "w") as fh:
        fh.write("prophage_id\tstatus\n")
        for pid, status in written:
            fh.write(f"{pid}\t{status}\n")

    print(f"[report] wrote {args.report} and {args.out}.manifest.tsv", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())