#!/usr/bin/env python3
"""
NTM v3 — extract full_prophages.fa from the v3 canonical objects using the
unified prophage manifest (build_v3_prophage_manifest.py output).

Prior art: ntm/v2/scripts/extract_full_ntm_prophages.py (v2), ntm/scripts/
extract_full_ntm_prophages.py (v1). Conventions carried over:

  * coordinates are used as given (no caller re-run); begin/end are 1-based
    inclusive, so the samtools faidx region is begin-end with off-by-one
    delta 0 — EXCEPT begin==0 rows (52 in the manifest, 21 extractable),
    where begin is 0-based while end stays 1-based inclusive: extracted as
    [1, end], extracted length = length - 1 (deliberate documented
    off-by-one, same convention and same assemblies as the v2 report table).
  * sequences are written single-line (per_clade_alignment_pipeline.py's
    offset-index reader format).
  * the transposable flag is metadata only — it does not shift coordinates
    (verified: every extracted length equals end-begin+1 off the manifest).

FASTA headers are PanSN: {canonical_acc}#1#{prophage_id} — the canonical
object name (the sample the sequence was physically read from), not the
source-manifest genome_id, so twin/BV-BRC variant genomes collapse onto one
sample namespace with no duplicate headers.

Scaffold -> contig resolution against the object .fai (v2 + acquisition
conventions, normalized both sides): exact bare name -> strip `accn|` ->
unversioned -> NZ_/NW_ prefix ADDED (RefSeq-renamed objects) or STRIPPED
(GCA twin objects). Every extractable row resolves (0 unresolved, 0
end-beyond-contig in the final run).

Outputs (NVMe):
  {work}/ntm/v3/full_prophages.fa              single-line FASTA
  {work}/ntm/v3/full_prophages.fa.manifest.tsv per-record status
Outputs (repo):
  ntm/v3/extract_report.md
  ntm/v3/inputs/coverage.tsv

Usage:
  python3 ntm/v3/scripts/extract_v3_prophages.py \
      [--work-dir /mnt/nvme3n1/erikg/phind-genome-work]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import random
import re
import statistics
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
MANIFEST_GZ = REPO / "ntm/v3/inputs/v3_prophage_manifest.tsv.gz"
STATS_JSON_NAME = "ntm/v3/scratch/v3_prophage_reconciliation.json"

UNVERSION_RE = re.compile(r"\.\d+$")
RUN_RE = re.compile(r"^(ERR|SRR|DRR)\d+$")
LENGTH_MIN, LENGTH_MAX = 1000, 100000


def unversion(tok: str) -> str:
    return UNVERSION_RE.sub("", tok)


def load_fai(pansn: str) -> dict[str, dict]:
    """resolution maps from the .fai sidecar (bare contig name keyed)."""
    fai = pansn + ".fai"
    exact: dict[str, tuple[str, int]] = {}
    unv: dict[str, list[str]] = defaultdict(list)
    with open(fai) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            full = parts[0]
            bare = full.split("#")[-1]
            try:
                ln = int(parts[1])
            except ValueError:
                continue
            exact[bare] = (full, ln)
            unv[unversion(bare)].append(bare)
    return {"exact": exact, "unv": unv}


def resolve_scaffold(fai: dict, scaffold: str):
    """scaffold -> (full_pansn_contig, contig_len, method) or (None,0,reason).

    Candidates tried in order: the bare scaffold as-is, then NZ_/NW_
    PREFIXED (the local RefSeq object renames contigs the export lists
    unprefixed — measured: the 72 AMRA/KB…/AP…/CP… scaffolds), then the
    NZ_/NW_ prefix STRIPPED (the export lists NZ_ contigs the local GCA
    twin object stores unprefixed). Exact name first, unversioned second.
    """
    s0 = scaffold
    s = s0[len("accn|"):] if s0.startswith("accn|") else s0
    suffix = "+accn_strip" if s0.startswith("accn|") else ""
    for tag, c in (("", s), ("nz_add", "NZ_" + s), ("nw_add", "NW_" + s)):
        if c in fai["exact"]:
            return (*fai["exact"][c], (tag or "exact") + suffix)
        u = fai["unv"].get(unversion(c))
        if u:
            pick = sorted(u)[0]
            return (*fai["exact"][pick],
                    (tag + "+version_suffix" if tag else "version_suffix") + suffix)
    for pfx in ("NZ_", "NW_"):
        if s.startswith(pfx):
            t = s[len(pfx):]
            if t in fai["exact"]:
                return (*fai["exact"][t], f"{pfx.rstrip('_').lower()}_strip" + suffix)
            u = fai["unv"].get(unversion(t))
            if u:
                pick = sorted(u)[0]
                return (*fai["exact"][pick],
                        f"{pfx.rstrip('_').lower()}_strip+version_suffix" + suffix)
    return None, 0, "unresolved"


def faidx(pansn: str, region: str) -> str | None:
    res = subprocess.run(["samtools", "faidx", pansn, region],
                         capture_output=True, text=True)
    if res.returncode != 0 or not res.stdout.strip():
        return None
    seq = "".join(l for l in res.stdout.splitlines()[1:]
                  if l and not l.startswith(">"))
    return seq or None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", default="/mnt/nvme3n1/erikg/phind-genome-work")
    args = ap.parse_args()
    work = Path(args.work_dir)
    genomes = work / "ntm/v3/genomes/canonical_objects"
    out_fa = work / "ntm/v3/full_prophages.fa"
    out_status = work / "ntm/v3/full_prophages.fa.manifest.tsv"
    report_path = REPO / "ntm/v3/extract_report.md"
    coverage_path = REPO / "ntm/v3/inputs/coverage.tsv"
    recon = {"work": str(work)}  # replaced below if present
    stats_file = work / STATS_JSON_NAME
    if stats_file.exists():
        import json
        recon = json.loads(stats_file.read_text())

    t0 = time.time()
    with gzip.open(MANIFEST_GZ, "rt", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    extractable = [r for r in rows if r["genome_status"] == "extractable"]
    blocked = [r for r in rows if r["genome_status"] == "blocked_run"]
    assert len(rows) == len(extractable) + len(blocked)

    by_obj: dict[str, list[dict]] = defaultdict(list)
    for r in extractable:
        by_obj[r["canonical_acc"]].append(r)

    n_ok = n_err = 0
    err_reasons: Counter = Counter()
    resolve_methods: Counter = Counter()
    lens: list[int] = []
    begin_zero: list[tuple] = []
    out_of_range: list[dict] = []
    end_beyond: list[dict] = []
    written: list[tuple[str, str]] = []   # (header, status)

    with open(out_fa, "w") as fout:
        for i, (obj, rs) in enumerate(sorted(by_obj.items()), 1):
            pansn = f"{genomes}/{obj}/{obj}.pansn.fa.gz"
            fai = load_fai(pansn)
            for r in sorted(rs, key=lambda r: (r["scaffold"], int(r["begin"]),
                                               int(r["end"]), r["prophage_id"])):
                header = f"{obj}#1#{r['prophage_id']}"
                begin, end, length = int(r["begin"]), int(r["end"]), int(r["length"])
                contig, clen, method = resolve_scaffold(fai, r["scaffold"])
                if contig is None:
                    err_reasons[f"scaffold_{method}"] += 1
                    n_err += 1
                    written.append((header, f"error:scaffold_{method}"))
                    continue
                resolve_methods[method] += 1
                if end > clen:
                    end_beyond.append({**r, "contig_len": clen})
                    err_reasons["end_beyond_contig"] += 1
                    n_err += 1
                    written.append((header, "error:end_beyond_contig"))
                    continue
                if begin == 0:
                    begin_zero.append((header, r["genome_id"], r["prophage_id"],
                                       r["scaffold"], end, length))
                seq = faidx(pansn, f"{contig}:{max(1, begin)}-{end}")
                if not seq:
                    err_reasons["faidx_empty"] += 1
                    n_err += 1
                    written.append((header, "error:faidx_empty"))
                    continue
                want = length - (1 if begin == 0 else 0)
                if len(seq) != want:
                    err_reasons["length_mismatch"] += 1
                    n_err += 1
                    written.append((header, "error:length_mismatch"))
                    continue
                fout.write(f">{header}\n{seq}\n")
                written.append((header, "ok"))
                n_ok += 1
                lens.append(len(seq))
                if not (LENGTH_MIN <= len(seq) <= LENGTH_MAX):
                    out_of_range.append({**r, "extracted_len": len(seq)})
            if i % 500 == 0:
                print(f"  {i}/{len(by_obj)} objects, {n_ok} extracted "
                      f"({time.time()-t0:.0f}s)", flush=True)

    print(f"[extract] ok={n_ok} err={n_err} of {len(extractable)} extractable "
          f"({time.time()-t0:.0f}s)", flush=True)
    assert n_ok == len(extractable) - n_err
    with open(out_status, "w") as fst:
        fst.write("header\tstatus\n")
        for header, status in written:
            fst.write(f"{header}\t{status}\n")
    assert len(written) == len(extractable)

    # ------------------------------------------------ independent spot check
    rng = random.Random(42)
    spot_objs = rng.sample(sorted(by_obj), 10)
    spot_rows = []
    fasta_records: dict[str, str] = {}
    with open(out_fa) as fh:
        hdr = None
        for line in fh:
            if line.startswith(">"):
                hdr = line[1:].rstrip("\n")
                fasta_records[hdr] = ""
            elif hdr:
                fasta_records[hdr] += line.rstrip("\n")
    spot_ok = 0
    for obj in spot_objs:
        r = rng.choice(by_obj[obj])
        pansn = f"{genomes}/{obj}/{obj}.pansn.fa.gz"
        fai = load_fai(pansn)
        contig, clen, method = resolve_scaffold(fai, r["scaffold"])
        begin, end = int(r["begin"]), int(r["end"])
        seq = faidx(pansn, f"{contig}:{max(1, begin)}-{end}") if contig else None
        header = f"{obj}#1#{r['prophage_id']}"
        match = (seq is not None and fasta_records.get(header) == seq
                 and len(seq) == int(r["length"]) - (1 if begin == 0 else 0))
        spot_rows.append((header, r["scaffold"], f"{contig}", begin, end,
                          len(seq or ""), "PASS" if match else "FAIL"))
        spot_ok += bool(match)
    assert spot_ok == 10, spot_rows

    # ------------------------------------------------------------ sha256
    h = hashlib.sha256()
    with open(out_fa, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    sha = h.hexdigest()
    n_records = len(fasta_records)
    assert n_records == n_ok, (n_records, n_ok)
    assert n_records == len(extractable) - n_err

    # blocked rows must all be run assemblies
    blocked_objs = {r["canonical_acc"] for r in blocked}
    assert all(RUN_RE.match(o) for o in blocked_objs), "non-run blocked row"

    total_bp = sum(lens)
    stats = {
        "min": min(lens), "median": int(statistics.median(lens)),
        "mean": int(statistics.mean(lens)), "max": max(lens),
    }

    # ---------------------------------------------------------- report
    def m(v: int) -> str:
        return f"{v:,}"

    by_src = Counter(r["source"] for r in rows)
    ext_src = Counter(r["source"] for r in extractable)
    blk_src = Counter(r["source"] for r in blocked)
    v2_run_manifest = sum(1 for r in rows
                          if r["source"] == "V2" and RUN_RE.match(r["accession"]))
    lines = [
        "# NTM v3 — unified prophage manifest + full_prophages.fa extraction report",
        "",
        f"Generated: {time.strftime('%Y-%m-%dT%H:%M:%SZ')}",
        "",
        "## Design",
        "",
        "v3 is all-inclusive: BV-BRC phigaro coordinates are primary for genomes",
        "present in the 2026-09-09 BV-BRC phigaro QC-passed export (any row",
        "namespace: ASSEMBLY runs, NCBI assemblies, BV-BRC taxid.version); v2",
        "collaborator calls are kept only for genomes absent from the export",
        "(v2-only NCBI assemblies + the 135 v2-only runs). No genome carries",
        "calls from both callers — on overlap the v2 set is superseded (counted",
        "below, not emitted). Coordinates are used as given; the caller is never",
        "re-run (v2 convention).",
        "",
        "## Manifest reconciliation",
        "",
        "| metric | count |",
        "|---|---:|",
        f"| coordinates CSV rows | {m(recon.get('coordinates_csv_rows', len(rows)))} |",
        f"| stub rows (empty prophage_id) in CSV | {m(recon.get('stub_rows', 0))} |",
        f"| stub rows dropped as exact duplicates of named rows | {m(recon.get('stub_rows_dropped_as_duplicates', 0))} |",
        f"| stub rows kept with synthesized `<scaffold>_prophage<N>` id | {m(recon.get('stub_rows_kept_synthesized_pid', 0))} |",
        f"| unique prophages after stub dedup | {m(recon.get('unique_prophages_after_stub_dedup', 0))} |",
        f"| export genomes whose summary prophage_count double-counts | {m(recon.get('genomes_double_counted_in_summary', 0))} |",
        f"| canonical assemblies with >1 export genome/object | {m(recon.get('export_internal_multi_genome_groups', 0))} (GCA/GCF twin objects: {m(recon.get('export_internal_twin_assembly_groups', 0))}) |",
        f"| export genomes superseded export-internally | {m(recon.get('export_internal_superseded_genomes', 0))} |",
        f"| export rows superseded export-internally (one genome per assembly kept, GCA-preferred) | {m(recon.get('export_internal_superseded_rows', 0))} |",
        f"| export genomes kept (phigaro calls) | {m(recon.get('export_genomes_kept', 0))} |",
        f"| v2 manifest rows / has_prophage / kept after GCA-dedup | {m(recon.get('v2_manifest_rows', 0))} / {m(recon.get('v2_has_prophage_rows', 0))} / {m(recon.get('v2_kept_after_gca_dedup', 0))} |",
        f"| v2 numeric rows superseded (genome in export) | {m(recon.get('v2_superseded_numeric_rows', 0))} ({m(recon.get('v2_superseded_numeric_genomes', 0))} genomes) |",
        f"| v2 run rows superseded (run in export) | {m(recon.get('v2_superseded_run_rows', 0))} |",
        f"| v2 rows kept (genomes absent from export) | {m(recon.get('v2_kept_manifest_rows', 0))} ({m(recon.get('v2_kept_numeric_rows', 0))} NCBI + {m(recon.get('v2_kept_run_rows', 0))} runs) |",
        f"| **unified manifest rows** | **{m(len(rows))}** |",
        f"| — source=BV-BRC (caller=phigaro) | {m(by_src['BV-BRC'])} |",
        f"| — source=V2 (caller=v2-collab) | {m(by_src['V2'])} |",
        "",
        "Genomes in both cohorts (reconciliation basis): every v2 kept",
        f"numeric assembly whose 9-digit numeric is in the export ({m(recon.get('v2_superseded_numeric_genomes', 0))} prophage-bearing",
        "genomes) and every v2 run accession present in the export",
        f"({m(recon.get('v2_superseded_run_genomes', 0))}) had its v2 prophage set superseded by the phigaro set;",
        f"{m(recon.get('v2_kept_numeric_rows', 0))} v2 NCBI rows and {m(recon.get('v2_kept_run_rows', 0))} v2 run rows are the v2-only retention.",
        "",
        "## Extraction counts",
        "",
        "| metric | count |",
        "|---|---:|",
        f"| unified manifest rows | {m(len(rows))} |",
        f"| extractable rows (genome local) | {m(len(extractable))} (BV-BRC {m(ext_src['BV-BRC'])} + V2 {m(ext_src['V2'])}) |",
        f"| extracted into FASTA | {m(n_ok)} |",
        f"| errors (scaffold/end/length) | {m(n_err)} |",
        f"| blocked rows (run assemblies, collaborator pending) | {m(len(blocked))} (BV-BRC {m(blk_src['BV-BRC'])} + V2 {m(blk_src['V2'])}) |",
        f"| distinct canonical objects read | {m(len(by_obj))} |",
        f"| FASTA records written | {m(n_records)} == extractable − errors |",
        "",
        "Run-assembly accounting (all-inclusive, none dropped): every v2 run",
        f"prophage row is present in the union — {m(recon.get('v2_superseded_run_rows', 0))} via the export phigaro",
        f"rows (blocked_run: FASTAs pending collaborator delivery, ntm/v2/run_assemblies/",
        f"REQUEST.md PENDING) and {v2_run_manifest} as source=V2 rows (the 36",
        "prophage-bearing v2-only runs; the other 99 v2-only runs carry no",
        "prophages). ENA substitute assemblies are content-mismatched and are",
        "not used (v2 ena_backfill study).",
        "",
        "## Length distribution (extracted prophages)",
        "",
        "| stat | bp |",
        "|---|---:|",
        f"| min | {m(stats['min'])} |",
        f"| median | {m(stats['median'])} |",
        f"| mean | {m(stats['mean'])} |",
        f"| max | {m(stats['max'])} |",
        f"| total | {m(total_bp)} |",
        "",
        "## Coordinate conventions",
        "",
        "Coordinates are 1-based inclusive: `length == end − begin + 1` holds for",
        f"all {m(len(rows))} manifest rows (asserted at build time), and every extracted",
        "sequence length equals its manifest length — except the begin==0 rows",
        f"({len(begin_zero)} extractable of 52 manifest-wide, same assemblies and",
        "convention as the v2 extract report): there begin is 0-based while end",
        "stays 1-based inclusive, so the region is extracted as [1, end] and the",
        "extracted length is length − 1. The transposable flag is metadata only",
        "and does not shift coordinates (verified across every row: extracted",
        "length == end − max(1,begin) + 1).",
        "",
        "Scaffold→contig resolution methods (normalized both sides:",
        "accn\| stripped, unversioned, NZ_/NW_ prefix added or stripped — both",
        "directions occur: RefSeq objects rename contigs the export lists",
        "unprefixed, and GCA twin objects store NZ_ contigs unprefixed):",
        "",
        "| method | rows |",
        "|---|---:|",
    ]
    for meth, c in sorted(resolve_methods.items(), key=lambda kv: -kv[1]):
        lines.append(f"| {meth} | {m(c)} |")
    lines += [
        "",
        "## Flagged warnings (kept, not dropped)",
        "",
        f"Prophages outside [{LENGTH_MIN:,}, {LENGTH_MAX:,}] bp: "
        f"**{len(out_of_range)}** — all extracted and retained:",
        "",
        "| prophage | source | length |",
        "|---|---|---:|",
    ]
    for r in out_of_range[:50]:
        lines.append(f"| {r['canonical_acc']}#1#{r['prophage_id']} | {r['source']} | {m(r['extracted_len'])} |")
    if len(out_of_range) > 50:
        lines.append(f"| … {len(out_of_range) - 50} more (see manifest; filter length column) | | |")
    lines += [
        "",
        f"end > contig length rows: **{len(end_beyond)}**",
        "",
    ]
    if end_beyond:
        lines += ["| prophage | end | contig_len |", "|---|---:|---:|"]
        for r in end_beyond[:50]:
            lines.append(f"| {r['canonical_acc']}#1#{r['prophage_id']} | {m(int(r['end']))} | {m(r['contig_len'])} |")
    else:
        lines.append("None — every prophage fits inside its resolved contig.")
    lines += [
        "",
        "## begin==0 rows (extracted as [1, end], extracted length = length − 1)",
        "",
        "| header | genome | scaffold | end | length |",
        "|---|---|---|---:|---:|",
    ]
    for header, gid, pid, scaf, end, length in begin_zero:
        lines.append(f"| {header} | {gid} | {scaf} | {m(end)} | {m(length)} |")
    lines += [
        "",
        "## Independent spot check (10 genomes, seed 42)",
        "",
        "Each row re-queried via a fresh `samtools faidx` call and byte-compared",
        "to the record parsed back out of the written FASTA:",
        "",
        "| header | scaffold | contig | begin | end | len | result |",
        "|---|---|---|---:|---:|---:|---|",
    ]
    for header, scaf, contig, begin, end, ln, res in spot_rows:
        lines.append(f"| {header} | {scaf} | {contig} | {m(begin)} | {m(end)} | {m(ln)} | {res} |")
    lines += [
        "",
        "## Output",
        "",
        f"- FASTA: `{out_fa}` — {m(n_records)} records, {m(total_bp)} bp, sha256 `{sha}`",
        f"- per-record status: `{out_status}`",
        f"- manifest (repo): `ntm/v3/inputs/v3_prophage_manifest.tsv.gz` "
        f"({m(len(rows))} rows; plain copy on NVMe at "
        "`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/inputs/v3_prophage_manifest.tsv`)",
        "",
        "FASTA stays on NVMe (repo holds only code, manifests and small",
        "reports — v1/v2 rule).",
        "",
    ]
    report_path.write_text("\n".join(lines) + "\n")

    # ------------------------------------------------------- coverage.tsv
    cov = [
        ("metric", "value", "detail"),
        ("manifest_rows", len(rows), "unified v3 prophage manifest (BV-BRC ∪ V2)"),
        ("manifest_rows_bvbrc", by_src["BV-BRC"], "caller=phigaro (2026-09-09 export)"),
        ("manifest_rows_v2", by_src["V2"], "caller=v2-collab (genomes absent from export)"),
        ("coords_csv_rows", recon.get("coordinates_csv_rows", 0), "ntm_qc_passed_phigaro_coordinates_20260909.csv"),
        ("stub_rows", recon.get("stub_rows", 0), "empty prophage_id rows in CSV"),
        ("stub_rows_dup_dropped", recon.get("stub_rows_dropped_as_duplicates", 0), "exact (genome,scaffold,begin,end) duplicate of a named row"),
        ("stub_rows_kept_synthesized", recon.get("stub_rows_kept_synthesized_pid", 0), "only representation; pid = <scaffold>_prophage<N>"),
        ("unique_prophages_after_stub_dedup", recon.get("unique_prophages_after_stub_dedup", 0), "42,067 expected"),
        ("export_internal_multi_genome_groups", recon.get("export_internal_multi_genome_groups", 0), "assemblies with >1 export genome/object (GCA/GCF twins, BV-BRC taxid variants)"),
        ("export_internal_twin_assembly_groups", recon.get("export_internal_twin_assembly_groups", 0), "assemblies held as both GCA_ and GCF_ canonical objects"),
        ("export_internal_superseded_genomes", recon.get("export_internal_superseded_genomes", 0), "export genomes dropped (one kept per assembly, GCA-preferred v2 convention)"),
        ("export_internal_superseded_rows", recon.get("export_internal_superseded_rows", 0), "prophage rows on superseded export genomes"),
        ("v2_rows_kept_total", recon.get("v2_kept_manifest_rows", 0), "v2 genomes absent from export (536 NCBI + 83 runs)"),
        ("v2_superseded_numeric_rows", recon.get("v2_superseded_numeric_rows", 0), "v2 rows whose genome is in the export (phigaro primary)"),
        ("v2_superseded_run_rows", recon.get("v2_superseded_run_rows", 0), "v2 run rows whose run is in the export (phigaro primary)"),
        ("extractable_rows", len(extractable), "genome local on NVMe"),
        ("extracted_records", n_ok, "FASTA records written == extractable − errors"),
        ("extraction_errors", n_err, "target 0"),
        ("blocked_run_rows", len(blocked), "ERR/SRR/DRR pending collaborator delivery; none dropped"),
        ("v2_run_rows_in_union", recon.get("v2_superseded_run_rows", 0) + v2_run_manifest, f"{recon.get('v2_superseded_run_rows', 0)} via export phigaro + {v2_run_manifest} source=V2"),
        ("distinct_objects_read", len(by_obj), "v3/genomes/canonical_objects"),
        ("total_bp", total_bp, "sum of extracted prophage lengths"),
        ("len_min", stats["min"], ""),
        ("len_median", stats["median"], ""),
        ("len_mean", stats["mean"], ""),
        ("len_max", stats["max"], ""),
        ("length_equals_end_minus_begin_plus_1", len(rows), "holds for all manifest rows (asserted)"),
        ("begin_zero_rows_extracted", len(begin_zero), "extracted as [1,end]; extracted = length - 1 (documented off-by-one)"),
        ("out_of_range_rows", len(out_of_range), f"outside [{LENGTH_MIN},{LENGTH_MAX}] bp; flagged, kept"),
        ("end_beyond_contig_rows", len(end_beyond), "target 0"),
        ("spot_check_genomes", 10, f"{spot_ok}/10 PASS (seed 42, independent faidx re-query)"),
        ("fasta_sha256", sha, str(out_fa)),
        ("fasta_records", n_records, "== extracted_records"),
        ("fasta_single_source_per_genome", len(rows), "no genome carries calls from two callers (asserted at build)"),
    ]
    with open(coverage_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for row in cov:
            w.writerow(row)
    print(f"[report] wrote {report_path}")
    print(f"[report] wrote {coverage_path}")
    print(f"[fasta] {out_fa} sha256={sha}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
