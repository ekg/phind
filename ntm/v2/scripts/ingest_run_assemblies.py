#!/usr/bin/env python3
"""Ingest collaborator run-assembly FASTAs into the v2 PanSN bgzip layout.

For every accession in the collaborator run-assembly cohort:
  1. locate  {staging}/{acc}/contigs.fasta
  2. (optional) verify sha256 against the collaborator-provided
     sha256sums.txt (transfer-integrity check)
  3. STRONG CONTENT CHECK: per-genome contig count and total length must
     equal the manifest columns `assembly_contigs` / `assembly_length_bp`
  4. convert to PanSN bgzip: {out}/{acc}/{acc}.pansn.fa.gz + .fai + .gzi
     with header {acc}#1#{exact_source_contig_token} (token preserved
     byte-for-byte — manifest prophage_contig values are exact source
     tokens, e.g. NODE_6_length_418114_cov_19.254028, and must remain
     addressable as {acc}#1#{prophage_contig} for extraction)
  5. samtools faidx validates the bgzip (also writes .fai/.gzi)

A per-genome report table is written (contigs + total bp vs manifest,
PASS/FAIL per genome) plus a machine summary JSON. Exit code 0 only if
every present accession passed all checks; mismatches are itemized.

Usage:
  python ntm/v2/scripts/ingest_run_assemblies.py \
      --manifest /mnt/.../NTM_QC_passed_prophage_master_manifest.tsv \
      --staging /mnt/.../ntm/v2/incoming/run_assemblies \
      --out /mnt/.../ntm/v2/genomes/run_assemblies \
      --accessions ntm/v2/run_assemblies/run_assemblies_needed.txt \
      [--checksums /mnt/.../ntm/v2/incoming/run_assemblies/sha256sums.txt] \
      [--report /mnt/.../ntm/v2/genomes/run_assemblies/ingest_report.md]
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

PANS_PREFIX_RE = re.compile(r"^>([^\s]+)")


def iter_records(fasta_path):
    """Yield (source_token, seq) from a plain or gzip FASTA."""
    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as fh:
        header = None
        seq = []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None and seq:
                    yield header, "".join(seq)
                m = PANS_PREFIX_RE.match(line)
                header = m.group(1) if m else line[1:].strip()
                seq = []
            else:
                seq.append(line.strip().replace(" ", ""))
        if header is not None and seq:
            yield header, "".join(seq)


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def load_manifest(manifest_path):
    """acc -> {assembly_contigs, assembly_length_bp, has_prophage} for runs."""
    out = {}
    with open(manifest_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            acc = r["accession"]
            if acc[:3] not in ("ERR", "SRR", "DRR"):
                continue
            key = int(r["assembly_contigs"])
            val = int(r["assembly_length_bp"])
            if acc in out:
                # all rows for one accession should agree; keep first, flag
                if (out[acc]["assembly_contigs"] != key
                        or out[acc]["assembly_length_bp"] != val):
                    out[acc]["conflict"] = True
            else:
                out[acc] = {
                    "assembly_contigs": key,
                    "assembly_length_bp": val,
                    "has_prophage": r["has_prophage"] == "True",
                }
    return out


def load_checksums(path):
    """sha256sums.txt -> {relpath: hexdigest}; normalizes './' and leading dirs."""
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # format: <hex> <sp> <path>  (path may be '*path' in binary mode)
            m = re.match(r"^([0-9a-f]{64})[ \t]\*?(.*)$", line)
            if not m:
                sys.stderr.write(f"unparseable checksum line: {line!r}\n")
                continue
            digest, rel = m.groups()
            rel = rel.lstrip("./")
            out[rel] = digest
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--staging", required=True,
                    help="dir with {acc}/contigs.fasta (post-fetch)")
    ap.add_argument("--out", required=True,
                    help="run_assemblies output dir (PanSN bgzip)")
    ap.add_argument("--accessions", required=True,
                    help="run_assemblies_needed.txt (order = priority)")
    ap.add_argument("--checksums", help="optional sha256sums.txt to verify")
    ap.add_argument("--report", help="markdown report output path")
    ap.add_argument("--summary", help="JSON summary output path")
    args = ap.parse_args()

    manifest = load_manifest(args.manifest)
    with open(args.accessions) as fh:
        accessions = [ln.strip() for ln in fh if ln.strip()]

    checksums = load_checksums(args.checksums) if args.checksums else {}

    os.makedirs(args.out, exist_ok=True)
    results = []
    n_ok = 0
    n_mismatch = 0
    n_missing = 0
    n_checksum_fail = 0

    for acc in accessions:
        if acc not in manifest:
            results.append({"acc": acc, "status": "NOT_IN_MANIFEST"})
            n_mismatch += 1
            continue
        exp = manifest[acc]
        fa = os.path.join(args.staging, acc, "contigs.fasta")
        if not os.path.exists(fa):
            results.append({"acc": acc, "status": "MISSING_FILE"})
            n_missing += 1
            continue

        # 1) checksum (transfer integrity) if a manifest was provided
        rel = f"{acc}/contigs.fasta"
        if checksums:
            want = checksums.get(rel)
            if want is None:
                results.append({"acc": acc, "status": "NO_CHECKSUM_ENTRY"})
                n_checksum_fail += 1
                continue
            got = sha256_file(fa)
            if got != want:
                results.append({"acc": acc, "status": "CHECKSUM_MISMATCH",
                                "detail": f"want {want} got {got}"})
                n_checksum_fail += 1
                continue

        # 2) content check: contig count + total length vs manifest
        n_contigs = 0
        total_bp = 0
        for token, seq in iter_records(fa):
            n_contigs += 1
            total_bp += len(seq)
        count_ok = n_contigs == exp["assembly_contigs"]
        bp_ok = total_bp == exp["assembly_length_bp"]
        if not (count_ok and bp_ok):
            results.append({
                "acc": acc, "status": "CONTENT_MISMATCH",
                "detail": (f"contigs got {n_contigs} want {exp['assembly_contigs']}; "
                           f"bp got {total_bp} want {exp['assembly_length_bp']}"),
            })
            n_mismatch += 1
            continue

        # 3) convert to PanSN bgzip
        out_dir = os.path.join(args.out, acc)
        os.makedirs(out_dir, exist_ok=True)
        out_stem = os.path.join(out_dir, f"{acc}.pansn")
        out_fa = out_stem + ".fa"
        with open(out_fa, "w") as out:
            for token, seq in iter_records(fa):
                out.write(f">{acc}#1#{token}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i + 80] + "\n")
        subprocess.run(["bgzip", "-f", out_fa], check=True)
        gz = out_fa + ".gz"  # bgzip -f replaces out_fa in place -> .gz
        subprocess.run(["samtools", "faidx", gz], check=True)
        results.append({"acc": acc, "status": "OK",
                        "contigs": n_contigs, "bp": total_bp})
        n_ok += 1

    # report
    report_lines = [
        "# Run-assembly ingestion report",
        "",
        f"- cohort: {len(accessions)} accessions",
        f"- OK (PanSN bgzip, content-verified): {n_ok}",
        f"- content/other mismatch: {n_mismatch}",
        f"- checksum failures: {n_checksum_fail}",
        f"- missing files (not yet uploaded): {n_missing}",
        "",
        "## Per-genome verification (contigs + total bp vs manifest)",
        "",
        "| accession | status | contigs | manifest_contigs | bp | manifest_bp | detail |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in sorted(results, key=lambda r: r["acc"]):
        report_lines.append(
            "| {acc} | {status} | {contigs} | {mcontigs} | {bp} | {mbp} | {detail} |".format(
                acc=r["acc"],
                status=r["status"],
                contigs=r.get("contigs", ""),
                mcontigs=manifest.get(r["acc"], {}).get("assembly_contigs", ""),
                bp=r.get("bp", ""),
                mbp=manifest.get(r["acc"], {}).get("assembly_length_bp", ""),
                detail=r.get("detail", ""),
            )
        )
    md = "\n".join(report_lines) + "\n"
    if args.report:
        os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
        with open(args.report, "w") as fh:
            fh.write(md)
    print(md)

    summary = {
        "cohort": len(accessions),
        "ok": n_ok,
        "mismatch": n_mismatch,
        "checksum_failures": n_checksum_fail,
        "missing": n_missing,
        "match_pct": round(100.0 * n_ok / len(accessions), 2) if accessions else 0.0,
    }
    if args.summary:
        with open(args.summary, "w") as fh:
            json.dump(summary, fh, indent=2)
    sys.exit(0 if (n_ok == len(accessions)) else 1)


if __name__ == "__main__":
    main()
