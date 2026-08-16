#!/usr/bin/env python3
"""ENA backfill feasibility: download ENA submitted assemblies for the
probed hit accessions, content-check against the v2 manifest, and ingest
into the PanSN bgzip layout with faidx validation.

Pipeline per accession:
  1. download https://{submitted_ftp} -> {staging}/{acc}/{acc}.ena.fa.gz
     (resumable: existing file is reused; gzip integrity re-checked)
  2. content check: contig count + total bp vs manifest
     (assembly_contigs / assembly_length_bp)
  3. PanSN conversion (independent of the content check result):
     {out}/{acc}/{acc}.pansn.fa.gz with headers {acc}#1#{source_token}
     (token preserved byte-for-byte from the ENA FASTA)
  4. samtools faidx validates the bgzip; .fai line count must equal the
     ENA contig count
  5. CONTENT_MISMATCH accessions are quarantined with a marker file
     (QUARANTINE.txt) in the out dir — they are NOT eligible to substitute
     for the collaborator file without re-verification.

A per-genome agreement table (ENA vs manifest) + machine summary JSON are
written. Download-only mode (--download-only) skips PanSN conversion.

Usage:
  python ntm/v2/scripts/ena_backfill_feasibility.py \
      --probe ntm/v2/run_assemblies/ena_probe_results.tsv \
      --manifest /mnt/.../NTM_QC_passed_prophage_master_manifest.tsv \
      --staging /mnt/.../ntm/v2/incoming/ena_backfill \
      --out /mnt/.../ntm/v2/genomes/ena_backfill \
      --max 50 \
      [--prefer-prophage] \
      [--download-only] \
      [--report ntm/v2/run_assemblies/ena_backfill_report.md] \
      [--summary ntm/v2/run_assemblies/ena_backfill_summary.json]
"""
import argparse
import csv
import gzip
import json
import os
import re
import subprocess
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

UA = "Mozilla/5.0 (phind NTM v2 ENA backfill feasibility; erik@hypervolu.me)"
PANS_PREFIX_RE = re.compile(r"^>([^\s]+)")


def iter_records(fasta_path):
    """Yield (source_token, seq) from a plain or gzip FASTA.

    gzip is detected by magic bytes (\x1f\x8b) so temp/cached files with
    non-.gz names are handled correctly.
    """
    with open(fasta_path, "rb") as fh:
        magic = fh.read(2)
    opener = gzip.open if magic == b"\x1f\x8b" else open
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


def check_gzip(path):
    """Verify gzip integrity + count records; returns (n_contigs, total_bp)."""
    n = 0
    bp = 0
    for _token, seq in iter_records(path):
        n += 1
        bp += len(seq)
    return n, bp


def download_one(acc, ftp_rel, staging, workers_ok):
    """Download one assembly; returns (acc, ok, detail)."""
    dest = os.path.join(staging, acc, f"{acc}.ena.fa.gz")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        try:
            n, bp = check_gzip(dest)
            return acc, True, f"cached (n={n} bp={bp})"
        except Exception as e:  # corrupt cache -> redownload
            os.unlink(dest)
    url = "https://" + ftp_rel
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    tmp = dest + ".part"
    try:
        with urllib.request.urlopen(req, timeout=180) as resp, open(tmp, "wb") as fh:
            while True:
                b = resp.read(1 << 20)
                if not b:
                    break
                fh.write(b)
        # integrity check before promoting
        n, bp = check_gzip(tmp)
        os.replace(tmp, dest)
        return acc, True, f"downloaded (n={n} bp={bp})"
    except Exception as e:
        if os.path.exists(tmp):
            os.unlink(tmp)
        return acc, False, f"{type(e).__name__}: {e}"


def load_manifest(manifest_path):
    out = {}
    with open(manifest_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            acc = r["accession"]
            if acc[:3] not in ("ERR", "SRR", "DRR"):
                continue
            out[acc] = {
                "assembly_contigs": int(r["assembly_contigs"]),
                "assembly_length_bp": int(r["assembly_length_bp"]),
                "has_prophage": r["has_prophage"] == "True",
            }
    return out


def to_pansn_bgzip(acc, src_gz, out_dir):
    """Convert ENA fasta.gz -> PanSN bgzip {acc}.pansn.fa.gz + .fai + .gzi."""
    os.makedirs(out_dir, exist_ok=True)
    out_stem = os.path.join(out_dir, f"{acc}.pansn")
    out_fa = out_stem + ".fa"
    n = 0
    with open(out_fa, "w") as out:
        for token, seq in iter_records(src_gz):
            n += 1
            out.write(f">{acc}#1#{token}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
    subprocess.run(["bgzip", "-f", out_fa], check=True)
    gz = out_fa + ".gz"
    subprocess.run(["samtools", "faidx", gz], check=True)
    # faidx sanity: .fai line count == contig count
    with open(gz + ".fai") as fh:
        fai_lines = sum(1 for _ in fh)
    return gz, n, fai_lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", required=True, help="ena_probe_results.tsv")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--staging", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max", type=int, default=50, help="batch size")
    ap.add_argument("--prefer-prophage", action="store_true",
                    help="prophage-bearing accessions first (priority order)")
    ap.add_argument("--download-only", action="store_true",
                    help="skip PanSN conversion (just download + content check)")
    ap.add_argument("--report", help="markdown report path")
    ap.add_argument("--summary", help="JSON summary path")
    args = ap.parse_args()

    with open(args.probe) as fh:
        probe = [r for r in csv.DictReader(fh, delimiter="\t") if r["submitted_ftp"]]

    hits = []
    for r in probe:
        ftp_rel = re.sub(r"^[a-z]+://", "", r["submitted_ftp"])
        hits.append({
            "run": r["run"],
            "sample_accession": r["sample_accession"],
            "ftp_rel": ftp_rel,
            "has_prophage": r["has_prophage"] == "True",
            "manifest_contigs": int(r["assembly_contigs"]),
            "manifest_bp": int(r["assembly_length_bp"]),
        })

    if args.prefer_prophage:
        hits.sort(key=lambda h: (not h["has_prophage"], h["run"]))
    else:
        hits.sort(key=lambda h: h["run"])
    batch = hits[: args.max]
    sys.stderr.write(f"batch: {len(batch)} accessions "
                     f"({sum(1 for h in batch if h['has_prophage'])} prophage-bearing)\n")

    manifest = load_manifest(args.manifest)

    # --- download ---
    results = {}
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(download_one, h["run"], h["ftp_rel"], args.staging, None):
                h for h in batch}
        done = 0
        for fut in as_completed(futs):
            acc, ok, detail = fut.result()
            h = futs[fut]
            rec = dict(h)
            rec["download_ok"] = ok
            rec["download_detail"] = detail
            results[acc] = rec
            done += 1
            sys.stderr.write(f"\rdownloaded {done}/{len(batch)} {acc} ok={ok}  ")
            sys.stderr.flush()
    sys.stderr.write("\n")

    # --- content check + PanSN ---
    os.makedirs(args.out, exist_ok=True)
    n_match = 0
    n_mismatch = 0
    for h in batch:
        acc = h["run"]
        rec = results[acc]
        src = os.path.join(args.staging, acc, f"{acc}.ena.fa.gz")
        if not rec["download_ok"] or not os.path.exists(src):
            rec["status"] = "DOWNLOAD_FAILED"
            rec["detail"] = rec["download_detail"]
            continue
        n_ena, bp_ena = check_gzip(src)
        rec["ena_contigs"] = n_ena
        rec["ena_bp"] = bp_ena
        match = (n_ena == rec["manifest_contigs"] and bp_ena == rec["manifest_bp"])
        rec["match"] = match
        if match:
            n_match += 1
        else:
            n_mismatch += 1
        if args.download_only:
            rec["status"] = "MATCH" if match else "CONTENT_MISMATCH"
            rec["detail"] = "" if match else (
                f"contigs got {n_ena} want {rec['manifest_contigs']}; "
                f"bp got {bp_ena} want {rec['manifest_bp']}")
            continue
        # PanSN conversion (all downloaded assemblies, match or not)
        out_dir = os.path.join(args.out, acc)
        gz, n_pan, fai_lines = to_pansn_bgzip(acc, src, out_dir)
        rec["pansn_gz"] = gz
        rec["fai_lines"] = fai_lines
        fai_ok = fai_lines == n_ena
        if not match:
            # quarantine marker
            with open(os.path.join(out_dir, "QUARANTINE.txt"), "w") as fh:
                fh.write(
                    f"ENA assembly differs from manifest\n"
                    f"ena_contigs={n_ena} manifest_contigs={rec['manifest_contigs']}\n"
                    f"ena_bp={bp_ena} manifest_bp={rec['manifest_bp']}\n"
                    f"reason: manifest prophage coordinates computed on collaborator "
                    f"contigs; ENA contig set is different\n")
        rec["status"] = "OK_INGESTED" if (match and fai_ok) else (
            "INGESTED_QUARANTINED" if fai_ok else "FAIDX_FAIL")
        rec["detail"] = "" if match else (
            f"contigs got {n_ena} want {rec['manifest_contigs']}; "
            f"bp got {bp_ena} want {rec['manifest_bp']}")

    # --- report ---
    lines = [
        "# ENA backfill feasibility — batch report",
        "",
        f"- probe input: `{args.probe}` (hits: {len(hits)}, batch: {len(batch)})",
        f"- downloads OK: {sum(1 for r in results.values() if r['download_ok'])}/{len(batch)}",
        f"- content MATCH vs manifest: {n_match}",
        f"- content MISMATCH vs manifest: {n_mismatch}",
        f"- PanSN bgzip ingested (faidx validated): "
        f"{sum(1 for r in results.values() if r.get('status','').startswith('OK_INGESTED')) + sum(1 for r in results.values() if r.get('status','').startswith('INGESTED_QUARANTINED'))}",
        "",
        "## Per-genome ENA-vs-manifest agreement",
        "",
        "| accession | sample | prophage | status | ena_contigs | manifest_contigs | ena_bp | manifest_bp | match | detail |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for acc in sorted(results):
        r = results[acc]
        lines.append(
            "| {run} | {sample} | {prophage} | {status} | {ena_c} | {man_c} | {ena_bp} | {man_bp} | {match} | {detail} |".format(
                run=r["run"],
                sample=r["sample_accession"],
                prophage="Y" if r["has_prophage"] else "N",
                status=r.get("status", ""),
                ena_c=r.get("ena_contigs", ""),
                man_c=r["manifest_contigs"],
                ena_bp=r.get("ena_bp", ""),
                man_bp=r["manifest_bp"],
                match=("MATCH" if r.get("match") else "MISMATCH") if "match" in r else "",
                detail=r.get("detail", ""),
            )
        )
    md = "\n".join(lines) + "\n"
    if args.report:
        os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
        with open(args.report, "w") as fh:
            fh.write(md)
    print(md)

    summary = {
        "batch": len(batch),
        "downloads_ok": sum(1 for r in results.values() if r["download_ok"]),
        "downloads_failed": sum(1 for r in results.values() if not r["download_ok"]),
        "content_match": n_match,
        "content_mismatch": n_mismatch,
        "prophage_in_batch": sum(1 for r in results.values() if r["has_prophage"]),
        "prophage_mismatch": sum(1 for r in results.values()
                                 if r["has_prophage"] and "match" in r and not r["match"]),
    }
    if args.summary:
        with open(args.summary, "w") as fh:
            json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
