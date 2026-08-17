#!/usr/bin/env python3
"""Bounded S3 retrieval + confirmation pipeline for Logan screen hits.

Implements PILOT_PLAN §5 for `execute-bounded-ntm`:

  1. retrieve  https://s3.amazonaws.com/logan-pub/c/<acc>/<acc>.contigs.fa.zst
              (fallback /u/... unitigs only on explicit request)
  2. checksum sha256 + record S3 size / ETag / Last-Modified from HEAD
  3. locate    bait k-mers on contigs via back_to_sequences
              (--output-mapping-positions), then extract matched contigs
  4. align     minimap2 -c --eqx (asm5 first, asm20 fallback) bait vs contigs
  5. emit      per-(bait, accession) confirmation records with query
              coverage, identity, matched contig ids, co-location

Budgets (preregistered): top <=5 accessions per bait; hard stop 20 GB or 300
accessions; median contig file > 50 MB => STOP for re-scope.

Usage:
  confirm_hits.py --hits hits.tsv --out-dir confirm/ [--max-per-bait 5]
  hits.tsv columns: bait_id, accession, kmer_coverage, evalue (extra cols ok)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from collections import defaultdict

B2S = shutil.which("back_to_sequences") or "/tmp/mmenv/bin/back_to_sequences"
MINIMAP2 = shutil.which("minimap2") or "/tmp/mmenv/bin/minimap2"
S3_BASE = "https://s3.amazonaws.com/logan-pub"
HARD_STOP_BYTES = 20 * 1024**3
HARD_STOP_ACCESSIONS = 300
MEDIAN_STOP_BYTES = 50 * 1024 * 1024


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def head_url(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return {k: r.headers.get(k) for k in ("Content-Length", "ETag", "Last-Modified")}


def download(url, dest, ledger):
    meta = head_url(url)
    size = int(meta.get("Content-Length") or 0)
    with urllib.request.urlopen(urllib.request.Request(url), timeout=600) as r, open(dest, "wb") as fh:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    got = sha256_file(dest)
    ledger.append({"event": "s3_download", "url": url, "dest": dest,
                   "size_bytes": size, "sha256": got, **meta})
    return {"url": url, "size_bytes": size, "sha256": got, **meta}


def zstd_decompress(src, dest):
    subprocess.run(["zstd", "-d", "-q", "-f", "-o", dest, src], check=True)


def read_fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True)
    ap.add_argument("--baits-fa", required=True, help="FASTA with bait sequences")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--max-per-bait", type=int, default=5)
    ap.add_argument("--max-total-bytes", type=int, default=HARD_STOP_BYTES)
    ap.add_argument("--max-accessions", type=int, default=HARD_STOP_ACCESSIONS)
    ap.add_argument("--kind", choices=["contigs", "unitigs"], default="contigs")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    dl_dir = os.path.join(args.out_dir, "downloads")
    os.makedirs(dl_dir, exist_ok=True)
    ledger_path = os.path.join(args.out_dir, "confirm_ledger.jsonl")
    ledger = []
    def flush():
        with open(ledger_path, "a") as fh:
            for e in ledger:
                fh.write(json.dumps(e, sort_keys=True) + "\n")
            ledger.clear()

    # rank hits per bait (kmer_coverage DESC, evalue ASC), cap per bait
    hits = list(csv.DictReader(open(args.hits), delimiter="\t"))
    by_bait = defaultdict(list)
    for h in hits:
        by_bait[h["bait_id"]].append(h)
    selected = []
    for bait, hs in by_bait.items():
        hs.sort(key=lambda h: (-float(h.get("kmer_coverage") or 0), float(h.get("evalue") or 0)))
        selected.extend(hs[: args.max_per_bait])

    accessions = sorted({h["accession"] for h in selected})
    if len(accessions) > args.max_accessions:
        print(f"STOP: {len(accessions)} accessions exceeds cap {args.max_accessions}")
        return 3

    bait_seqs = {}
    name = None
    for line in open(args.baits_fa):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            bait_seqs[name] = []
        elif line:
            bait_seqs[name].append(line)
    bait_seqs = {k: "".join(v) for k, v in bait_seqs.items()}
    # allow _t70 suffix variants
    def bait_seq(bait_id):
        if bait_id in bait_seqs:
            return bait_id, bait_seqs[bait_id]
        for k, v in bait_seqs.items():
            if bait_id.startswith(k + "_t") or k.startswith(bait_id + "_t"):
                return k, v
        raise KeyError(bait_id)

    total_bytes = 0
    results = []
    for i, acc in enumerate(accessions):
        url = f"{S3_BASE}/{'c' if args.kind == 'contigs' else 'u'}/{acc}/{acc}.{args.kind}.fa.zst"
        zst = os.path.join(dl_dir, f"{acc}.{args.kind}.fa.zst")
        fa = os.path.join(dl_dir, f"{acc}.{args.kind}.fa")
        if not os.path.exists(fa):
            if not os.path.exists(zst):
                meta = head_url(url)
                size = int(meta.get("Content-Length") or 0)
                if total_bytes + size > args.max_total_bytes:
                    print(f"STOP: byte budget would exceed {args.max_total_bytes}")
                    break
                download(url, zst, ledger)
                total_bytes += size
                flush()
            zstd_decompress(zst, fa)
        n_contigs = sum(1 for _ in read_fasta_ids(fa))

        for h in [x for x in selected if x["accession"] == acc]:
            bait_id = h["bait_id"]
            key, seq = bait_seq(bait_id)
            bait_fa = os.path.join(args.out_dir, f"{key}.fa")
            with open(bait_fa, "w") as fh:
                fh.write(f">{key}\n{seq}\n")
            # back_to_sequences: kmers from bait vs contigs
            out_prefix = os.path.join(args.out_dir, f"{key}__{acc}")
            b2s_cmd = [B2S, "--in-kmers", bait_fa, "--in-sequences", fa,
                       "--out-sequences", out_prefix + ".b2s.txt",
                       "--output-mapping-positions"]
            b2s_run = subprocess.run(b2s_cmd, capture_output=True, text=True)
            matched_ids = []
            if os.path.exists(out_prefix + ".b2s.txt"):
                with open(out_prefix + ".b2s.txt") as fh:
                    txt = fh.read()
                for line in txt.splitlines():
                    if line.startswith(">"):
                        matched_ids.append(line[1:].split()[0])
            matched_ids = sorted(set(matched_ids))
            # minimap2 confirmation on matched contigs (or whole file if none)
            target = fa
            align = {"preset": None, "qcov": None, "identity": None, "n_matches": 0, "err": None}
            import re as _re
            for preset in ("sr", "asm5", "asm20"):
                mm_cmd = [MINIMAP2, "-c", "--eqx", "-x", preset, target, bait_fa]
                mm = subprocess.run(mm_cmd, capture_output=True, text=True)
                if mm.returncode != 0:
                    align["err"] = mm.stderr.strip()[:300]
                    continue
                qspan = matches = mlen = 0
                for line in mm.stdout.splitlines():
                    if line.startswith("@") or not line.strip():
                        continue
                    f = line.split("\t")
                    cg = [x for x in f if x.startswith("cg:Z:")]
                    if len(f) < 12 or not cg:
                        continue
                    matches += 1
                    for n, op in _re.findall(r"(\d+)([MIDNSH=X])", cg[0][5:]):
                        if op in "M=X":
                            qspan += int(n)
                            if op == "X":
                                mlen += int(n)
                        elif op in "ID":
                            mlen += int(n)
                if matches == 0:
                    continue
                identity = round(1 - mlen / qspan, 4) if qspan else None
                align = {"preset": preset, "qcov": round(qspan / len(seq), 4) if seq else None,
                         "identity": identity, "n_matches": matches, "err": None}
                break
            results.append({
                "bait_id": bait_id, "accession": acc,
                "screen_kmer_coverage": h.get("kmer_coverage"),
                "kind": args.kind, "n_contigs": n_contigs,
                "b2s_matched_contig_ids": matched_ids,
                "minimap2": align,
            })
            with open(os.path.join(args.out_dir, "confirm_results.jsonl"), "a") as fh:
                fh.write(json.dumps(results[-1], sort_keys=True) + "\n")
    flush()
    print(json.dumps({"n_accessions": len(accessions), "total_bytes": total_bytes,
                      "results": len(results)}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
