#!/usr/bin/env python3
"""Confirmation of tier screen hits (RUN_PLAN §8).

Per (bait, accession) with kmer_coverage >= threshold:
  1. b2s --out-sequences + --output-mapping-positions -> matched contig ids
     + covered positions on each contig
  2. minimap2 -c --eqx preset ladder sr -> asm5 -> asm20 vs matched contigs
  3. per-accession co-location: baits x contigs matrix

Classification per PILOT_PLAN §5 (frozen):
  - accession-level screen hit: cov >= 0.7
  - module-only homology: >= 1 interior-module bait aligned (qcov >= 0.5)
  - junction/synteny supported: junction bait aligned >= 0.5 qcov on ONE
    contig (a full-length junction-bait alignment IS same-contig support by
    construction: the bait spans the host/phage partition junction)
  - near-complete analogue: >= 2 module baits + >= 1 junction bait, all
    aligned on the same contig, combined bait-coverage span >= 20 kb
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

S3 = "https://s3.amazonaws.com/logan-pub/c/{acc}/{acc}.contigs.fa.zst"


def load_bait_seqs(panel_fa):
    baits, name, buf = {}, None, []
    for line in open(panel_fa):
        if line.startswith(">"):
            if name:
                baits[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        baits[name] = "".join(buf)
    return baits


def contig_len_map(fa):
    lens, name, n = {}, None, 0
    for line in open(fa):
        if line.startswith(">"):
            if name:
                lens[name] = n
            name, n = line[1:].split()[0], 0
        else:
            n += len(line.strip())
    if name:
        lens[name] = n
    return lens


def parse_pos_positions(txt):
    """Parse b2s --out-sequences+mapping-positions header:
    >contig ka:f:x  nshared ratio p1 p2 ... (span)"""
    out = {}
    for line in open(txt):
        if not line.startswith(">"):
            continue
        parts = line[1:].split()
        cid = parts[0]
        nums = [p for p in parts[1:] if not p.startswith("ka:") and ":" not in p]
        try:
            nshared = int(nums[0])
        except (ValueError, IndexError):
            continue
        span = None
        if nums and nums[-1].startswith("(") and nums[-1].endswith(")"):
            try:
                span = int(nums[-1][1:-1])
            except ValueError:
                pass
        out[cid] = {"n_shared_kmers": nshared, "covered_positions": None,
                    "span": span}
    return out


def minimap2_ladder(minimap2, target_fa, query_fa, ladder=("sr", "asm5", "asm20")):
    """PAF: col0=query(bait) name, 2-3=query span, 4=strand, 5=target(contig)
    name, 7-8=target span. qcov = UNION of query intervals / query length."""
    for preset in ladder:
        r = subprocess.run([minimap2, "-c", "--eqx", "-x", preset, target_fa, query_fa],
                           capture_output=True, text=True)
        if r.returncode != 0:
            continue
        alns = []
        for line in r.stdout.splitlines():
            if line.startswith("@") or not line.strip():
                continue
            f = line.split("\t")
            cg = [x for x in f if x.startswith("cg:Z:")]
            if len(f) < 12 or not cg:
                continue
            m, mm = 0, 0
            for n, op in re.findall(r"(\d+)([MIDNSH=X])", cg[0][5:]):
                if op in "M=X":
                    m += int(n)
                    if op == "X":
                        mm += int(n)
                elif op in "ID":
                    mm += int(n)
            alns.append({"query": f[0], "q_start": int(f[2]), "q_end": int(f[3]),
                         "strand": f[4], "target": f[5],
                         "t_start": int(f[7]), "t_end": int(f[8]),
                         "aln_len": m, "identity": round(1 - mm / m, 4) if m else None,
                         "cigar": cg[0][5:]})
        if alns:
            best = max(alns, key=lambda a: a["aln_len"])
            # union of query intervals
            ivs = sorted((a["q_start"], a["q_end"]) for a in alns)
            union = []
            for s, e in ivs:
                if union and s <= union[-1][1]:
                    union[-1][1] = max(union[-1][1], e)
                else:
                    union.append([s, e])
            qcov_span = sum(e - s for s, e in union)
            return {"preset": preset, "n_alignments": len(alns),
                    "qcov_union_span": qcov_span, "n_alns": alns,
                    "best": {"target": best["target"],
                             "t_start": best["t_start"], "t_end": best["t_end"],
                             "identity": best["identity"],
                             "aln_len": best["aln_len"], "strand": best["strand"]}}
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="tier results jsonl")
    ap.add_argument("--baits", required=True)
    ap.add_argument("--downloads", required=True)
    ap.add_argument("--b2s", required=True)
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threshold", type=float, default=0.7)
    args = ap.parse_args()

    baits = load_bait_seqs(args.baits)
    rows = [json.loads(l) for l in open(args.results)]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    scratch = os.path.join(os.path.dirname(args.out), "scratch")
    os.makedirs(scratch, exist_ok=True)

    out_fh = open(args.out, "a")
    n_conf = 0
    for r in rows:
        acc = r["accession"]
        cov = r["per_bait_kmer_cov"]
        hit_baits = [b for b, v in cov.items()
                     if v >= args.threshold and b in baits and not b.startswith("NTMBAIT_SHUF")]
        if not hit_baits:
            continue
        zst = os.path.join(args.downloads, f"{acc}.contigs.fa.zst")
        fa = zst[:-len(".zst")]
        if not os.path.exists(fa):
            subprocess.run(["zstd", "-d", "-q", "-f", "-o", fa, zst], check=True)
        lens = contig_len_map(fa)
        # localise each hit bait independently
        bait_contig = {}
        for b in hit_baits:
            bait_fa = os.path.join(scratch, f"{b}.fa")
            with open(bait_fa, "w") as fh:
                fh.write(f">{b}\n{baits[b]}\n")
            seqs_out = os.path.join(scratch, f"{b}__{acc}.matched.fa")
            cmd = [args.b2s, "--in-kmers", bait_fa, "--in-sequences", fa,
                   "--out-sequences", seqs_out, "--output-mapping-positions"]
            run = subprocess.run(cmd, capture_output=True, text=True)
            matched = parse_pos_positions(seqs_out) if os.path.exists(seqs_out) else {}
            os.remove(seqs_out)
            aln = minimap2_ladder(args.minimap2, fa, bait_fa)
            qlen = len(baits[b])
            bait_contig[b] = {
                "kmer_cov": cov[b],
                "matched_contigs": {c: m for c, m in matched.items()},
                "minimap2": aln,
                "qcov": round(aln["qcov_union_span"] / qlen, 4) if aln else 0.0,
                "identity": aln["best"]["identity"] if aln else None,
                "best_contig": aln["best"]["target"] if aln else None,
                "contig_len": lens.get(aln["best"]["target"]) if aln else None,
            }
        # classification (frozen ladder)
        modules = [b for b in hit_baits if "INTERIOR" in b]
        junctions = [b for b in hit_baits if "JUNCTION" in b]
        aligned_mods = [b for b in modules if bait_contig[b]["qcov"] >= 0.5]
        aligned_juncs = [b for b in junctions if bait_contig[b]["qcov"] >= 0.5]
        per_contig = {}
        for b, info in bait_contig.items():
            if info["best_contig"]:
                per_contig.setdefault(info["best_contig"], []).append(b)
        same_contig_pairs = {c: sorted(bs) for c, bs in per_contig.items()
                             if len(bs) >= 2}
        cat = "accession-level screen hit"
        if aligned_mods or aligned_juncs:
            cat = "module-only homology"
        if aligned_juncs:
            cat = "junction/synteny supported"
        if aligned_juncs and len(aligned_mods) >= 2:
            contigs_with_both = [c for c, bs in per_contig.items()
                                 if any(x in aligned_mods for x in bs)
                                 and any(x in aligned_juncs for x in bs)]
            span = sum(bait_contig[b]["qcov"] * len(baits[b]) for b in hit_baits)
            if contigs_with_both and span >= 20000:
                cat = "near-complete analogue candidate"
        rec = {"accession": acc, "tier": r["tier"],
               "scientific_name": r.get("scientific_name"),
               "n_baits_ge_thr": len(hit_baits), "classification": cat,
               "baits": bait_contig, "same_contig_baits": same_contig_pairs}
        out_fh.write(json.dumps(rec, sort_keys=True) + "\n")
        out_fh.flush()
        n_conf += 1
        print(f"CONFIRMED {acc} [{cat}] {len(hit_baits)} baits", flush=True)
    out_fh.close()
    print(f"{n_conf} accessions confirmed -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
