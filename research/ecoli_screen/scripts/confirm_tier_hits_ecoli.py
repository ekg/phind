#!/usr/bin/env python3
"""Confirmation of tier screen hits (RUN_PLAN §8) — E. coli generalization.

Generalized from ntm local_screen/confirm_tier_hits.py (validated on 97/97
NTM hits); diffs: ECBAIT tokens, control exclusions include HOST controls,
near-complete-analogue rule = task definition (>= 2 non-control baits >= 0.7
in one accession, each alignment-confirmed qcov >= 0.5); the NTM
stricter rule (>=1 junction + >=2 modules same contig, span >= 20 kb) is
reported alongside as a secondary, non-promoting statistic.

Per (bait, accession) with kmer_coverage >= threshold:
  1. b2s --out-sequences + --output-mapping-positions -> matched contig ids
  2. minimap2 -c --eqx preset ladder sr -> asm5 -> asm20 vs matched contigs
  3. per-accession co-location: baits x contigs matrix
Classification ladder (frozen): accession-level screen hit -> module-only
homology -> junction/synteny supported (full-length alignment on ONE contig;
a full-length junction-bait alignment IS same-contig support by construction)
-> near-complete analogue (>= 2 non-control baits aligned).
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys

S3 = "https://s3.amazonaws.com/logan-pub/c/{acc}/{acc}.contigs.fa.zst"
CONTROL_PREFIX = ("ECBAIT_SHUF", "ECBAIT_HOST")


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
                             "aln_len": best["aln_len"], "strand": best["strand"],
                             "cigar": best["cigar"]}}
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
    ap.add_argument("--max-pairs", type=int, default=40000,
                    help="RUN_PLAN §8 volume guard: max (bait, accession) pairs")
    ap.add_argument("--seed", type=int, default=20260822)
    ap.add_argument("--availability", help="availability.tsv for sha256 provenance")
    ap.add_argument("--runs-tsv", nargs="*",
                    help="frozen tier TSVs to join bioproject per run")
    args = ap.parse_args()

    proj_by_run = {}
    for p in (args.runs_tsv or []):
        import csv as _csv0
        for r in _csv0.DictReader(open(p), delimiter="\t"):
            proj_by_run.setdefault(r["run"], r.get("bioproject", ""))

    baits = load_bait_seqs(args.baits)
    sha_by_acc = {}
    if args.availability and os.path.exists(args.availability):
        import csv as _csv
        for r in _csv.DictReader(open(args.availability), delimiter="\t"):
            if r.get("sha256") and r.get("status") == "200":
                sha_by_acc.setdefault(r["accession"], r["sha256"])

    rows = [json.loads(l) for l in open(args.results)]
    pairs = []
    for r in rows:
        acc = r["accession"]
        for b, v in r["per_bait_kmer_cov"].items():
            if v >= args.threshold and b in baits \
                    and not b.startswith(CONTROL_PREFIX) \
                    and b not in ("ECBAIT_POSCON_LAMBDA_01", "ECBAIT_POSCON_HK97_02"):
                pairs.append((acc, b))
    pairs.sort()
    truncated = False
    if len(pairs) > args.max_pairs:
        import random as _random
        n_before = len(pairs)
        rng = _random.Random(args.seed)
        keep = set(rng.sample(range(n_before), args.max_pairs))
        pairs = [p for i, p in enumerate(pairs) if i in keep]
        truncated = True
        print(f"WARNING: volume guard engaged, sampled {args.max_pairs}/{n_before} pairs "
              f"(seed {args.seed}) — MUST be logged as an amendment")

    rows_by_acc = {r["accession"]: r for r in rows}
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    scratch = os.path.join(os.path.dirname(args.out), "scratch")
    os.makedirs(scratch, exist_ok=True)

    out_fh = open(args.out, "a")
    n_conf = 0
    accs = sorted({a for a, _ in pairs})
    for acc in accs:
        r = rows_by_acc[acc]
        cov = r["per_bait_kmer_cov"]
        hit_baits = [b for b, v in cov.items()
                     if v >= args.threshold and b in baits
                     and not b.startswith(CONTROL_PREFIX)
                     and b not in ("ECBAIT_POSCON_LAMBDA_01", "ECBAIT_POSCON_HK97_02")]
        if not hit_baits:
            continue
        zst = os.path.join(args.downloads, f"{acc}.contigs.fa.zst")
        fa = zst[:-len(".zst")]
        if not os.path.exists(fa):
            subprocess.run(["zstd", "-d", "-q", "-f", "-o", fa, zst], check=True)
        lens = contig_len_map(fa)
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
                "target_span": (f'{aln["best"]["t_start"]}-{aln["best"]["t_end"]}'
                                if aln else None),
                "cigar": (aln["best"]["cigar"][:200] if aln else None),
            }
        # classification (frozen ladder; RUN_PLAN §8 task-definition rule)
        modules = [b for b in hit_baits if "INTERIOR" in b or "MEMBER" in b]
        junctions = [b for b in hit_baits if "JUNCTION" in b]
        aligned_mods = [b for b in modules if bait_contig[b]["qcov"] >= 0.5]
        aligned_juncs = [b for b in junctions if bait_contig[b]["qcov"] >= 0.5]
        aligned_all = [b for b in hit_baits if bait_contig[b]["qcov"] >= 0.5]
        per_contig = {}
        for b, info in bait_contig.items():
            if info["best_contig"]:
                per_contig.setdefault(info["best_contig"], []).append(b)
        same_contig_baits = {c: sorted(bs) for c, bs in per_contig.items()
                             if len(bs) >= 2}
        cat = "accession-level screen hit"
        if aligned_mods or aligned_juncs:
            cat = "module-only homology"
        if aligned_juncs:
            cat = "junction/synteny supported"
        if len(aligned_all) >= 2:
            cat = "near-complete analogue"
        # NTM secondary (non-promoting): >=1 junction + >=2 modules, same
        # contig has both, combined aligned span >= 20 kb
        ntm_rule = False
        if aligned_juncs and len(aligned_mods) >= 2:
            contigs_with_both = [c for c, bs in per_contig.items()
                                 if any(x in aligned_mods for x in bs)
                                 and any(x in aligned_juncs for x in bs)]
            span = sum(bait_contig[b]["qcov"] * len(baits[b]) for b in hit_baits)
            ntm_rule = bool(contigs_with_both and span >= 20000)
        rec = {"accession": acc, "tier": r.get("tier"),
               "scientific_name": r.get("scientific_name"),
               "bioproject": r.get("bioproject") or proj_by_run.get(acc, ""),
               "s3_url": S3.format(acc=acc),
               "zst_sha256": sha_by_acc.get(acc),
               "n_baits_ge_thr": len(hit_baits), "classification": cat,
               "ntm_strict_rule_near_complete": ntm_rule,
               "baits": bait_contig, "same_contig_baits": same_contig_baits}
        out_fh.write(json.dumps(rec, sort_keys=True) + "\n")
        out_fh.flush()
        n_conf += 1
        print(f"CONFIRMED {acc} [{cat}] {len(hit_baits)} baits"
              + (" [ntm-strict too]" if ntm_rule else ""), flush=True)
    out_fh.close()
    print(f"{n_conf} accessions confirmed -> {args.out}"
          + (" [TRUNCATED]" if truncated else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
