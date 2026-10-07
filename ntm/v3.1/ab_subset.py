#!/usr/bin/env python3
"""v3.1 A/B harness — same v3 clade + partitions, free vs observed reconstruction.

Isolates the walk change: clustering and alignment/partition inputs are IDENTICAL
to v3 (reused from $NVME/ntm/v3/clades/<cid>/). Only --path-mode varies.

Writes to $NVME/ntm/v3.1/ab/<cid>/{free,observed}/ and emits ab_summary.tsv.
"""
from __future__ import annotations
import argparse, csv, json, os, statistics, subprocess, sys

NV = "/mnt/nvme3n1/erikg/phind-genome-work"
REPO = "/home/erikg/phind"
CL = f"{NV}/ntm/v3/clades"
OUT = f"{NV}/ntm/v3.1/ab"


def read_fasta_lengths(path):
    lens, name, n = [], None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lens.append(n)
                name, n = line[1:].strip(), 0
            else:
                n += len(line.strip())
    if name is not None:
        lens.append(n)
    return lens


def read_single_fasta_len(path):
    if not os.path.exists(path):
        return None
    n = 0
    for line in open(path):
        if line.startswith(">"):
            continue
        n += len(line.strip())
    return n


def run(cid, mode, outdir_prefix, n_samples=25, target=None):
    cl = f"{CL}/{cid}"
    cmd = ["python3", f"{REPO}/scripts/traverse_partitions.py",
           "--partitions-dir", f"{cl}/partitions",
           "--bed", f"{cl}/partitions.bed",
           "--output", outdir_prefix,
           "--mode", "ml", "--n-samples", str(n_samples), "--seed", "42"]
    if mode == "observed":
        cmd += ["--path-mode", "observed"]
        if target:
            cmd += ["--target-length", str(int(target))]
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    return r.returncode, (r.stderr or "")[-400:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clades", required=True, help="file with one clade id per line")
    ap.add_argument("--modes", default="free,observed")
    ap.add_argument("--n-samples", type=int, default=25)
    ap.add_argument("--targets", default=None, help="TSV clade_id<TAB>target_len")
    args = ap.parse_args()
    targets = {}
    if args.targets:
        for r in csv.DictReader(open(args.targets), delimiter="\t"):
            try:
                targets[r["clade_id"]] = int(float(r["target_len"]))
            except Exception:
                pass
    cids = [c.strip() for c in open(args.clades) if c.strip()]
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for i, cid in enumerate(cids, 1):
        cl = f"{CL}/{cid}"
        mem = read_fasta_lengths(f"{cl}/sequences.fa")
        med = int(statistics.median(mem)) if mem else None
        rec = {"clade_id": cid, "n_members": len(mem), "median_member_len": med}
        for mode in args.modes.split(","):
            d = f"{OUT}/{cid}/{mode}"
            os.makedirs(d, exist_ok=True)
            prefix = f"{d}/{mode}"
            rc, err = run(cid, mode, prefix, args.n_samples, targets.get(cid))
            ml = read_single_fasta_len(f"{prefix}.ml.fa")
            st = {}
            sp = f"{prefix}.stats.json"
            if os.path.exists(sp):
                try:
                    st = json.load(open(sp))
                except Exception:
                    st = {}
            ps = st.get("path_stats") or {}
            rec[f"{mode}_rc"] = rc
            rec[f"{mode}_len"] = ml
            rec[f"{mode}_ratio"] = round(ml / med, 3) if (ml and med) else None
            if mode == "observed":
                rec["obs_terminated_by"] = ps.get("terminated_by")
                rec["obs_budget_hit"] = ps.get("budget_hit")
                rec["obs_best_by"] = ps.get("best_by")
                rec["obs_best_len"] = ps.get("best_len")
                rec["obs_median_member_len"] = ps.get("median_member_len")
                rec["obs_path_steps"] = ps.get("n_path_steps")
            if rc != 0:
                rec[f"{mode}_err"] = err
        rows.append(rec)
        print(f"[{i}/{len(cids)}] {cid} n={rec['n_members']} med={med} "
              f"free={rec['free_len']}({rec['free_ratio']}x) "
              f"obs={rec['observed_len']}({rec['observed_ratio']}x) "
              f"{rec.get('obs_terminated_by')}", flush=True)
    cols = ["clade_id", "n_members", "median_member_len", "free_len", "free_ratio",
            "observed_len", "observed_ratio", "obs_terminated_by", "obs_budget_hit",
            "obs_best_by", "obs_best_len", "obs_median_member_len", "obs_path_steps",
            "free_rc", "observed_rc"]
    with open(f"{OUT}/ab_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nwrote {OUT}/ab_summary.tsv")


if __name__ == "__main__":
    main()
