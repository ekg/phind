#!/usr/bin/env python3
"""v3.1 full ML-only re-run — genome-path reconstruction over the EXISTING v3 clades.

The mosaic defect is created by the ML step alone, so no re-alignment or
re-clustering is needed: this reuses $NVME/ntm/v3/clades/<cid>/{partitions,partitions.bed,sequences.fa}
and re-runs traverse_partitions.py with --path-mode observed.

- alignable clades (n>=2): traverse_partitions --mode ml --path-mode observed [--target-length T]
- singletons (n==1): pass through the member sequence as-is (v1/v2 convention)
- ancestral mode: optional (--ancestral) for the same clades

Targets (--targets TSV: clade_id<TAB>target_len) come from CheckV aai_expected_length
medians over a clade's members; clades absent from the TSV fall back to rule v2
(no --target-length passed).

Writes $NVME/ntm/v3.1/ml/<cid>/{ml.ml.fa,anc.ancestral.genome.fa,...} and a manifest.
"""
from __future__ import annotations
import argparse, csv, json, os, subprocess, sys
from concurrent.futures import ProcessPoolExecutor, as_completed

NV = "/mnt/nvme3n1/erikg/phind-genome-work"
REPO = "/home/erikg/phind"
CL = f"{NV}/ntm/v3/clades"
OUT = f"{NV}/ntm/v3.1/ml"


def read_clades(path):
    with open(path) as fh:
        return json.load(fh)


def read_fasta(path):
    recs, name, buf = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def do_clade(cid, members, target, ancestral, n_samples, seed):
    cl = f"{CL}/{cid}"
    od = f"{OUT}/{cid}"
    os.makedirs(od, exist_ok=True)
    res = {"clade_id": cid, "n_members": len(members), "target_len": target,
           "status": "", "ml_len": None, "anc_len": None, "err": ""}
    if len(members) == 1:
        m = members[0]
        recs = read_fasta(f"{cl}/sequences.fa")
        seq = {}
        for n, s in recs:
            seq[n.split()[0]] = s
        s = seq.get(m) or (list(seq.values())[0] if seq else None)
        if s is None:
            res["status"] = "singleton_no_seq"
            return res
        with open(f"{od}/ml.ml.fa", "w") as fh:
            fh.write(f">{cid}_ML status=singleton n_members=1 length={len(s)}\n{s}\n")
        res.update(status="singleton_passthrough", ml_len=len(s))
        return res
    if not os.path.exists(f"{cl}/partitions.bed"):
        res["status"] = "no_partitions"
        return res
    modes = [("ml", f"{od}/ml")]
    if ancestral:
        modes.append(("ancestral", f"{od}/anc"))
    for mode, prefix in modes:
        cmd = ["python3", f"{REPO}/scripts/traverse_partitions.py",
               "--partitions-dir", f"{cl}/partitions", "--bed", f"{cl}/partitions.bed",
               "--output", prefix, "--mode", mode,
               "--n-samples", str(n_samples), "--seed", str(seed),
               "--path-mode", "observed"]
        if target:
            cmd += ["--target-length", str(int(target))]
        p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
        if p.returncode != 0:
            res["status"] = "failed"
            res["err"] = (p.stderr or "")[-300:]
            return res
        fa = f"{prefix}.ml.fa" if mode == "ml" else f"{prefix}.ancestral.genome.fa"
        n = 0
        if os.path.exists(fa):
            for line in open(fa):
                if not line.startswith(">"):
                    n += len(line.strip())
        if mode == "ml":
            res["ml_len"] = n
        else:
            res["anc_len"] = n
    res["status"] = "ok"
    return res


def main():
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--clades", default="/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/clades/0/tight_clades.json")
    ap.add_argument("--targets", default=None, help="TSV clade_id<TAB>target_len")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--jobs", type=int, default=24)
    ap.add_argument("--n-samples", type=int, default=25)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--ancestral", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--clade", action="append")
    args = ap.parse_args()
    OUT = args.out
    os.makedirs(OUT, exist_ok=True)

    clades = read_clades(args.clades)
    targets = {}
    if args.targets:
        for r in csv.DictReader(open(args.targets), delimiter="\t"):
            try:
                targets[r["clade_id"]] = int(float(r["target_len"]))
            except Exception:
                pass
    cids = args.clade or sorted(clades)
    if args.limit:
        cids = cids[: args.limit]

    rows = []
    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {ex.submit(do_clade, c, clades[c], targets.get(c), args.ancestral,
                          args.n_samples, args.seed): c for c in cids}
        done = 0
        for f in as_completed(futs):
            done += 1
            try:
                rows.append(f.result())
            except Exception as e:
                rows.append({"clade_id": futs[f], "status": "exception", "err": str(e)[:200]})
            if done % 100 == 0:
                print(f"  {done}/{len(cids)}", flush=True)
    cols = ["clade_id", "n_members", "target_len", "status", "ml_len", "anc_len", "err"]
    with open(f"{OUT}/v31_ml_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["clade_id"]):
            w.writerow(r)
    ok = sum(1 for r in rows if r["status"] in ("ok", "singleton_passthrough"))
    print(f"V31_ML_DONE ok={ok} total={len(rows)} -> {OUT}/v31_ml_summary.tsv")


if __name__ == "__main__":
    main()
