#!/usr/bin/env python3
"""
resume_per_clade_partition.py — drive the per-clade alignment pipeline to
completion, resuming any partial run and recording explicit failures.

Given the tight-clades file and pipeline inputs, this script:

  1. Determines which clades are incomplete (no `manifest.json` with a
     finished `pipeline.total_runtime_s` in <outdir>/<clade_id>/).
  2. Re-invokes `scripts/per_clade_alignment_pipeline.py` for exactly those
     clades (its built-in resume makes re-runs idempotent and cheap).
  3. If the pipeline exits non-zero (a worker crashed / clade error), retries
     the still-incomplete set — completed clades are never re-run.
  4. When the pipeline completes cleanly (exit 0) but some clades are STILL
     incomplete, those are genuine per-clade failures: after `--max-retries`
     clean runs they get an explicit `<outdir>/<clade_id>/FAILED.log` (last
     command + stderr from commands.log) and are dropped from retries.

This mirrors the v1 practice of capping concurrent allwave processes via the
pipeline's `--jobs` (each worker runs one clade end-to-end: extract -> allwave
-> segment -> partition).

Usage:
  python3 ntm/scripts/resume_per_clade_partition.py \
      --outdir /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/clades \
      --clades /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/clades/0/tight_clades.json \
      --fasta /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/full_prophages.fa \
      --index /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/mash_clades/full_prophages.idx.json \
      --threads 8 --jobs 16 \
      [--max-attempts 5] [--dry-run] [--clade 0_0000 ...]
"""
import argparse
import json
import os
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PIPELINE = os.path.join(REPO, "scripts", "per_clade_alignment_pipeline.py")


def is_complete(outdir, cid):
    mpath = os.path.join(outdir, cid, "manifest.json")
    if not os.path.exists(mpath):
        return False
    try:
        m = json.load(open(mpath))
    except Exception:
        return False
    return m.get("pipeline", {}).get("total_runtime_s") is not None


def is_failed(outdir, cid):
    return os.path.exists(os.path.join(outdir, cid, "FAILED.log"))


def incomplete_clades(outdir, clades):
    inc = []
    for cid in clades:
        if is_complete(outdir, cid) or is_failed(outdir, cid):
            continue
        inc.append(cid)
    return inc


def failure_reason(outdir, cid):
    """Best-effort reason: tail of the clade's commands.log."""
    log = os.path.join(outdir, cid, "commands.log")
    if os.path.exists(log):
        with open(log) as f:
            lines = f.read().strip().splitlines()
        return "\n".join(lines[-8:])[-1500:]
    return "no commands.log; no output produced"


def mark_failed(outdir, cid, attempt):
    d = os.path.join(outdir, cid)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "FAILED.log"), "w") as f:
        f.write(f"clade {cid}: failed after {attempt} clean pipeline runs\n")
        f.write("reason (commands.log tail):\n")
        f.write(failure_reason(outdir, cid) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--clades", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--max-attempts", type=int, default=5,
                    help="max total pipeline invocations (default 5)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--clade", action="append", help="restrict to these clades")
    args = ap.parse_args()

    tc = json.load(open(args.clades))
    cids = args.clade or list(tc.keys())
    os.makedirs(args.outdir, exist_ok=True)

    attempt = 0
    t0 = time.time()
    while True:
        inc = incomplete_clades(args.outdir, cids)
        if not inc:
            print(f"[resume] all {len(cids)} clades complete or failed "
                  f"({time.time()-t0:.0f}s total)", flush=True)
            return 0
        attempt += 1
        print(f"[resume] attempt {attempt}: {len(inc)} incomplete clades "
              f"(of {len(cids)})", flush=True)
        if args.dry_run:
            print("  would run:", " ".join(inc[:10]),
                  ("..." if len(inc) > 10 else ""), flush=True)
            return 0
        cmd = [sys.executable, PIPELINE,
               "--community", "0",
               "--outdir", args.outdir,
               "--clades", args.clades,
               "--fasta", args.fasta,
               "--index", args.index,
               "--threads", str(args.threads),
               "--jobs", str(args.jobs)] + \
              [x for c in inc for x in ("--clade", c)]
        print("+ " + " ".join(cmd[:12]) + f" ... ({len(inc)} clades)",
              flush=True)
        r = subprocess.run(cmd)
        if r.returncode != 0:
            print(f"[resume] pipeline exited {r.returncode}; will retry "
                  f"remaining", flush=True)
            if attempt >= args.max_attempts:
                print(f"[resume] attempt budget exhausted ({args.max_attempts})"
                      f" with pipeline crashes; stopping", flush=True)
                return 1
            continue  # crash, not a per-clade failure: no FAILED.log yet
        # clean exit: whatever remains is a genuine per-clade failure
        still = incomplete_clades(args.outdir, cids)
        if not still:
            print(f"[resume] all clades done after attempt {attempt} "
                  f"({time.time()-t0:.0f}s)", flush=True)
            return 0
        for cid in still:
            mark_failed(args.outdir, cid, attempt)
        print(f"[resume] {len(still)} clades marked FAILED after clean run "
              f"(attempt {attempt}):", flush=True)
        for cid in still:
            print(f"  {cid}: {failure_reason(args.outdir, cid)[:200]}",
                  flush=True)
        if attempt >= args.max_attempts:
            print(f"[resume] attempt budget exhausted ({args.max_attempts})",
                  flush=True)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
