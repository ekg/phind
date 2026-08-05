#!/usr/bin/env python3
"""
per_clade_alignment_pipeline.py — extract → allwave → segment → partition
for tight clades (research/clades/<community>/<clade_id>/).

Per clade (task per-clade-fastga):
  1. Extract sequences.fa from the full prophage FASTA (offset-indexed:
     build_research/clades/full_prophages.idx.json once with --build-index).
  2. Align with allwave:
       - n == 1        : skipped (no pairs; strategy "none (singleton)")
       - 2 <= n <= 50  : -p none (all-pairs; small clades)
       - n > 50        : -p tree:5:0:0.05 (k-nearest=5, k-farthest=0,
                         small random fraction; NO stranger-joining)
  3. Segment PAF (CIGAR-aware, window 500, gap 1, max span 1000 both axes)
     -> allwave.segmented.paf
  4. impg partition -w 500 -d 0 --min-boundary-distance 0 --min-missing-size 0
     --no-rehome-singletons
       -o bed -> partitions.bed (combined 4-col)
       -o maf --separate-files -> partitions/partition<N>.maf
  5. manifest.json + commands.log (reproducible).

Outputs are written under <outdir>/<clade_id>/.

Usage:
  per_clade_alignment_pipeline.py --build-index --fasta <full_prophages.fa> \
      --index <idx.json>
  per_clade_alignment_pipeline.py \
      --community 10 \
      --outdir research/clades/10 \
      --clades research/clades/10/tight_clades.json \
      --fasta prophage_homology_survey/full_prophages.fa \
      --index research/clades/full_prophages.idx.json \
      --threads 8 --jobs 16 \
      [--clade 10_0000 --clade 10_0134 ... | --all]
"""

import argparse
import json
import os
import statistics
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

ALLWAVE = "/home/erikg/.cargo/bin/allwave"
IMPG = "/home/erikg/.cargo/bin/impg"
SEGMENT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "segment_paf.py")

WINDOW = 500
GAP = 0          # contiguous chunks (matches validated run + README: no inter-chunk gaps)
MAX_SPAN = 2 * WINDOW
ALL_PAIRS_MAX_N = 30
TREE_STRAT = "tree:5:0:0.0"      # k-nearest only, k-farthest=0
BIG_TREE_STRAT = "tree:10:0:0.0"  # k-nearest=10 for large clades (linear 10n pairs,
                                   # better connectivity than k=5; still k-farthest=0)
BIG_CLADE_N = 200
PARTITION_FLAGS = ["-w", str(WINDOW), "-d", "0",
                   "--min-boundary-distance", "0", "--min-missing-size", "0",
                   "-m", "1",   # max transitive depth 1: keeps intervals ~500 bp
                   "--no-rehome-singletons"]


def build_index(fasta_path, index_path):
    """Scan the FASTA once, record byte offset + length per sequence."""
    idx = {}
    with open(fasta_path, "rb") as f:
        name = None
        start = None
        length = 0
        for line in f:
            if line.startswith(b">"):
                if name is not None:
                    idx[name] = [start, length]
                name = line[1:].strip().split()[0].decode()
                start = f.tell()
                length = 0
            else:
                length += len(line.strip())
        if name is not None:
            idx[name] = [start, length]
    tmp = index_path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(idx, f)
    os.replace(tmp, index_path)
    print(f"indexed {len(idx)} sequences -> {index_path}", flush=True)


def run_cmd(cmd, log_path, timeout=None):
    with open(log_path, "a") as f:
        f.write("$ " + " ".join(cmd) + "\n")
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    dt = time.time() - t0
    with open(log_path, "a") as f:
        if r.returncode != 0:
            f.write(f"[exit {r.returncode}] stderr:\n{r.stderr[-4000:]}\n")
        f.write(f"[{dt:.1f}s]\n")
    if r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {' '.join(cmd)}\n"
                           f"{r.stderr[-2000:]}")
    return r, dt


def extract_fasta(fasta_path, idx, members, out_path):
    """Seek-based extraction using the offset index."""
    with open(fasta_path, "rb") as fin, open(out_path, "w") as fout:
        for m in members:
            off, ln = idx[m]
            fin.seek(off)
            data = fin.read(ln)
            fout.write(">" + m + "\n")
            # line-wrap at 80 like the original (data has no newlines stripped)
            for i in range(0, len(data), 80):
                fout.write(data[i:i + 80].decode() + "\n")
    return len(members)


def seq_lengths(fasta_path):
    lens = []
    with open(fasta_path) as f:
        cur = 0
        for line in f:
            if line.startswith(">"):
                if cur:
                    lens.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
        if cur:
            lens.append(cur)
    return lens


def paf_stats(paf_path):
    nlines = 0
    names = set()
    pairs = set()
    with open(paf_path) as f:
        for line in f:
            c = line.split("\t")
            nlines += 1
            names.add(c[0])
            names.add(c[5])
            pairs.add(tuple(sorted((c[0], c[5]))))
    return nlines, names, pairs


def bed_interval_stats(bed_path):
    spans = []
    pids = set()
    with open(bed_path) as f:
        for line in f:
            c = line.split("\t")
            if len(c) < 4:
                continue
            spans.append(int(c[2]) - int(c[1]))
            pids.add(c[3].strip())
    if not spans:
        return {"n_intervals": 0, "n_partitions": 0, "min": None,
                "median": None, "max": None, "n_gt_1000": 0, "n_lt_100": 0}
    ss = sorted(spans)
    return {
        "n_intervals": len(spans),
        "n_partitions": len(pids),
        "min": ss[0],
        "median": statistics.median(spans),
        "max": ss[-1],
        "n_gt_1000": sum(1 for s in spans if s > 1000),
        "n_lt_100": sum(1 for s in spans if s < 100),
    }


def process_clade(args, cid, members):
    outdir = os.path.join(args.outdir, cid)
    os.makedirs(outdir, exist_ok=True)
    log_path = os.path.join(outdir, "commands.log")

    # resume: skip clades whose outputs are already complete
    manifest_path = os.path.join(outdir, "manifest.json")
    if os.path.exists(manifest_path):
        try:
            existing = json.load(open(manifest_path))
        except Exception:
            existing = None
        if existing and existing.get("pipeline", {}).get("total_runtime_s") is not None:
            return cid, existing

    start = time.time()
    manifest = {
        "clade_id": cid,
        "community": args.community,
        "n_members": len(members),
        "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "strategy": None,
        "params": {},
        "pipeline": {},
    }

    # 1. extract (seek-based)
    fasta = os.path.join(outdir, "sequences.fa")
    idx = _IDX_CACHE
    n = extract_fasta(args.fasta, idx, members, fasta)
    lens = seq_lengths(fasta)
    manifest["sequence_lengths"] = {
        "min": min(lens), "median": statistics.median(lens),
        "max": max(lens), "total_bp": sum(lens),
    }

    # 2. align
    paf = os.path.join(outdir, "allwave.paf")
    if n == 1:
        manifest["strategy"] = "none (singleton, no alignment)"
        manifest["pipeline"]["allwave"] = {"pairs_aligned": 0,
                                           "possible_pairs": 0,
                                           "alignment_rate": None,
                                           "sequences_in_paf": 0,
                                           "runtime_s": 0.0}
        with open(paf, "w") as f:
            pass
    else:
        if n <= ALL_PAIRS_MAX_N:
            strat = "none"
        elif n <= BIG_CLADE_N:
            strat = TREE_STRAT
        else:
            strat = BIG_TREE_STRAT
        k_nearest = 5 if strat == TREE_STRAT else (10 if strat == BIG_TREE_STRAT else None)
        manifest["strategy"] = strat
        manifest["params"]["allwave"] = {
            "sparsification": strat, "k_nearest": k_nearest, "k_farthest": 0,
            "random_fraction": None,
            "threads": args.threads, "scores": "0,5,8,2,24,1",
            "note": "k-nearest only, k-farthest=0 (no stranger-joining); "
                    "k=10 for n>200 to keep the k-nearest graph connected at "
                    "linear pair count (10n)",
        }
        cmd = [ALLWAVE, "-i", fasta, "-o", paf, "-p", strat,
               "-t", str(args.threads), "--no-progress"]
        _, dt = run_cmd(cmd, log_path)
        nlines, names, pairs = paf_stats(paf)
        possible_pairs = n * (n - 1) // 2
        manifest["pipeline"]["allwave"] = {
            "pairs_aligned": len(pairs),
            "directed_records": nlines,
            "possible_pairs": possible_pairs,
            "alignment_rate": round(len(pairs) / possible_pairs, 6),
            "sequences_in_paf": len(names),
            "runtime_s": round(dt, 1),
        }

    # 3. segment
    seg_paf = os.path.join(outdir, "allwave.segmented.paf")
    if n > 1:
        cmd = ["python3", SEGMENT, "-i", paf, "-o", seg_paf,
               "-w", str(WINDOW), "--gap", str(GAP), "--max-span", str(MAX_SPAN),
               "--log", log_path]
        _, dt = run_cmd(cmd, log_path)
        nseg = sum(1 for _ in open(seg_paf))
        manifest["params"]["segment"] = {
            "window": WINDOW, "gap": GAP, "max_span": MAX_SPAN,
            "n_chunks": nseg, "runtime_s": round(dt, 1),
        }
    else:
        with open(seg_paf, "w") as f:
            pass
        manifest["params"]["segment"] = {"window": WINDOW, "gap": GAP,
                                         "max_span": MAX_SPAN, "n_chunks": 0}

    # 4. partition
    if n > 1:
        bed = os.path.join(outdir, "partitions.bed")
        mafdir = os.path.join(outdir, "partitions")
        os.makedirs(mafdir, exist_ok=True)
        cmd = [IMPG, "partition", "-a", seg_paf] + PARTITION_FLAGS + \
              ["-o", "bed", "--output-folder", outdir, "--temp-dir", outdir]
        _, dt = run_cmd(cmd, log_path)
        manifest["pipeline"]["partition_bed"] = {
            "runtime_s": round(dt, 1),
            "interval_stats": bed_interval_stats(bed),
        }
        cmd = [IMPG, "partition", "-a", seg_paf] + PARTITION_FLAGS + \
              ["-o", "maf", "--separate-files", "--output-folder", mafdir,
               "--sequence-files", fasta, "--temp-dir", outdir]
        _, dt = run_cmd(cmd, log_path)
        nmaf = sum(1 for x in os.listdir(mafdir) if x.endswith(".maf"))
        manifest["pipeline"]["partition_maf"] = {
            "runtime_s": round(dt, 1), "n_maf_partitions": nmaf,
        }
        # clean up impg index
        for junk in (seg_paf + ".impg",):
            if os.path.exists(junk):
                os.remove(junk)
    else:
        manifest["pipeline"]["partition"] = "skipped (singleton)"

    manifest["pipeline"]["total_runtime_s"] = round(time.time() - start, 1)
    with open(os.path.join(outdir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    return cid, manifest


_IDX_CACHE = None


def _init_worker(args):
    global _IDX_CACHE
    with open(args.index) as f:
        _IDX_CACHE = json.load(f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-index", action="store_true")
    ap.add_argument("--community", default=None)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--clades", default=None, help="tight_clades.json")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--clade", action="append")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()

    if args.build_index:
        build_index(args.fasta, args.index)
        return 0

    if not (args.clades and args.outdir and args.community is not None):
        print("need --community, --outdir, --clades (or --build-index)",
              file=sys.stderr)
        return 1

    tc = json.load(open(args.clades))
    cids = args.clade or (list(tc.keys()) if args.all else [])
    if not cids:
        print("no clades selected (use --all or --clade)", file=sys.stderr)
        return 1

    results = {}
    if args.jobs > 1 and len(cids) > 1:
        with ProcessPoolExecutor(max_workers=args.jobs,
                                 initializer=_init_worker,
                                 initargs=(args,)) as ex:
            futs = {ex.submit(process_clade, args, cid, tc[cid]): cid
                    for cid in cids}
            for fut in futs:
                cid, manifest = fut.result()
                results[cid] = manifest
                print(f"done {cid}: n={manifest['n_members']} "
                      f"strategy={manifest['strategy']}", flush=True)
    else:
        _init_worker(args)
        for cid in cids:
            _, manifest = process_clade(args, cid, tc[cid])
            results[cid] = manifest
            print(f"done {cid}: n={manifest['n_members']} "
                  f"strategy={manifest['strategy']}", flush=True)

    with open(os.path.join(args.outdir, "pipeline_results.json"), "w") as f:
        json.dump(results, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
