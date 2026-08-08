#!/usr/bin/env python3
"""
resume_per_clade_partition.py — resume the per-clade NTM pipeline from the
allwave.paf stage.

The full per-clade run (scripts/per_clade_alignment_pipeline.py over all 913
tight prophage clades) completed for 898 clades but was interrupted during the
remaining 15 after producing sequences.fa + allwave.paf (before segment /
partition / manifest).

This resumes those clades from the existing allwave.paf (no re-alignment),
running exactly the same segment -> impg partition (bed + maf) -> manifest
steps as process_clade(), then writes the aggregate pipeline_results.json over
ALL clades (reading each clade's manifest.json).

Output: ntm/v1/clades/<clade>/ (symlinked to mash_clades/clades/), plus
        ntm/v1/clades/pipeline_results.json.

Usage:
  resume_per_clade_partition.py \
      --clades <tight_clades.json> --outdir <clades_root> --community 0 \
      --fasta <full_prophages.fa> --threads 8 \
      [--clade 0_0097 --clade 0_0101 ... | --all-resume]
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "scripts"))
from per_clade_alignment_pipeline import (  # noqa: E402
    IMPG, SEGMENT, WINDOW, GAP, MAX_SPAN, PARTITION_FLAGS,
    run_cmd, seq_lengths, bed_interval_stats, paf_stats,
)


def resume_clade(args, cid, members):
    outdir = os.path.join(args.outdir, cid)
    os.makedirs(outdir, exist_ok=True)
    log_path = os.path.join(outdir, "commands.log")

    manifest_path = os.path.join(outdir, "manifest.json")
    if os.path.exists(manifest_path):
        try:
            existing = json.load(open(manifest_path))
        except Exception:
            existing = None
        if existing and existing.get("pipeline", {}).get("total_runtime_s") is not None:
            return cid, existing

    outdir_abs = os.path.abspath(outdir)
    paf = os.path.join(outdir_abs, "allwave.paf")
    if not os.path.exists(paf) or os.path.getsize(paf) == 0:
        raise RuntimeError(
            f"clade {cid}: allwave.paf missing/empty — cannot resume from "
            f"allwave stage; rerun per_clade_alignment_pipeline.py instead")

    start = time.time()
    n = len(members)
    fasta = os.path.join(outdir_abs, "sequences.fa")
    lens = seq_lengths(fasta)
    ALL_PAIRS_MAX_N = 30
    BIG_CLADE_N = 200
    if n <= ALL_PAIRS_MAX_N:
        strat = "none"
    elif n <= BIG_CLADE_N:
        strat = "tree:5:0:0.0"
    else:
        strat = "tree:10:0:0.0"
    k_nearest = 5 if strat == "tree:5:0:0.0" else (10 if strat == "tree:10:0:0.0" else None)
    manifest = {
        "clade_id": cid,
        "community": args.community,
        "n_members": n,
        "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "strategy": strat,
        "params": {},
        "pipeline": {},
    }
    manifest["sequence_lengths"] = {
        "min": min(lens), "median": float(statistics_median(lens)),
        "max": max(lens), "total_bp": sum(lens),
    }

    # allwave stage already done -> fill from paf
    nlines, names, pairs = paf_stats(paf)
    possible_pairs = n * (n - 1) // 2
    manifest["params"]["allwave"] = {
        "sparsification": manifest["strategy"],
        "k_nearest": k_nearest,
        "k_farthest": 0,
        "random_fraction": None,
        "threads": args.threads,
        "scores": "0,5,8,2,24,1",
        "note": "k-nearest only, k-farthest=0 (no stranger-joining); "
                "k=10 for n>200 to keep the k-nearest graph connected at "
                "linear pair count (10n)",
    }
    manifest["pipeline"]["allwave"] = {
        "pairs_aligned": len(pairs),
        "directed_records": nlines,
        "possible_pairs": possible_pairs,
        "alignment_rate": round(len(pairs) / possible_pairs, 6),
        "sequences_in_paf": len(names),
        "runtime_s": None,  # alignment ran in the prior session
    }

    # segment
    seg_paf = os.path.join(outdir_abs, "allwave.segmented.paf")
    cmd = ["python3", SEGMENT, "-i", paf, "-o", seg_paf,
           "-w", str(WINDOW), "--gap", str(GAP), "--max-span", str(MAX_SPAN),
           "--log", log_path]
    _, dt = run_cmd(cmd, log_path)
    nseg = sum(1 for _ in open(seg_paf))
    manifest["params"]["segment"] = {
        "window": WINDOW, "gap": GAP, "max_span": MAX_SPAN,
        "n_chunks": nseg, "runtime_s": round(dt, 1),
    }

    # partition -> bed
    bed = os.path.join(outdir_abs, "partitions.bed")
    cmd = [IMPG, "partition", "-a", seg_paf] + PARTITION_FLAGS + \
          ["-o", "bed", "--output-folder", outdir_abs, "--temp-dir", outdir_abs]
    _, dt = run_cmd(cmd, log_path)
    manifest["pipeline"]["partition_bed"] = {
        "runtime_s": round(dt, 1),
        "interval_stats": bed_interval_stats(bed),
    }

    # partition -> maf (separate files)
    mafdir = os.path.join(outdir_abs, "partitions")
    os.makedirs(mafdir, exist_ok=True)
    cmd = [IMPG, "partition", "-a", seg_paf] + PARTITION_FLAGS + \
          ["-o", "maf", "--separate-files", "--output-folder", mafdir,
           "--sequence-files", fasta, "--temp-dir", outdir_abs]
    _, dt = run_cmd(cmd, log_path)
    nmaf = sum(1 for x in os.listdir(mafdir) if x.endswith(".maf"))
    manifest["pipeline"]["partition_maf"] = {
        "runtime_s": round(dt, 1), "n_maf_partitions": nmaf,
    }
    for junk in (seg_paf + ".impg",):
        if os.path.exists(junk):
            os.remove(junk)

    manifest["pipeline"]["total_runtime_s"] = round(time.time() - start, 1)
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=1)
    return cid, manifest


def statistics_median(xs):
    import statistics
    return statistics.median(xs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--community", default=None)
    ap.add_argument("--outdir", required=True, help="clades root dir")
    ap.add_argument("--clades", required=True, help="tight_clades.json")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--clade", action="append", dest="clade_ids")
    ap.add_argument("--all-resume", action="store_true",
                    help="resume every clade lacking a complete manifest")
    args = ap.parse_args()

    tc = json.load(open(args.clades))
    cids = args.clade_ids
    if args.all_resume:
        cids = [c for c in tc if not _complete(args.outdir, c)]
    if not cids:
        print("nothing to resume", file=sys.stderr)
        return 1

    results = {}
    for cid in cids:
        _, manifest = resume_clade(args, cid, tc[cid])
        results[cid] = manifest
        print(f"done {cid}: n={manifest['n_members']} "
              f"strategy={manifest['strategy']} "
              f"median_span={manifest['pipeline']['partition_bed']['interval_stats']['median']}")
        sys.stdout.flush()

    # aggregate over ALL clades from manifests
    agg = {}
    for cid in tc:
        mp = os.path.join(args.outdir, cid, "manifest.json")
        if os.path.exists(mp):
            try:
                agg[cid] = json.load(open(mp))
            except Exception:
                agg[cid] = {"clade_id": cid, "error": "unreadable manifest"}
        else:
            agg[cid] = {"clade_id": cid, "error": "missing manifest"}
    with open(os.path.join(args.outdir, "pipeline_results.json"), "w") as f:
        json.dump(agg, f, indent=1)
    print(f"wrote {os.path.join(args.outdir, 'pipeline_results.json')} "
          f"({len(agg)} clades)")
    return 0


def _complete(outdir, cid):
    mp = os.path.join(outdir, cid, "manifest.json")
    if not os.path.exists(mp):
        return False
    try:
        return json.load(open(mp)).get("pipeline", {}).get("total_runtime_s") is not None
    except Exception:
        return False


if __name__ == "__main__":
    sys.exit(main())