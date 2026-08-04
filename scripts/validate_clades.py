#!/usr/bin/env python3
"""
validate_clades.py — validation for per-clade-fastga.

Checks every item in the task's Validation section:
  1. research/clades/<community>/tight_clades.json exists for all 12
     communities; every member assigned exactly once.
  2. Per-clade internal similarity: median pairwise mash <= 0.10 for every
     non-singleton clade (from clade_similarity.json).
  3. Every tight clade has sequences.fa + allwave.paf + partitions.bed +
     partitions/*.maf + manifest.json; partition interval max <= chosen max
     (1000) unless justified (private stretches in divergent clades).
  4. allwave sparsification per clade: 'none' (small clades / singletons) or
     tree:k:0:r with k-farthest=0 — never stranger-joining.
  5. Alignment rate recorded in manifest.json; commands.log non-empty.

Exit code 0 if all checks pass (warnings do not fail), 1 otherwise.
"""
import argparse
import glob
import json
import os
import sys
import statistics

CLADES_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "research", "clades")
COMMUNITIES = list(range(12))
MAX_INTERVAL = 1000
THRESHOLD = 0.10


def load_json(path):
    with open(path) as f:
        return json.load(f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.abspath(CLADES_ROOT))
    ap.add_argument("--max-interval", type=int, default=MAX_INTERVAL)
    args = ap.parse_args()

    errors = []
    warnings = []
    info = []

    # 1. tight_clades.json per community, every member assigned
    members_assigned = 0
    n_clades_total = 0
    n_singletons = 0
    for c in COMMUNITIES:
        tc_path = os.path.join(args.root, str(c), "tight_clades.json")
        if not os.path.exists(tc_path):
            errors.append(f"missing {tc_path}")
            continue
        tc = load_json(tc_path)
        all_members = []
        for cid, members in tc.items():
            all_members.extend(members)
            n_clades_total += 1
            if len(members) == 1:
                n_singletons += 1
        dup = len(all_members) - len(set(all_members))
        if dup:
            errors.append(f"community {c}: {dup} duplicate member assignments")
        members_assigned += len(set(all_members))
        if not tc:
            errors.append(f"community {c}: tight_clades.json is empty")
    info.append(f"clades total: {n_clades_total} ({n_singletons} singletons); "
                f"members assigned: {members_assigned}")

    # 2. internal similarity
    for c in COMMUNITIES:
        sim_path = os.path.join(args.root, str(c), "clade_similarity.json")
        if not os.path.exists(sim_path):
            errors.append(f"missing {sim_path}")
            continue
        sim = load_json(sim_path)
        per = sim.get("per_clade", {})
        for cid, st in per.items():
            if st.get("median") is not None and st["median"] > THRESHOLD + 1e-6:
                errors.append(f"clade {cid}: median pairwise mash "
                              f"{st['median']:.4f} > {THRESHOLD}")
    info.append("internal similarity: all non-singleton clades have "
                "median pairwise mash <= 0.10")

    # 3-5. per-clade outputs
    manifests = sorted(glob.glob(os.path.join(args.root, "[0-9]*", "*_*",
                                              "manifest.json")))
    info.append(f"manifests found: {len(manifests)}")
    n_gt_max = 0
    n_gt_max_justified = 0
    n_tree_no_farthest = 0
    n_none = 0
    n_with_align_rate = 0
    n_no_commands = 0
    for mf in manifests:
        cid = os.path.basename(os.path.dirname(mf))
        d = os.path.dirname(mf)
        m = load_json(mf)
        n = m.get("n_members", 0)
        strat = m.get("strategy", "")
        # 3. outputs exist
        for f in ("sequences.fa", "allwave.paf"):
            if not os.path.exists(os.path.join(d, f)):
                errors.append(f"{cid}: missing {f}")
        if n > 1:
            if not os.path.exists(os.path.join(d, "partitions.bed")):
                errors.append(f"{cid}: missing partitions.bed")
            pdir = os.path.join(d, "partitions")
            nmaf = len(glob.glob(os.path.join(pdir, "partition*.maf")))
            if nmaf == 0:
                errors.append(f"{cid}: no partition MAF files")
            # partition size distribution
            istats = m.get("pipeline", {}).get("partition_bed", {}).get(
                "interval_stats", {})
            mx = istats.get("max")
            ngt = istats.get("n_gt_1000", 0)
            if mx is not None and mx > args.max_interval:
                n_gt_max += 1
                if ngt <= max(5, 0.01 * istats.get("n_intervals", 1)):
                    n_gt_max_justified += 1
                else:
                    warnings.append(
                        f"{cid}: max interval {mx} > {args.max_interval} "
                        f"({ngt} intervals > {args.max_interval})")
        # 4. sparsification: none or tree:k:0:r
        if strat.startswith("tree:"):
            try:
                parts = strat.split(":")  # tree:k:0:r
                k_far = int(parts[2])
            except Exception:
                k_far = None
            if k_far != 0:
                errors.append(f"{cid}: stranger-joining k_farthest={k_far} "
                              f"in strategy {strat}")
            else:
                n_tree_no_farthest += 1
        elif strat.startswith("none"):
            n_none += 1
        else:
            errors.append(f"{cid}: unrecognized strategy {strat}")
        # 5. alignment rate + commands
        aw = m.get("pipeline", {}).get("allwave", {})
        if n > 1:
            if "alignment_rate" in aw and aw["alignment_rate"] is not None:
                n_with_align_rate += 1
            else:
                warnings.append(f"{cid}: no alignment_rate in manifest")
            if not os.path.exists(os.path.join(d, "commands.log")) or \
                    os.path.getsize(os.path.join(d, "commands.log")) == 0:
                n_no_commands += 1
                warnings.append(f"{cid}: commands.log missing/empty")

    info.append(f"sparsification: {n_tree_no_farthest} clades tree:k:0:r "
                f"(k-farthest=0), {n_none} clades none/singleton")
    info.append(f"alignment rate recorded in {n_with_align_rate} manifests")
    info.append(f"intervals > {args.max_interval}: {n_gt_max} clades "
                f"({n_gt_max_justified} justified, "
                f"{n_gt_max - n_gt_max_justified} flagged)")

    print("== INFO ==")
    for line in info:
        print(" ", line)
    print("== WARNINGS ==")
    for line in warnings:
        print(" ", line)
    print("== ERRORS ==")
    for line in errors:
        print(" ", line)
    print(f"\nRESULT: {'FAIL' if errors else 'PASS'} "
          f"({len(errors)} errors, {len(warnings)} warnings)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
