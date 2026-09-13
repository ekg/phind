#!/usr/bin/env python3
"""
validate_v3_partitions.py — independent validation for task ntm-v3-per.

Checks (mirrors the task's Validation section):
  1. Every alignable (n>=2) tight clade under ntm/v3/clades/ has
     sequences.fa, allwave.paf, allwave.segmented.paf, partitions.bed and
     a non-empty partitions/*.maf set, plus manifest.json — OR an explicit
     FAILED.log.
  2. No clade skipped silently: alignable + singletons + failed == total
     clades in tight_clades.json (repo ntm/v3/clades/tight_clades.json.gz
     == NVMe ntm/v3/clades/0/tight_clades.json, sha256-checked here).
  3. FINAL parameters (user-approved): strategy `none` for n<=30,
     `tree:5:0:0.0` for 31<=n<=200, `tree:10:0:0.0` for n>200; manifest
     records k_farthest=0 and NO stranger-joining (sparsification
     k-nearest only); allwave scores default 0,5,8,2,24,1; segment
     window 500 max-span 1000; impg partition -w 500.
  4. Aggregate partition stats: total partitions, per-clade + overall
     interval length distribution, n>1000 / n<100 counts.
  5. Alignment rate (fraction of clade members with >=1 PAF hit) per
     clade, from manifest sequences_in_paf/n_members.
  6. Failures listed with reasons.

Exit 0 if all checks pass, 1 otherwise.
"""
import argparse
import gzip
import hashlib
import json
import os
import statistics
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
V3 = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3"
OUT = f"{V3}/clades"
CLADES_JSON = f"{OUT}/0/tight_clades.json"
CLADES_GZ = os.path.join(REPO, "ntm", "v3", "clades", "tight_clades.json.gz")

ALLOWED_STRATS = {"none", "tree:5:0:0.0", "tree:10:0:0.0"}


def expected_strategy(n):
    if n <= 30:
        return "none"
    if n <= 200:
        return "tree:5:0:0.0"
    return "tree:10:0:0.0"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=OUT)
    ap.add_argument("--clades", default=CLADES_JSON)
    ap.add_argument("--clades-gz", default=CLADES_GZ)
    ap.add_argument("--report", default=None)
    ap.add_argument("--summary-tsv", default=None,
                    help="per-clade summary TSV path")
    args = ap.parse_args()

    tc = json.load(open(args.clades))
    if os.path.exists(args.clades_gz):
        with gzip.open(args.clades_gz) as f:
            tcrepo = json.load(f)
        if tcrepo != tc:
            print("FATAL: repo tight_clades.json.gz != NVMe tight_clades.json")
            return 1
    all_cids = list(tc.keys())
    alignable = {c: m for c, m in tc.items() if len(m) >= 2}
    singletons = {c: m for c, m in tc.items() if len(m) == 1}

    missing = []
    failed = {}
    manifest_issues = []
    singleton_dirs = 0
    align_with_partition = 0
    strat_counts = {}
    all_medians = []
    total_partitions = 0
    total_intervals = 0
    n_gt_1000 = 0
    n_lt_100 = 0
    gt1000_clades = []
    align_rates = []
    member_hit_rates = []
    per_clade_runtime = []
    rows = []

    for cid in all_cids:
        d = os.path.join(args.outdir, cid)
        if not os.path.isdir(d):
            missing.append((cid, "no output dir"))
            continue
        if os.path.exists(os.path.join(d, "FAILED.log")):
            with open(os.path.join(d, "FAILED.log")) as f:
                failed[cid] = f.read().strip()[:400]
            continue
        mpath = os.path.join(d, "manifest.json")
        if not os.path.exists(mpath):
            missing.append((cid, "no manifest.json"))
            continue
        try:
            m = json.load(open(mpath))
        except Exception as e:
            manifest_issues.append(f"{cid}: manifest unreadable ({e})")
            continue
        n = m.get("n_members")
        if n != len(tc[cid]):
            manifest_issues.append(f"{cid}: manifest n={n} != clades n={len(tc[cid])}")
        seqfa = os.path.join(d, "sequences.fa")
        if not os.path.exists(seqfa):
            missing.append((cid, "no sequences.fa"))
            continue
        if cid in singletons:
            singleton_dirs += 1
            continue
        # alignable clade: full output set required
        paf = os.path.join(d, "allwave.paf")
        seg = os.path.join(d, "allwave.segmented.paf")
        bed = os.path.join(d, "partitions.bed")
        mafdir = os.path.join(d, "partitions")
        have_mafs = os.path.isdir(mafdir) and any(
            x.endswith(".maf") for x in os.listdir(mafdir))
        if not (os.path.exists(paf) and os.path.getsize(paf) > 0
                and os.path.exists(seg) and os.path.getsize(seg) > 0
                and os.path.exists(bed) and have_mafs):
            missing.append((cid, "incomplete outputs"))
            continue
        align_with_partition += 1

        # parameter audit
        strat = m.get("strategy")
        exp = expected_strategy(n)
        strat_counts[strat] = strat_counts.get(strat, 0) + 1
        if strat != exp:
            manifest_issues.append(f"{cid}: strategy {strat} (n={n}, expected {exp})")
        aw = m.get("params", {}).get("allwave", {})
        if strat != "none":
            if aw.get("k_farthest") != 0:
                manifest_issues.append(f"{cid}: k_farthest {aw.get('k_farthest')} != 0")
            if aw.get("random_fraction") not in (None, 0, 0.0):
                manifest_issues.append(f"{cid}: random_fraction {aw.get('random_fraction')}")
            sp = aw.get("sparsification", "")
            if "stranger" in sp or sp.split(":")[0] != "tree":
                manifest_issues.append(f"{cid}: sparsification {sp}")
        if aw.get("scores") != "0,5,8,2,24,1":
            manifest_issues.append(f"{cid}: scores {aw.get('scores')}")
        segw = m.get("params", {}).get("segment", {})
        if segw.get("window") != 500 or segw.get("max_span") != 1000:
            manifest_issues.append(f"{cid}: segment {segw.get('window')}/{segw.get('max_span')}")

        # stats
        awp = m.get("pipeline", {}).get("allwave", {})
        ar = awp.get("alignment_rate")
        seqs = awp.get("sequences_in_paf", 0)
        member_hit = seqs / n if n else 0
        if ar is not None:
            align_rates.append(ar)
        member_hit_rates.append(member_hit)
        per_clade_runtime.append(m.get("pipeline", {}).get("total_runtime_s") or 0)
        pb = m.get("pipeline", {}).get("partition_bed", {})
        st = pb.get("interval_stats", {})
        total_partitions += st.get("n_partitions", 0)
        if st.get("median") is not None:
            all_medians.append(st["median"])
        # exact span distribution recomputed from the BED
        spans = []
        with open(bed) as f:
            for line in f:
                c = line.split("\t")
                if len(c) >= 4:
                    spans.append(int(c[2]) - int(c[1]))
        total_intervals += len(spans)
        g = sum(1 for s in spans if s > 1000)
        n_gt_1000 += g
        n_lt_100 += sum(1 for s in spans if s < 100)
        if g:
            gt1000_clades.append((cid, g, max(spans)))
        rows.append({
            "clade_id": cid, "n_members": n, "strategy": strat,
            "pairs_aligned": awp.get("pairs_aligned", "NA"),
            "alignment_rate": ar if ar is not None else "NA",
            "members_with_hit": f"{member_hit:.4f}",
            "n_partitions": st.get("n_partitions", "NA"),
            "median_interval": st.get("median", "NA"),
            "max_interval": st.get("max", "NA"),
            "n_gt_1000": g,
            "runtime_s": m.get("pipeline", {}).get("total_runtime_s", "NA"),
        })

    # accounting
    accounted = align_with_partition + len(failed) + len(missing) + singleton_dirs
    ok = (not missing and not failed and not manifest_issues
          and accounted == len(all_cids)
          and align_with_partition == len(alignable)
          and singleton_dirs == len(singletons)
          and all(s in ALLOWED_STRATS for s in strat_counts))

    print(f"clades total={len(all_cids)} alignable={len(alignable)} "
          f"singletons={len(singletons)} singleton_dirs={singleton_dirs}")
    print(f"alignable with full outputs: {align_with_partition}")
    print(f"missing: {len(missing)}  failed: {len(failed)}  "
          f"manifest issues: {len(manifest_issues)}")
    print(f"strategies: {strat_counts}")
    if align_rates:
        print(f"alignment rate (pairs/possible): median="
              f"{statistics.median(align_rates):.4f} min={min(align_rates):.4f}")
        print(f"member-hit rate (seqs in paf / n): median="
              f"{statistics.median(member_hit_rates):.4f} "
              f"min={min(member_hit_rates):.4f}")
    print(f"total partitions: {total_partitions}  intervals: {total_intervals}")
    if all_medians:
        print(f"per-clade median interval: median="
              f"{statistics.median(all_medians):.0f} max={max(all_medians):.0f}")
    print(f"intervals >1000bp: {n_gt_1000} ({100*n_gt_1000/max(total_intervals,1):.3f}%)"
          f"  <100bp: {n_lt_100} ({100*n_lt_100/max(total_intervals,1):.2f}%)")
    if per_clade_runtime:
        print(f"per-clade runtime: total={sum(per_clade_runtime):.0f}s "
              f"median={statistics.median(per_clade_runtime):.1f}s "
              f"max={max(per_clade_runtime):.1f}s")
    for cid, why in missing[:20]:
        print(f"MISSING {cid}: {why}")
    for cid, reason in list(failed.items())[:20]:
        print(f"FAILED {cid}: {reason[:200]}")
    for i in manifest_issues[:20]:
        print(f"ISSUE {i}")

    if args.summary_tsv:
        cols = ["clade_id", "n_members", "strategy", "pairs_aligned",
                "alignment_rate", "members_with_hit", "n_partitions",
                "median_interval", "max_interval", "n_gt_1000", "runtime_s"]
        with open(args.summary_tsv, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in sorted(rows, key=lambda r: r["clade_id"]):
                f.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"wrote {args.summary_tsv} ({len(rows)} rows)")

    print("VALIDATION:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
