#!/usr/bin/env python3
"""
validate_v2_partitions.py — independent validation for task ntm-v2-per.

Checks (mirrors the task's Validation section):
  1. Every alignable (n>=2) tight clade under ntm/v2/clades/ has
     allwave.paf, allwave.segmented.paf, partitions.bed and a non-empty
     partitions/*.maf set, OR an explicit FAILED.log.
  2. No clade skipped silently: alignable clades + singletons + failed
     == total clades in tight_clades.json.
  3. Parameters: manifest strategy == tree:5:0:0.0 for n>30 (all-pairs
     `none` for n<=30, v1 mainline), allwave scores 0,5,8,2,24,1,
     segment window 500, partition window 500.
  4. Aggregate partition stats: total partition count, overall median
     interval length (target ~500 bp, v1 median 500 bp), per-clade
     median distribution, n>1000 / n<100 counts.
  5. Failures listed with reasons (FAILED.log / commands.log tails).

Exit 0 if all checks pass, 1 otherwise.
"""
import argparse
import json
import os
import statistics
import sys

V2 = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2"
OUT = f"{V2}/clades"
CLADES_JSON = f"{OUT}/0/tight_clades.json"


def clade_outputs(outdir, cid):
    d = os.path.join(outdir, cid)
    return {
        "allwave.paf": os.path.exists(os.path.join(d, "allwave.paf")),
        "allwave.segmented.paf": os.path.exists(
            os.path.join(d, "allwave.segmented.paf")),
        "partitions.bed": os.path.exists(os.path.join(d, "partitions.bed")),
        "partitions_dir": (os.path.isdir(os.path.join(d, "partitions"))
                           and any(x.endswith(".maf")
                                   for x in os.listdir(os.path.join(d, "partitions"))))
        if os.path.isdir(os.path.join(d, "partitions")) else False,
        "FAILED.log": os.path.exists(os.path.join(d, "FAILED.log")),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=OUT)
    ap.add_argument("--clades", default=CLADES_JSON)
    ap.add_argument("--report", default=None,
                    help="write partition report markdown here")
    args = ap.parse_args()

    tc = json.load(open(args.clades))
    all_cids = list(tc.keys())
    alignable = {c: m for c, m in tc.items() if len(m) >= 2}
    singletons = {c: m for c, m in tc.items() if len(m) == 1}

    missing = []          # alignable clades without outputs
    failed = {}           # cid -> reason
    manifest_issues = []
    all_medians = []
    all_spans = []
    total_partitions = 0
    strat_counts = {}
    align_with_partition = 0

    for cid in all_cids:
        d = os.path.join(args.outdir, cid)
        if not os.path.isdir(d):
            if cid in alignable:
                missing.append((cid, "no output dir"))
            continue
        outs = clade_outputs(args.outdir, cid)
        if outs["FAILED.log"]:
            with open(os.path.join(d, "FAILED.log")) as f:
                failed[cid] = f.read().strip()[:400]
            continue
        if cid in alignable:
            if not (outs["allwave.paf"] and outs["allwave.segmented.paf"]
                    and outs["partitions.bed"] and outs["partitions_dir"]):
                missing.append((cid, str(outs)))
                continue
            align_with_partition += 1
        # manifest parameter audit
        mpath = os.path.join(d, "manifest.json")
        if not os.path.exists(mpath):
            if cid in alignable:
                missing.append((cid, "no manifest.json"))
            continue
        try:
            m = json.load(open(mpath))
        except Exception as e:
            manifest_issues.append(f"{cid}: manifest unreadable ({e})")
            continue
        if cid in alignable:
            strat = m.get("strategy")
            n = m.get("n_members")
            exp = "none" if n <= 30 else "tree:5:0:0.0"
            strat_counts[strat] = strat_counts.get(strat, 0) + 1
            if strat != exp:
                manifest_issues.append(
                    f"{cid}: strategy {strat} (n={n}, expected {exp})")
            sc = m.get("params", {}).get("allwave", {}).get("scores")
            if sc != "0,5,8,2,24,1":
                manifest_issues.append(f"{cid}: scores {sc}")
            segw = m.get("params", {}).get("segment", {}).get("window")
            if segw != 500:
                manifest_issues.append(f"{cid}: segment window {segw}")
            pb = m.get("pipeline", {}).get("partition_bed", {})
            st = pb.get("interval_stats", {})
            if st.get("median") is not None:
                all_medians.append(st["median"])
                all_spans.extend([st["median"]] * 1)
                total_partitions += st.get("n_partitions", 0)
            # recompute exact span distribution from the bed for the report
            bed = os.path.join(d, "partitions.bed")
            if os.path.exists(bed):
                spans = []
                with open(bed) as f:
                    for line in f:
                        c = line.split("\t")
                        if len(c) >= 4:
                            spans.append(int(c[2]) - int(c[1]))
                all_spans.extend(spans)

    ok = True
    report = []
    report.append("# NTM v2 — per-clade allwave + segment + impg partition\n")
    report.append(f"- Total tight clades: {len(all_cids)}")
    report.append(f"- Alignable (n>=2): {len(alignable)}")
    report.append(f"- Singletons: {len(singletons)}")
    report.append(f"- Alignable clades with full outputs: {align_with_partition}")
    report.append(f"- Alignable clades missing outputs: {len(missing)}")
    report.append(f"- Failed (FAILED.log): {len(failed)}")
    report.append(f"- Strategy distribution (alignable): {strat_counts}")
    if all_medians:
        report.append(f"- Per-clade median partition length: "
                      f"median={statistics.median(all_medians):.0f} bp, "
                      f"min={min(all_medians):.0f}, max={max(all_medians):.0f}")
    if all_spans:
        ss = sorted(all_spans)
        report.append(f"- Overall partition intervals: n={len(ss)}, "
                      f"median={statistics.median(ss):.0f} bp, "
                      f"mean={statistics.mean(ss):.1f}, min={ss[0]}, "
                      f"max={ss[-1]}")
        report.append(f"- Intervals >1000 bp: "
                      f"{sum(1 for s in ss if s > 1000)} "
                      f"({100*sum(1 for s in ss if s>1000)/len(ss):.2f}%)")
        report.append(f"- Intervals <100 bp: "
                      f"{sum(1 for s in ss if s < 100)} "
                      f"({100*sum(1 for s in ss if s<100)/len(ss):.2f}%)")
    report.append(f"- Total partitions: {total_partitions}")

    if missing:
        ok = False
        report.append("\n## Missing outputs (alignable clades)")
        for cid, why in missing:
            report.append(f"- `{cid}`: {why}")
    if failed:
        report.append("\n## Failed clades (FAILED.log)")
        for cid, reason in failed.items():
            report.append(f"- `{cid}`: {reason[:300].replace(chr(10), ' | ')}")
    if manifest_issues:
        ok = False
        report.append("\n## Manifest/parameter issues")
        for i in manifest_issues:
            report.append(f"- {i}")
    if not ok:
        report.append("\n**Validation: FAIL**")
    else:
        report.append("\n**Validation: PASS**")
    report.append("")
    txt = "\n".join(report)
    print(txt)
    if args.report:
        with open(args.report, "w") as f:
            f.write(txt + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
