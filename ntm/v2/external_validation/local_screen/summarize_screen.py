#!/usr/bin/env python3
"""Summarize tier results + merge with pilot external-evidence categories.

Outputs:
  results/<tier>_hits.tsv          per-bait screen summary
  results/genome_external_evidence_update.tsv  per-bait best evidence across
    tiers, merged with ntm/v2/pilot/validation/bait_external_evidence.tsv
    SRA column (category ladder from the frozen RUN_PLAN §8).
  results/screen_summary.json      run-level receipts
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict

ROLES = {
    "POSCON_D29_AF022214_2_MID1200": "positive_control",
    "POSCON_L5_NC_001335_1_MID1200": "positive_control",
    "NTMBAIT_SHUF_CONTROL_01": "negative_control",
    "NTMBAIT_SHUF_CONTROL_02": "negative_control",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", nargs="+", required=True, help="tier jsonl files")
    ap.add_argument("--tiers", nargs="+", required=True)
    ap.add_argument("--pilot-bait-tsv", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--threshold", type=float, default=0.7)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    rows = []
    for res, tier in zip(args.results, args.tiers):
        for line in open(res):
            r = json.loads(line)
            r["tier"] = r.get("tier") or tier
            rows.append(r)
    print(f"{len(rows)} screened accessions")

    # negative gate over ALL screened accessions
    negs = ["NTMBAIT_SHUF_CONTROL_01", "NTMBAIT_SHUF_CONTROL_02"]
    neg_viol = [(r["accession"], b, r["per_bait_kmer_cov"][b]) for r in rows
                for b in negs if r["per_bait_kmer_cov"].get(b, 0) >= args.threshold]

    # per-bait aggregation
    baits = sorted(rows[0]["per_bait_kmer_cov"])
    agg = {}
    for b in baits:
        covs = sorted(((r["per_bait_kmer_cov"][b], r["accession"]) for r in rows), reverse=True)
        agg[b] = {
            "role": ROLES.get(b, "ntm_bait"),
            "n_accessions": len(rows),
            "n_cov_gt_0": sum(1 for c, _ in covs if c > 0),
            "n_cov_ge_0.5": sum(1 for c, _ in covs if c >= 0.5),
            "n_cov_ge_thr": sum(1 for c, _ in covs if c >= args.threshold),
            "max_cov": covs[0][0],
            "best_accession": covs[0][1],
            "top5": [(c, a) for c, a in covs[:5] if c > 0],
        }

    with open(os.path.join(args.out_dir, "screen_hits_by_bait.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["bait_id", "role", "n_accessions_screened", "n_cov_gt_0",
                    "n_cov_ge_0.5", "n_cov_ge_0.7", "max_kmer_cov",
                    "best_accession", "top5_accessions_cov"])
        for b in sorted(agg, key=lambda x: -agg[x]["n_cov_ge_thr"]):
            a = agg[b]
            w.writerow([b, a["role"], a["n_accessions"], a["n_cov_gt_0"],
                        a["n_cov_ge_0.5"], a["n_cov_ge_thr"], f"{a['max_cov']:.4f}",
                        a["best_accession"],
                        ";".join(f"{x[1]}:{x[0]:.3f}" for x in a["top5"])])

    # merge with pilot bait external evidence (SRA axis)
    pilot = {}
    with open(args.pilot_bait_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pilot[r["bait_id"]] = r
    cols = ["bait_id", "role", "n_accessions_screened_A_B",
            "n_screen_hits_ge_0.7", "best_accession", "best_kmer_cov",
            "sra_evidence_category", "pilot_public_axis_best_ref",
            "pilot_public_axis_cov_identity"]
    with open(os.path.join(args.out_dir, "genome_external_evidence_update.tsv"),
              "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, lineterminator="\n")
        w.writeheader()
        for b in sorted(agg):
            a = agg[b]
            p = pilot.get(b, {})
            if a["role"] == "negative_control":
                cat = "negative_control_clean" if not neg_viol else "NEGATIVE_CONTROL_VIOLATION"
            elif a["role"] == "positive_control":
                cat = "positive_control_validated"
            elif a["n_cov_ge_thr"] == 0:
                cat = "no hit (not detected in screened subset)"
            elif a["n_cov_ge_thr"] > 0:
                cat = ("accession-level screen hit"
                       + (" (multi-accession)" if a["n_cov_ge_thr"] > 1 else ""))
            w.writerow({
                "bait_id": b, "role": a["role"],
                "n_accessions_screened_A_B": a["n_accessions"],
                "n_screen_hits_ge_0.7": a["n_cov_ge_thr"],
                "best_accession": a["best_accession"],
                "best_kmer_cov": f"{a['max_cov']:.4f}",
                "sra_evidence_category": cat,
                "pilot_public_axis_best_ref": f"{p.get('best_public_ref','')}|{p.get('best_public_acc','')}",
                "pilot_public_axis_cov_identity": f"{p.get('cov_frac','')}/{p.get('identity','')}",
            })

    # per-accession hit list
    hits = []
    for r in rows:
        hb = {b: v for b, v in r["per_bait_kmer_cov"].items()
              if v >= args.threshold and b not in ROLES}
        if hb:
            hits.append({"accession": r["accession"], "tier": r["tier"],
                         "scientific_name": r.get("scientific_name", ""),
                         "baits": hb})
    with open(os.path.join(args.out_dir, "hit_accessions.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "tier", "scientific_name", "n_baits_ge_0.7",
                    "baits_ge_0.7"])
        for h in hits:
            w.writerow([h["accession"], h["tier"], h["scientific_name"],
                        len(h["baits"]),
                        ";".join(f"{b}:{v:.3f}" for b, v in sorted(h["baits"].items()))])

    summary = {
        "n_screened": len(rows),
        "n_hit_accessions": len(hits),
        "negative_gate_violations": neg_viol,
        "per_bait": {b: {k: v for k, v in a.items() if k != "top5"}
                     for b, a in agg.items()},
        "generated_utc": __import__("time").strftime("%Y-%m-%dT%H:%M:%SZ",
                                                     __import__("time").gmtime()),
    }
    with open(os.path.join(args.out_dir, "screen_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"hit accessions: {len(hits)}; negative violations: {len(neg_viol)}")
    for h in hits[:20]:
        print(f"  {h['accession']} [{h['tier']}] {h['scientific_name']}: "
              f"{sorted(h['baits'])}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
