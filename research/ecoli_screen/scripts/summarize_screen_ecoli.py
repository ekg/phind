#!/usr/bin/env python3
"""Summarize tier results for the E. coli local screen + per-clade merge.

Generalized from ntm local_screen/summarize_screen.py; diffs:
  - ROLES -> ECBAIT controls (poscon / shuffled-neg / host-neg diagnostic)
  - pilot merge axis -> per-genome functional report
    (research/phage_annotation/per_genome_annotation_qc.tsv, source=ecoli_ml)
  - effective-independence audit: hits collapsed by BioProject and by
    identical download sha256 (clone/identical-assembly detection)
Outputs:
  results/screen_hits_by_bait.tsv
  results/genome_external_evidence_update.tsv   (per-clade merge)
  results/hit_accessions.tsv
  results/screen_summary.json
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
from collections import defaultdict

ROLES = {
    "ECBAIT_POSCON_LAMBDA_01": "positive_control",
    "ECBAIT_POSCON_HK97_02": "positive_control",
    "ECBAIT_SHUF_CONTROL_01": "negative_control",
    "ECBAIT_SHUF_CONTROL_02": "negative_control",
    "ECBAIT_HOST_CONTROL_01": "host_control_diagnostic",
    "ECBAIT_HOST_CONTROL_02": "host_control_diagnostic",
}
NEGS = ["ECBAIT_SHUF_CONTROL_01", "ECBAIT_SHUF_CONTROL_02"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", nargs="+", required=True, help="tier jsonl files")
    ap.add_argument("--tiers", nargs="+", required=True)
    ap.add_argument("--functional-qc-tsv", required=True,
                    help="research/phage_annotation/per_genome_annotation_qc.tsv")
    ap.add_argument("--bait-manifest", required=True,
                    help="research/ecoli_bait/panel/bait_manifest.tsv")
    ap.add_argument("--availability", help="downloads/availability.tsv (sha256)")
    ap.add_argument("--runs-tsv", nargs="*",
                    help="frozen tier TSVs to join bioproject per run")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--threshold", type=float, default=0.7)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    proj_by_run = {}
    for p in (args.runs_tsv or []):
        for r in csv.DictReader(open(p), delimiter="\t"):
            proj_by_run.setdefault(r["run"], r.get("bioproject", ""))

    rows = []
    for res, tier in zip(args.results, args.tiers):
        for line in open(res):
            r = json.loads(line)
            r["tier"] = r.get("tier") or tier
            rows.append(r)
    print(f"{len(rows)} screened accessions")

    neg_viol = [(r["accession"], b, r["per_bait_kmer_cov"][b]) for r in rows
                for b in NEGS if r["per_bait_kmer_cov"].get(b, 0) >= args.threshold]

    baits = sorted(rows[0]["per_bait_kmer_cov"])
    agg = {}
    for b in baits:
        covs = sorted(((r["per_bait_kmer_cov"][b], r["accession"]) for r in rows), reverse=True)
        agg[b] = {
            "role": ROLES.get(b, "ecoli_bait"),
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
            w.writerow([x for x in [b, a["role"], a["n_accessions"], a["n_cov_gt_0"],
                        a["n_cov_ge_0.5"], a["n_cov_ge_thr"], f"{a['max_cov']:.4f}",
                        a["best_accession"],
                        ";".join(f"{x[1]}:{x[0]:.3f}" for x in a["top5"]) if a["top5"] else "-"]])

    # per-accession hit list (non-control baits only)
    hits = []
    for r in rows:
        hb = {b: v for b, v in r["per_bait_kmer_cov"].items()
              if v >= args.threshold and b not in ROLES}
        if hb:
            hits.append({"accession": r["accession"], "tier": r["tier"],
                         "scientific_name": r.get("scientific_name", ""),
                         "bioproject": r.get("bioproject") or proj_by_run.get(r["accession"], ""),
                         "baits": hb})
    with open(os.path.join(args.out_dir, "hit_accessions.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "tier", "scientific_name", "bioproject",
                    "n_baits_ge_0.7", "baits_ge_0.7"])
        for h in hits:
            w.writerow([h["accession"], h["tier"], h["scientific_name"],
                        h["bioproject"], len(h["baits"]),
                        ";".join(f"{b}:{v:.3f}" for b, v in sorted(h["baits"].items()))])

    # ---- effective-independence audit (BioProject + identical sha256 collapse)
    sha_by_acc = {}
    if args.availability and os.path.exists(args.availability):
        for r in csv.DictReader(open(args.availability), delimiter="\t"):
            if r.get("sha256") and r.get("status") == "200":
                sha_by_acc.setdefault(r["accession"], r["sha256"])
    proj_of = {h["accession"]: (h["bioproject"] or "(none)") for h in hits}
    proj_coll = defaultdict(set)
    sha_coll = defaultdict(set)
    for h in hits:
        proj_coll[proj_of[h["accession"]]].add(h["accession"])
        if h["accession"] in sha_by_acc:
            sha_coll[sha_by_acc[h["accession"]]].add(h["accession"])
    dup_sha = {s: sorted(a) for s, a in sha_coll.items() if len(a) > 1}
    independence = {
        "n_hit_accessions": len(hits),
        "n_distinct_bioprojects": len(proj_coll),
        "top_bioprojects": sorted(((p, len(a)) for p, a in proj_coll.items()),
                                  key=lambda kv: -kv[1])[:15],
        "n_distinct_download_sha256": len(sha_coll),
        "n_sha256_groups_with_multiple_runs": len(dup_sha),
        "identical_assembly_groups": {s[:16] + "…": a[:10]
                                      for s, a in list(dup_sha.items())[:20]},
        "note": "per-clade/per-bait 'n accessions' is an upper bound on "
                "independent observations (clone complexes + re-sequenced "
                "isolates collapse; identical zst sha256 = identical file)",
    }

    # ---- per-clade external-evidence merge with the functional report
    bait_meta = {}
    for r in csv.DictReader(open(args.bait_manifest), delimiter="\t"):
        bait_meta[r["bait_id"]] = r
    qc = {}
    for r in csv.DictReader(open(args.functional_qc_tsv), delimiter="\t"):
        if r.get("source") == "ecoli_ml":
            qc[r["genome_id"]] = r
    cols = ["genome_id", "clade_id", "length", "flag", "completeness_pct",
            "host_genes", "functional_tier_note",
            "n_accessions_screened", "n_screen_hits_ge_0.7_any_bait",
            "n_baits_with_hits", "best_kmer_cov", "best_accession",
            "sra_evidence_category"]
    with open(os.path.join(args.out_dir, "genome_external_evidence_update.tsv"),
              "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, lineterminator="\n")
        w.writeheader()
        # clade-level aggregation from bait-level agg
        clade_agg = defaultdict(lambda: {"n_hits": 0, "baits": set(),
                                         "best": (0.0, "")})
        for b, a in agg.items():
            if a["role"] != "ecoli_bait":
                continue
            m = bait_meta.get(b, {})
            cid = m.get("clade_id", "")
            if not cid:
                continue
            qc_id = f"clade_{cid}_ML"  # functional-QC genome_id / clade_id form
            ca = clade_agg[qc_id]
            ca["n_hits"] += a["n_cov_ge_thr"]
            if a["n_cov_ge_thr"]:
                ca["baits"].add(b)
            if a["max_cov"] > ca["best"][0]:
                ca["best"] = (a["max_cov"], a["best_accession"])
        for gid, r in sorted(qc.items()):
            cid = r["clade_id"]
            ca = clade_agg.get(cid, {"n_hits": 0, "baits": set(), "best": (0.0, "")})
            if ca["n_hits"] == 0:
                cat = "no hit (not detected in screened subset)"
            else:
                # refine with confirm_tier_hits_ecoli classifications if run
                cat = ("SRA_screen_hit (" + str(len(ca["baits"])) +
                       " bait(s); see confirm classifications)")
            w.writerow({
                "genome_id": gid, "clade_id": cid, "length": r["length"],
                "flag": r["flag"], "completeness_pct": r["completeness_pct"],
                "host_genes": r["host_genes"], "functional_tier_note": "",
                "n_accessions_screened": len(rows),
                "n_screen_hits_ge_0.7_any_bait": ca["n_hits"],
                "n_baits_with_hits": len(ca["baits"]),
                "best_kmer_cov": f"{ca['best'][0]:.4f}",
                "best_accession": ca["best"][1],
                "sra_evidence_category": cat,
            })

    summary = {
        "n_screened": len(rows),
        "n_hit_accessions": len(hits),
        "negative_gate_violations": neg_viol,
        "effective_independence": independence,
        "per_bait": {b: {k: v for k, v in a.items() if k != "top5"}
                     for b, a in agg.items()},
        "generated_utc": __import__("time").strftime("%Y-%m-%dT%H:%M:%SZ",
                                                     __import__("time").gmtime()),
    }
    with open(os.path.join(args.out_dir, "screen_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"hit accessions: {len(hits)}; negative violations: {len(neg_viol)}; "
          f"distinct bioprojects among hits: {independence['n_distinct_bioprojects']}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
