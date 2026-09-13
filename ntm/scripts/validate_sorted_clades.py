#!/usr/bin/env python3
"""
validate_sorted_clades.py — deterministic re-validation of the corrected
(sorted-ids) NTM v1/v2 tight-clade re-derivation (task verify-v2-clade-order).

Independently re-checks, from the on-disk artifacts only:
  1. v2 corrected clades: sum of clade sizes == 8,502; every prophage in
     exactly one clade (all unique); every non-singleton internal median
     MASH distance <= 0.25 (recomputed from distances.npz, not trusted from
     the run's JSON); corrected matrix has 0 NaN pairs.
  2. corrected v2 clade count is reported alongside the published number with
     the delta, in ntm/v2/clades_sorted_report.md, and the same numbers as in
     clades_sorted/sorted_vs_published.json.
  3. frozen v1/v2 outputs untouched: every frozen input hashed by the driver
     BEFORE the run still hashes to the same value now, and the driver's own
     before-vs-after gate passed (recorded in sorted_vs_published.json); the
     corrected outputs live only in the new *_sorted outdirs.
  4. the v2 report notes the impact on the downstream artifacts that consumed
     the clade count (ml / ml_validation / partition / release reports).
  5. v1: corrected counts == published counts (913 / 472 / 441) and the
     corrected matrix has 0 NaN.

Prints one PASS/FAIL line per check; exits non-zero on any failure.

Usage: python3 ntm/scripts/validate_sorted_clades.py
"""
import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
NVME = "/mnt/nvme3n1/erikg/phind-genome-work/ntm"
THRESHOLD = 0.25

failures = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}: {name}" + (f" — {detail}" if detail else ""))
    if not ok:
        failures.append(name)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def load(v):
    base = f"{NVME}/{v}"
    tc = json.load(open(f"{base}/clades_sorted/0/tight_clades.json"))
    sim = json.load(open(f"{base}/clades_sorted/0/clade_similarity.json"))["per_clade"]
    d = np.load(f"{base}/clades_sorted/0/distances.npz", allow_pickle=True)
    return tc, sim, d["D"], [str(x) for x in d["members"]], \
        json.load(open(f"{base}/clades_sorted/sorted_vs_published.json"))


def main():
    # ---------- 1. v2 corrected clade invariants ----------------------------
    tc2, sim2, D2, mem2, cmp2 = load("v2")
    n2 = cmp2["n_prophages"]
    sizes = [len(ms) for ms in tc2.values()]
    allm = [m for ms in tc2.values() for m in ms]
    check("v2 sum of clade sizes == 8,502", sum(sizes) == n2 == 8502,
          f"sum={sum(sizes)}")
    check("v2 every prophage in exactly one clade (ids unique)",
          len(set(allm)) == n2, f"unique={len(set(allm))}")
    offd = ~np.eye(n2, dtype=bool)
    check("v2 corrected matrix has 0 NaN pairs",
          int(np.isnan(D2[offd]).sum()) == 0)
    pos = {s: i for i, s in enumerate(mem2)}
    worst, bad = 0.0, []
    for cid, ms in tc2.items():
        if len(ms) < 2:
            continue
        idx = np.array([pos[s] for s in ms])
        vals = D2[np.ix_(idx, idx)][np.triu_indices(len(ms), 1)]
        med = float(np.median(vals))
        worst = max(worst, med)
        if med > THRESHOLD:
            bad.append((cid, med))
    check("v2 every non-singleton median <= 0.25 (recomputed from matrix)",
          not bad, f"worst={worst:.6f}, {len(bad)} violations")

    # ---------- 2. corrected vs published with delta ------------------------
    pub, cor = cmp2["published_fragmented"], cmp2["corrected_sorted"]
    rep = open(os.path.join(REPO, "ntm/v2/clades_sorted_report.md")).read()
    check("v2 corrected count 767 vs published 2388 with delta -1621 recorded",
          cor["n_clades"] == 767 and pub["n_clades"] == 2388
          and cmp2["delta_clades"] == -1621,
          f"{pub['n_clades']} -> {cor['n_clades']} (delta {cmp2['delta_clades']})")
    check("v2 report quotes both counts and the delta",
          "2,388" in rep and "767" in rep and "1,621" in rep)
    check("v2 report alignable/singleton split matches artifacts",
          cor["n_alignable"] == 413 and cor["n_singletons"] == 354
          and f"{cor['n_alignable']}" in rep and f"{cor['n_singletons']}" in rep,
          f"{cor['n_alignable']} alignable + {cor['n_singletons']} singletons")

    # ---------- 3. frozen outputs untouched ---------------------------------
    for v in ("v1", "v2"):
        cmpv = json.load(open(f"{NVME}/{v}/clades_sorted/sorted_vs_published.json"))
        before = cmpv["frozen_input_sha256_before"]
        check(f"{v} driver before-vs-after frozen hash gate recorded",
              cmpv["checks"].get("frozen_inputs_unchanged_before_vs_after") is True
              and cmpv["checks"].get("frozen_inputs_hashed_before") == len(before),
              f"{len(before)} files")
        drift = [p for p, h in before.items() if sha256(p) != h]
        check(f"{v} frozen inputs still hash to pre-run values", not drift,
              f"{len(before)} files re-hashed now"
              + (f"; DRIFTED: {drift}" if drift else ""))
        inside = set(os.listdir(f"{NVME}/{v}/mash_clades")) & \
                 {"ids_sorted.txt", "prophages_mash_sorted.dist"}
        check(f"{v} no corrected outputs inside frozen mash_clades/", not inside)
        check(f"{v} corrected outputs only in new outdirs",
              os.path.isdir(f"{NVME}/{v}/clades_sorted")
              and os.path.isdir(f"{NVME}/{v}/mash_clades_sorted"))

    # ---------- 4. report notes downstream impact ---------------------------
    need = ["ml_report.md", "ml_validation_report.md", "partition_report.md",
            "RELEASE.md", "v1_v2_comparison.md"]
    missing = [f for f in need if f not in rep]
    check("v2 report notes downstream artifacts (ml/ml_validation/partition/release)",
          not missing, f"missing: {missing}" if missing else f"all {len(need)} noted")

    # ---------- 5. v1 --------------------------------------------------------
    tc1, sim1, D1, mem1, cmp1 = load("v1")
    n1 = cmp1["n_prophages"]
    offd1 = ~np.eye(n1, dtype=bool)
    check("v1 corrected matrix has 0 NaN pairs",
          int(np.isnan(D1[offd1]).sum()) == 0)
    sizes1 = [len(ms) for ms in tc1.values()]
    check("v1 sum of clade sizes == 10,438 and ids unique",
          sum(sizes1) == n1 == 10438
          and len({m for ms in tc1.values() for m in ms}) == n1)
    pub1, cor1 = cmp1["published_fragmented"], cmp1["corrected_sorted"]
    check("v1 corrected counts == published counts (913/472/441)",
          cor1["n_clades"] == pub1["n_clades"] == 913
          and cor1["n_alignable"] == pub1["n_alignable"] == 472
          and cor1["n_singletons"] == pub1["n_singletons"] == 441)

    print()
    if failures:
        print(f"VALIDATION FAILED: {len(failures)} failing check(s): {failures}")
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
