#!/usr/bin/env python3
"""Controls gate for the local E. coli screen (RUN_PLAN §7).

Gate 1 (positive): ECBAIT_POSCON_LAMBDA_01 vs full NC_001416.1 (lambda) and
  ECBAIT_POSCON_HK97_02 vs full NC_002167.1 (HK97) must each reach
  kmer_coverage >= 0.7 AND minimap2 identity >= 0.7 on the matched sequence.
Gate 2 (negative): ECBAIT_SHUF_CONTROL_{01,02} must have kmer_coverage
  < 0.7 vs both control references (and vs every screened accession --
  checked by summarize_screen.py over the tier results).
Host controls (ECBAIT_HOST_CONTROL_{01,02}) are diagnostics only: E. coli
  host-derived baits are EXPECTED to hit E. coli accessions; excluded from
  gates and classification (RUN_PLAN §7.3).

Exit 0 = PASS, 1 = FAIL (writes out-dir/controls_gate.json either way).
Generalized from ntm local_screen/controls_gate.py: bait ids parameterized,
host-control rows appended as diagnostic (pass=None).
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from screen_accessions import canon, load_baits  # noqa: E402

LAMBDA_BAIT = "ECBAIT_POSCON_LAMBDA_01"
HK97_BAIT = "ECBAIT_POSCON_HK97_02"
NEGS = ["ECBAIT_SHUF_CONTROL_01", "ECBAIT_SHUF_CONTROL_02"]
HOSTS = ["ECBAIT_HOST_CONTROL_01", "ECBAIT_HOST_CONTROL_02"]


def b2s_kmer_cov(b2s, panel_fa, target_fa, sets, index, scratch, tag):
    out = os.path.join(scratch, f"ctrl_{tag}.kmers.txt")
    subprocess.run([b2s, "--in-kmers", panel_fa, "--in-sequences", target_fa,
                    "--out-kmers", out, "--output-kmer-positions"],
                   check=True, capture_output=True)
    found = {}
    for line in open(out):
        if "(" in line:
            k = canon(line.split()[0])
            for b in index.get(k, []):
                found[b] = found.get(b, 0) + 1
    os.remove(out)
    return {b: round(found.get(b, 0) / len(ks), 4) for b, ks in sets.items()}


def minimap2_identity(minimap2, target_fa, bait_fa, preset_ladder=("sr", "asm5", "asm20")):
    for preset in preset_ladder:
        r = subprocess.run([minimap2, "-c", "--eqx", "-x", preset, target_fa, bait_fa],
                           capture_output=True, text=True)
        if r.returncode != 0:
            continue
        qspan = mlen = matches = 0
        qlen = 0
        for line in open(bait_fa):
            if not line.startswith(">"):
                qlen += len(line.strip())
        for line in r.stdout.splitlines():
            if line.startswith("@") or not line.strip():
                continue
            f = line.split("\t")
            cg = [x for x in f if x.startswith("cg:Z:")]
            if len(f) < 12 or not cg:
                continue
            matches += 1
            for n, op in re.findall(r"(\d+)([MIDNSH=X])", cg[0][5:]):
                if op in "M=X":
                    qspan += int(n)
                    if op == "X":
                        mlen += int(n)
                elif op in "ID":
                    mlen += int(n)
        if matches:
            return {"preset": preset, "qcov": round(qspan / qlen, 4) if qlen else None,
                    "identity": round(1 - mlen / qspan, 4) if qspan else None,
                    "n_matches": matches}
    return {"preset": None, "qcov": None, "identity": None, "n_matches": 0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baits", required=True)
    ap.add_argument("--lambda-fa", required=True)
    ap.add_argument("--hk97-fa", required=True)
    ap.add_argument("--b2s", required=True)
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    scratch = os.path.join(args.out_dir, "scratch")
    os.makedirs(scratch, exist_ok=True)

    sets, index = load_baits(args.baits)
    results = {"generated_utc": __import__("time").strftime("%Y-%m-%dT%H:%M:%SZ",
               __import__("time").gmtime()), "rows": [], "gates": {}}

    def bait_seq(bid):
        seq, name, buf = None, None, []
        for line in open(args.baits):
            if line.startswith(">"):
                if name == bid:
                    seq = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if name == bid:
            seq = "".join(buf)
        return seq

    for bait, ref_fa, tag in ((LAMBDA_BAIT, args.lambda_fa, "LAMBDA"),
                              (HK97_BAIT, args.hk97_fa, "HK97")):
        cov = b2s_kmer_cov(args.b2s, args.baits, ref_fa, sets, index, scratch, tag)
        bait_fa = os.path.join(scratch, f"{bait}.fa")
        with open(bait_fa, "w") as fh:
            fh.write(f">{bait}\n{bait_seq(bait)}\n")
        aln = minimap2_identity(args.minimap2, ref_fa, bait_fa)
        kcov = cov.get(bait)
        row = {"bait": bait, "reference": os.path.basename(ref_fa),
               "kmer_coverage": kcov, "minimap2": aln,
               "pass": bool(kcov is not None and kcov >= 0.7
                            and (aln.get("identity") or 0) >= 0.7)}
        results["rows"].append(row)

    # negatives (gate) + hosts (diagnostic) vs both references
    for ref_fa, tag in ((args.lambda_fa, "LAMBDA"), (args.hk97_fa, "HK97")):
        cov = b2s_kmer_cov(args.b2s, args.baits, ref_fa, sets, index, scratch, "neg_" + tag)
        for nb in NEGS:
            results["rows"].append({"bait": nb, "reference": os.path.basename(ref_fa),
                                    "kmer_coverage": cov.get(nb), "minimap2": None,
                                    "pass": bool((cov.get(nb) or 0) < 0.7)})
        for hb in HOSTS:
            results["rows"].append({"bait": hb, "reference": os.path.basename(ref_fa),
                                    "kmer_coverage": cov.get(hb), "minimap2": None,
                                    "pass": None,
                                    "diagnostic": "host_control_expected_host_like"})

    pos_ok = all(r["pass"] for r in results["rows"] if r["bait"] in (LAMBDA_BAIT, HK97_BAIT))
    neg_ok = all(r["pass"] for r in results["rows"] if r["bait"] in NEGS)
    results["gates"] = {"positive_controls": "PASS" if pos_ok else "FAIL",
                        "negative_controls_refs": "PASS" if neg_ok else "FAIL",
                        "overall": "PASS" if (pos_ok and neg_ok) else "FAIL"}
    with open(os.path.join(args.out_dir, "controls_gate.json"), "w") as fh:
        json.dump(results, fh, indent=1)
    print(json.dumps(results["gates"]))
    for r in results["rows"]:
        print(f"  {r['bait']:26s} vs {r['reference']:14s} kcov={r['kmer_coverage']} "
              f"aln={r['minimap2']} pass={r['pass']}")
    return 0 if results["gates"]["overall"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
