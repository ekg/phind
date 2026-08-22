#!/usr/bin/env python3
"""Known-phage calibration (RUN_PLAN §10.6): which confirmed hits match
lambda / HK97 / P2-class references vs none (novelty estimate).

For every confirmed (accession, bait) pair in the confirm JSONL:
  align the bait to a known-phage reference set with the minimap2 ladder
  (sr -> asm5 -> asm20); known-phage-matched iff qcov >= 0.5 AND
  identity >= 0.7 on the matched reference.

References: lambda NC_001416.1, HK97 NC_002167.1 (panel public_refs),
P2 NC_001895.1 (fetched from NCBI EFetch, sha256-logged here).
Poscon baits are excluded (they are slices of the references by design).

Outputs results/known_phage_calibration.tsv + json summary.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import urllib.request

REF_URLS = {
    "P2_NC_001895.1": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                      "?db=nuccore&id=NC_001895.1&rettype=fasta&retmode=text",
}


def http_get(url, retries=3, backoff=(5, 20, 60)):
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "local-ecoli-logan/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            if attempt < retries:
                time.sleep(backoff[min(attempt, len(backoff) - 1)])
    raise RuntimeError(f"GET failed after {retries} retries: {url}: {last}")


def fetch_p2(out_dir, log):
    fa = os.path.join(out_dir, "NC_001895.1.fa")
    if os.path.exists(fa):
        return fa
    body = http_get(REF_URLS["P2_NC_001895.1"])
    with open(fa, "wb") as fh:
        fh.write(body)
    with open(log, "a") as lg:
        lg.write(json.dumps({"event": "efetch_p2", "url": REF_URLS["P2_NC_001895.1"],
                             "size": len(body),
                             "sha256": hashlib.sha256(body).hexdigest(),
                             "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}) + "\n")
    return fa


def minimap2_best(minimap2, target_fa, query_fa, ladder=("sr", "asm5", "asm20")):
    for preset in ladder:
        r = subprocess.run([minimap2, "-c", "--eqx", "-x", preset, target_fa, query_fa],
                           capture_output=True, text=True)
        if r.returncode != 0:
            continue
        qspan = mlen = matches = 0
        qlen = 0
        for line in open(query_fa):
            if not line.startswith(">"):
                qlen += len(line.strip())
        best = None
        for line in r.stdout.splitlines():
            if line.startswith("@") or not line.strip():
                continue
            f = line.split("\t")
            cg = [x for x in f if x.startswith("cg:Z:")]
            if len(f) < 12 or not cg:
                continue
            matches += 1
            m, mm = 0, 0
            for n, op in re.findall(r"(\d+)([MIDNSH=X])", cg[0][5:]):
                if op in "M=X":
                    m += int(n)
                    if op == "X":
                        mm += int(n)
                elif op in "ID":
                    mm += int(n)
            qspan += m
            mlen += mm
            if best is None or m > best[0]:
                best = (m, round(1 - mm / m, 4) if m else None,
                        int(f[7]), int(f[8]), f[5])
        if matches:
            return {"preset": preset, "qcov": round(qspan / qlen, 4) if qlen else 0.0,
                    "identity": round(1 - mlen / qspan, 4) if qspan else None,
                    "best_target": best[4] if best else None,
                    "best_span": f"{best[2]}-{best[3]}" if best else None}
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--confirm-jsonl", required=True)
    ap.add_argument("--baits", required=True)
    ap.add_argument("--lambda-fa", required=True)
    ap.add_argument("--hk97-fa", required=True)
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--log", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    refs = {"lambda_NC_001416.1": args.lambda_fa,
            "HK97_NC_002167.1": args.hk97_fa,
            "P2_NC_001895.1": fetch_p2(args.out_dir, args.log)}
    for name, p in refs.items():
        print(f"ref {name}: {p} sha256="
              f"{hashlib.sha256(open(p,'rb').read()).hexdigest()[:16]}…")

    baits = {}
    name, buf = None, []
    for line in open(args.baits):
        if line.startswith(">"):
            if name:
                baits[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        baits[name] = "".join(buf)

    scratch = os.path.join(args.out_dir, "scratch")
    os.makedirs(scratch, exist_ok=True)
    cache = {}

    def bait_vs_refs(b):
        if b in cache:
            return cache[b]
        res = {}
        bait_fa = os.path.join(scratch, f"{b}.fa")
        with open(bait_fa, "w") as fh:
            fh.write(f">{b}\n{baits[b]}\n")
        for rname, rfa in refs.items():
            aln = minimap2_best(args.minimap2, rfa, bait_fa)
            if aln and (aln["qcov"] or 0) >= 0.5 and (aln["identity"] or 0) >= 0.7:
                res[rname] = aln
        cache[b] = res
        return res

    rows, out = [], open(os.path.join(args.out_dir, "known_phage_calibration.tsv"),
                         "w", newline="")
    import csv as _csv
    w = _csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "tier", "bioproject", "bait", "bait_class",
                "kmer_cov", "known_phage_match", "match_ref", "qcov", "identity",
                "novelty"])
    n_pairs = n_known = n_novel = 0
    for line in open(args.confirm_jsonl):
        rec = json.loads(line)
        acc, tier = rec["accession"], rec.get("tier", "")
        proj = rec.get("bioproject", "")
        for b, info in rec["baits"].items():
            if b.startswith(("ECBAIT_POSCON", "ECBAIT_SHUF", "ECBAIT_HOST")):
                continue
            m = bait_vs_refs(b)
            bclass = ("junction" if "JUNCTION" in b else
                      "member_interior" if "MEMBER" in b else "interior_module")
            n_pairs += 1
            if m:
                n_known += 1
                for rname, aln in m.items():
                    w.writerow([acc, tier, proj, b, bclass, info["kmer_cov"],
                                "yes", rname, aln["qcov"], aln["identity"],
                                "known-phage-matched"])
            else:
                n_novel += 1
                w.writerow([acc, tier, proj, b, bclass, info["kmer_cov"],
                            "no", "-", "-", "-", "novel (no known-phage match)"])
    out.close()
    summary = {
        "n_confirmed_pairs": n_pairs, "n_known_phage_matched": n_known,
        "n_novel_no_known_match": n_novel,
        "novelty_frac": round(n_novel / n_pairs, 4) if n_pairs else None,
        "rule": "known-phage-matched iff bait aligns to lambda/HK97/P2 with "
                "qcov >= 0.5 AND identity >= 0.7 (RUN_PLAN §10.6)",
        "references": {k: os.path.basename(v) for k, v in refs.items()},
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    with open(os.path.join(args.out_dir, "known_phage_calibration.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
