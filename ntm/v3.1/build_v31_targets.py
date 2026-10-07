#!/usr/bin/env python3
"""Build v3.1 per-clade genome-length targets from CheckV member AAI estimates.

Input : $NVME/ntm/v3.1/member_checkv/completeness.tsv  (CheckV on full_prophages.fa)
        $NVME/ntm/v3/clades/0/tight_clades.json        (clade_id -> member prophage ids)
        $NVME/ntm/v3/clades/<cid>/sequences.fa         (fallback: member sequence lengths)
Output: $NVME/ntm/v3.1/targets.tsv  (clade_id  target_len  source  n_aai  n_members)

target_len = median CheckV aai_expected_length over the clade's members that have a
usable estimate (aai_confidence in {high, medium} and a non-null expected length).
Where fewer than 2 usable estimates exist, fall back to the median member SEQUENCE
length (source = member_seq_median) -- which we measured to run ~0.78x low, so those
clades are flagged for review rather than silently trusted.
"""
from __future__ import annotations
import csv, json, os, statistics, sys

NV = "/mnt/nvme3n1/erikg/phind-genome-work"
CKV = f"{NV}/ntm/v3.1/member_checkv/completeness.tsv"
CLADES = f"{NV}/ntm/v3/clades/0/tight_clades.json"
OUT = f"{NV}/ntm/v3.1/targets.tsv"


def member_seq_lengths(cid):
    p = f"{NV}/ntm/v3/clades/{cid}/sequences.fa"
    L, name, n = {}, None, 0
    if not os.path.exists(p):
        return L
    for line in open(p):
        if line.startswith(">"):
            if name is not None:
                L[name] = n
            name, n = line[1:].strip().split()[0], 0
        else:
            n += len(line.strip())
    if name is not None:
        L[name] = n
    return L


def main():
    aai = {}
    with open(CKV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                exp = float(r.get("aai_expected_length", ""))
            except Exception:
                exp = None
            conf = (r.get("aai_confidence") or "").strip().lower()
            aai[r["contig_id"].strip()] = (exp, conf)
    clades = json.load(open(CLADES))
    rows = []
    n_target = n_fallback = 0
    for cid in sorted(clades):
        members = clades[cid]
        est = []
        for m in members:
            key = m.strip()
            exp, conf = aai.get(key, (None, ""))
            # tolerate PanSN-formatted keys in the CheckV output
            if exp is None:
                for cand in (key.split("#")[-1], key.replace("#", "_")):
                    if cand in aai and aai[cand][0] is not None:
                        exp, conf = aai[cand]
                        break
            if exp and conf in ("high", "medium") and 5000 <= exp <= 300000:
                est.append(exp)
        if len(est) >= 2:
            rows.append((cid, int(statistics.median(est)), "checkv_aai_median",
                         len(est), len(members)))
            n_target += 1
        else:
            L = member_seq_lengths(cid)
            if L:
                rows.append((cid, int(statistics.median(L.values())),
                             "member_seq_median", len(est), len(members)))
                n_fallback += 1
            else:
                rows.append((cid, 0, "none", 0, len(members)))
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["clade_id", "target_len", "source", "n_aai", "n_members"])
        w.writerows(rows)
    print(f"wrote {OUT}: checkv_aai_median={n_target} member_seq_median_fallback={n_fallback}")
    fb = [r for r in rows if r[2] == "member_seq_median"]
    if fb:
        print(f"NOTE: {len(fb)} clades used the (biased-low) member median fallback: "
              f"{[r[0] for r in fb[:12]]}")


if __name__ == "__main__":
    main()
