#!/usr/bin/env python3
"""
compare_v1_v2.py — quantify the v1 -> v2 redo of the NTM ML phage genome
release.

Compares:
  1. prophage calls: v1 geNomad (ntm/v1/prophages/ntm_prophages.csv) vs the
     collaborator's manifest calls restricted to the v2 canonical cohort
     (NTM_QC_passed_prophage_master_manifest.tsv). Prophages are matched by
     accession (numeric, GCA/GCF-insensitive) + contig + reciprocal interval
     overlap >= 50%.
  2. clades: v1 tight clades (913) vs v2 (2388) matched by MEMBER overlap of
     the matched prophages (each v1 clade -> modal v2 clade).
  3. host clades represented in each release (v1 per_ntm_clade_summary vs v2).
  4. ML genome length / GC distributions (from the respective qc_table.tsv).

Outputs (ntm/v2/release/):
  v1_v2_clade_map.tsv    per v1 clade: matched members, modal v2 clade, fraction
  v1_v2_prophage_stats.tsv  per shared accession: n_v1, n_v2, n_matched
  v1_v2_summary.json     all headline numbers
"""
import csv
import json
import os
import statistics
from collections import Counter, defaultdict

V1 = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1"
V2 = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2"
REL = f"{V2}/release"
RECIP = 0.5  # reciprocal overlap fraction of BOTH intervals


def norm_acc(a):
    """GCA_000123.4 / GCF_000123.4 -> 000123 (numeric base)."""
    base = a.split(".")[0]
    return base.split("_", 1)[1] if "_" in base else base


def norm_ctg(c):
    """PanSN `acc#hap#ctig` -> contig; strip NZ_/NC_ prefixes."""
    c = c.split("#")[-1]
    for p in ("NZ_", "NC_"):
        if c.startswith(p):
            c = c[len(p):]
    return c


def load_prophages():
    # v1 geNomad calls
    v1 = defaultdict(list)  # (numeric_acc, ctig) -> [(start, end, prophage_id)]
    with open(f"{V1}/prophages/ntm_prophages.csv") as f:
        for row in csv.DictReader(f):
            v1[(norm_acc(row["genome"]), norm_ctg(row["scaffold"]))].append(
                (int(row["begin"]), int(row["end"]), row["prophage_id"]))
    # v2 collaborator calls (dedup GCA/GCF twins by prophage_id suffix on numeric)
    seen = set()
    v2 = defaultdict(list)
    with open("/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/inputs/"
              "NTM_QC_passed_prophage_master_manifest.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if not row["prophage_id"]:
                continue
            num = norm_acc(row["accession"])
            pid = row["prophage_id"]
            suffix = pid.split("_prophage")[-1]
            key = (num, suffix)
            if key in seen:
                continue
            seen.add(key)
            v2[(num, norm_ctg(row["prophage_contig"]))].append(
                (int(float(row["prophage_start"])), int(float(row["prophage_end"])), pid))
    return v1, v2


def match(v1, v2):
    """Return (v1id->v2id, v2id->count_of_v1_matches)."""
    fwd = {}
    rev = Counter()
    keys = set(v1) & set(v2)
    for k in keys:
        for s1, e1, pid1 in v1[k]:
            best, best_ov = None, 0
            for s2, e2, pid2 in v2[k]:
                ov = min(e1, e2) - max(s1, s2)
                if ov <= 0:
                    continue
                if ov >= RECIP * (e1 - s1) and ov >= RECIP * (e2 - s2):
                    if ov > best_ov:
                        best, best_ov = pid2, ov
            if best:
                fwd[pid1] = best
                rev[best] += 1
    return fwd, rev


def dist_stats(qc_path, status_filter=None):
    rows = [r for r in csv.DictReader(open(qc_path), delimiter="\t")]
    if status_filter:
        rows = [r for r in rows if r["status"] == status_filter]
    lens = sorted(int(r["length_bp"]) for r in rows)
    gc = sorted(float(r["gc_content"]) for r in rows)
    nr = [int(r["max_n_run"]) for r in rows]
    comp = [float(r["completeness"]) for r in rows if r["completeness"] != ""]
    buckets = Counter()
    for L in lens:
        buckets["<10kb" if L < 10_000 else
                "10-50kb" if L < 50_000 else
                "50-100kb" if L < 100_000 else
                "100-150kb" if L < 150_000 else ">=150kb"] += 1
    return {
        "n": len(lens),
        "len_min": lens[0], "len_median": int(statistics.median(lens)),
        "len_mean": int(statistics.mean(lens)), "len_max": lens[-1],
        "len_buckets": dict(buckets),
        "gc_median": round(statistics.median(gc), 4),
        "gc_min": round(gc[0], 4), "gc_max": round(gc[-1], 4),
        "n_run_median": int(statistics.median(nr)),
        "n_ge100": sum(1 for n in nr if n >= 100),
        "n_max": max(nr),
        "comp_median": round(statistics.median(comp), 3) if comp else None,
        "comp_lt05": sum(1 for c in comp if c < 0.5),
        "comp_n": len(comp),
    }


def main():
    os.makedirs(REL, exist_ok=True)
    out = {}

    v1p, v2p = load_prophages()
    fwd, rev = match(v1p, v2p)
    n_v1 = sum(len(v) for v in v1p.values())
    n_v2 = sum(len(v) for v in v2p.values())
    v1_accs = {k[0] for k in v1p}
    v2_accs = {k[0] for k in v2p}
    shared = v1_accs & v2_accs
    n_v1_shared = sum(len(v) for k, v in v1p.items() if k[0] in shared)
    n_v2_shared = sum(len(v) for k, v in v2p.items() if k[0] in shared)
    out["prophages"] = {
        "v1_total": n_v1, "v2_total": n_v2,
        "v1_genomes": len(v1_accs), "v2_genomes": len(v2_accs),
        "shared_genomes": len(shared),
        "v1_on_shared": n_v1_shared, "v2_on_shared": n_v2_shared,
        "matched": len(fwd),
        "v1_unmatched": n_v1_shared - len(fwd),
        "v2_unmatched": n_v2_shared - len(rev),
    }

    # per-shared-accession table
    per_acc = defaultdict(lambda: [0, 0, 0])
    for k, v in v1p.items():
        if k[0] in shared:
            per_acc[k[0]][0] += len(v)
    for k, v in v2p.items():
        if k[0] in shared:
            per_acc[k[0]][1] += len(v)
    for pid1, pid2 in fwd.items():
        for k, v in v1p.items():
            pass
    # recompute matched per accession from fwd via v1 index
    pid2num = {}
    for (num, _), v in v2p.items():
        for _, _, pid in v:
            pid2num[pid] = num
    for pid1 in fwd:
        num = pid2num[fwd[pid1]]
        per_acc[num][2] += 1
    with open(f"{REL}/v1_v2_prophage_stats.tsv", "w") as f:
        f.write("numeric_accession\tn_v1_prophages\tn_v2_prophages\tn_matched\n")
        for num in sorted(per_acc):
            a, b, c = per_acc[num]
            f.write(f"{num}\t{a}\t{b}\t{c}\n")

    # clade matching by member overlap
    tc1 = json.load(open(f"{V1}/mash_clades/clades/0/tight_clades.json"))
    tc2 = json.load(open(f"{V2}/clades/0/tight_clades.json"))
    member2clade2 = {}
    for cid, ms in tc2.items():
        for m in ms:
            member2clade2[m] = cid
    clade_rows = []
    matched_clades = majority = split = vanished = 0
    for cid, members in sorted(tc1.items()):
        targets = Counter()
        for m in members:
            # v1 member id `ACC_prophage_N`; match key must carry accession+contig
            t = fwd.get(m)
            if t and t in member2clade2:
                targets[member2clade2[t]] += 1
        n_matched = sum(targets.values())
        if n_matched == 0:
            vanished += 1
            clade_rows.append((cid, len(members), 0, "NA", 0.0))
            continue
        matched_clades += 1
        modal, cnt = targets.most_common(1)[0]
        frac = cnt / n_matched
        if len(targets) > 1:
            split += 1
        if frac >= 0.5 and n_matched >= 2:
            majority += 1
            tag = "majority"
        elif frac >= 0.5:
            tag = "single_member"
        else:
            tag = "split"
        clade_rows.append((cid, len(members), n_matched, modal, round(frac, 3)))
    with open(f"{REL}/v1_v2_clade_map.tsv", "w") as f:
        f.write("v1_clade\tn_members_v1\tn_members_matched\tmodal_v2_clade\tmodal_fraction\n")
        for r in clade_rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    out["clades"] = {
        "v1_clades": len(tc1), "v2_clades": len(tc2),
        "v1_matched_ge1": matched_clades, "v1_majority_matched": majority,
        "v1_split": split, "v1_vanished": vanished,
    }

    # host clades represented
    def host_clades(path):
        return {r["host_clade_id"] for r in csv.DictReader(open(path), delimiter="\t")}
    hc1 = host_clades(f"{V1}/release/per_ntm_clade_summary.tsv")
    hc2 = host_clades(f"{REL}/per_ntm_clade_summary.tsv")
    out["host_clades"] = {"v1": len(hc1), "v2": len(hc2)}

    # distributions
    out["v1_ml"] = dist_stats(f"{REL}/v1_qc_baseline/qc_table.tsv")
    out["v2_ml"] = dist_stats(f"{REL}/qc_table.tsv")
    out["v2_ml_alignable_only"] = dist_stats(f"{REL}/qc_table.tsv", "ml")

    with open(f"{REL}/v1_v2_summary.json", "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
