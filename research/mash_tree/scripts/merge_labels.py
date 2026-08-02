#!/usr/bin/env python3
"""merge_labels.py — overlay existing cluster/community/MDS labels on tree leaves.

Inputs (existing, from the prophage homology survey):
  /home/erikg/phind/prophage_homology_survey/full_heatmap_clusters.csv
      sequence,community            (all 132,393 prophages, 12 communities)
  /home/erikg/phind/prophage_homology_survey/full_prophage_clusters.csv
      sequence,cluster,genome       (5,001 prophages, 428 clusters)
  /home/erikg/phind/prophage_homology_survey/full_prophage_mds_coords.csv
      sequence,MDS1,MDS2,genome     (5,001 prophages)
This task's:
  data/ids.txt                      (132,393 leaf ids in tree order)

Outputs:
  full_prophages_labels.csv         one row per leaf: sequence,community,
                                    cluster,genome,MDS1,MDS2,in_tree
  full_prophages_tree_labeled.nwk   Newick with leaf names
                                    GCF_..._prophage_N|C<community> (community
                                    from full_heatmap_clusters; "_" if none)
  label_merge_stats.json            merge rates
"""
import collections, csv, json, os, sys

SURVEY = "/home/erikg/phind/prophage_homology_survey"

def read_csv(path, cols):
    rows = {}
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            rows[row["sequence"]] = {c: row.get(c) for c in cols}
    return rows

def main():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ids = [l.strip() for l in open(os.path.join(base, "data", "ids.txt")) if l.strip()]
    n = len(ids)

    community = read_csv(os.path.join(SURVEY, "full_heatmap_clusters.csv"), ["community"])
    clusters  = read_csv(os.path.join(SURVEY, "full_prophage_clusters.csv"),
                         ["cluster", "genome"])
    mds       = read_csv(os.path.join(SURVEY, "full_prophage_mds_coords.csv"),
                         ["MDS1", "MDS2", "genome"])

    n_comm = n_clust = n_mds = 0
    out_rows = []
    # real (non-singleton) community ids
    real_comm = set()
    cc = collections.Counter(community[s]["community"] for s in ids if s in community)
    real_comm = {cid for cid, cnt in cc.items() if cnt > 1}
    print(f"[labels] communities: {len(cc)} total, {len(real_comm)} shared (>1 member), "
          f"{len(cc) - len(real_comm)} singletons")
    with open(os.path.join(base, "full_prophages_labels.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sequence", "community", "cluster", "genome", "MDS1", "MDS2", "in_tree"])
        for seq in ids:
            c = community.get(seq, {}).get("community", "")
            cl = clusters.get(seq, {}).get("cluster", "")
            g1 = clusters.get(seq, {}).get("genome", "")
            m = mds.get(seq, {})
            g2 = m.get("genome", "")
            genome = g1 or g2 or ""
            row = [seq, c, cl, genome,
                   m.get("MDS1", ""), m.get("MDS2", ""), "1"]
            w.writerow(row)
            out_rows.append(row)
            if c != "":
                n_comm += 1
            if cl != "":
                n_clust += 1
            if m.get("MDS1") not in (None, ""):
                n_mds += 1
    print(f"[labels] {n} leaves; community={n_comm} cluster={n_clust} mds={n_mds}")

    # labeled Newick: append |C<community> (single-pass regex, avoids O(n²)
    # string.replace scanning)
    import re
    nwk_path = os.path.join(base, "full_prophages_tree.nwk")
    with open(nwk_path) as f:
        nwk = f.read().strip()
    pat = re.compile(r"(?<=[,:(])GCF_\d+\.\d+_prophage_\d+(?=[:,)])")
    def repl(m):
        s = m.group(0)
        c = community.get(s, {}).get("community", "")
        if c == "":
            return f"{s}|C_"
        if c in real_comm:
            return f"{s}|C{c}"
        return f"{s}|iso"
    labeled = pat.sub(repl, nwk)
    labeled_path = os.path.join(base, "full_prophages_tree_labeled.nwk")
    with open(labeled_path, "w") as f:
        f.write(labeled + "\n")
    print(f"[labels] wrote {labeled_path}")

    stats = {
        "n_leaves": n,
        "with_community": n_comm,
        "community_merge_rate": round(n_comm / n, 4),
        "with_cluster": n_clust,
        "cluster_merge_rate": round(n_clust / n, 4),
        "with_mds": n_mds,
        "mds_merge_rate": round(n_mds / n, 4),
        "total_communities": len(cc),
        "shared_communities": len(real_comm),
        "singleton_communities": len(cc) - len(real_comm),
        "leaves_in_shared_communities": sum(1 for r in out_rows if r[1] in real_comm),
    }
    with open(os.path.join(base, "label_merge_stats.json"), "w") as f:
        json.dump(stats, f, indent=2)
    print(json.dumps(stats, indent=2))

if __name__ == "__main__":
    main()
