#!/usr/bin/env python3
"""Build the prophage -> host accession mapping for the full 132k prophage cohort.

Each prophage leaf in research/mash_tree/full_prophages_tree.nwk is named
  GCF_xxx_prophage_N
The host accession is the prefix before '_prophage_'.  It must exist in the
host tree (research/host_mash_tree/host_tree.nwk, 26,074 leaves).

Outputs (under research/cophylogeny/mapping/):
  prophage_host_map.tsv      - full 132,393-row mapping
  mapping_stats.json         - coverage / multiplicity summary
"""
import csv, json, re, sys, os

RES = "research"

def load_tree_leaves(nwk_path):
    txt = open(nwk_path).read()
    labels = set(re.findall(r'([A-Za-z][A-Za-z0-9._]+?)(?=[:,)])', txt))
    return labels

def host_of(prophage_id):
    return prophage_id.split("_prophage_")[0]

def main():
    host_tree = os.path.join(RES, "host_mash_tree", "host_tree.nwk")
    p_tree    = os.path.join(RES, "mash_tree",   "full_prophages_tree.nwk")
    labels    = os.path.join(RES, "mash_tree",   "full_prophages_labels.csv")

    hosts = load_tree_leaves(host_tree)
    pro  = load_tree_leaves(p_tree)
    pro  = {p for p in pro if "_prophage_" in p}
    print("host leaves:", len(hosts), "prophage leaves:", len(pro))

    # authoritative ordered list of prophage leaves from labels csv (same ids)
    order = []
    lab = csv.DictReader(open(labels))
    lab_ids = [r["sequence"] for r in lab]
    # use labels csv ordering; should equal tree leaf set
    order = lab_ids
    assert set(order) == pro, "labels csv ids != tree leaf ids: %d vs %d" % (len(set(order)), len(pro))

    out = os.path.join(RES, "cophylogeny", "mapping", "prophage_host_map.tsv")
    with open(out, "w") as f:
        f.write("prophage_id\thost_accession\thost_in_tree\n")
        missing = []
        mapped  = {}
        for p in order:
            h = host_of(p)
            in_tree = "yes" if h in hosts else "no"
            f.write("%s\t%s\t%s\n" % (p, h, in_tree))
            if in_tree == "no":
                missing.append(p)
            mapped[h] = mapped.get(h, 0) + 1

    n_covered = len(order) - len(missing)
    stats = {
        "n_prophage_leaves": len(order),
        "n_mapping_covered": n_covered,
        "coverage_pct": round(100.0 * n_covered / len(order), 4),
        "n_missing_host_in_tree": len(missing),
        "missing_sample": missing[:10],
        "n_unique_hosts_mapped": len(mapped),
        "host_multiplicity_min": min(mapped.values()),
        "host_multiplicity_max": max(mapped.values()),
        "host_multiplicity_median": sorted(mapped.values())[len(mapped)//2],
    }
    with open(os.path.join(RES, "cophylogeny", "mapping", "mapping_stats.json"), "w") as f:
        json.dump(stats, f, indent=2)
    print(json.dumps(stats, indent=2))

if __name__ == "__main__":
    main()
