#!/usr/bin/env python3
"""Prepare clade + host-group data for the cophylogeny analyses.

Reads:
  research/clades/<community>/tight_clades.json   -> clade_id -> [prophage ids]
  research/mash_tree/full_prophages_labels.csv     -> prophage -> community
  research/phylogroups/phylogroups.tsv             -> host accession -> phylogroup

Writes (research/cophylogeny/clade_resolution/):
  clade_memberships.tsv      prophage_id -> clade_id (19638 rows within tight clades)
  clade_meta.tsv             clade_id -> community, n_members, n_hosts
  host_phylogroup_map.tsv    host accession -> phylogroup (used everywhere)
  association_matrix.tsv     prophage_clade x host_phylogroup counts
  host_clade_matrix.tsv      host accession x prophage-clade presence (0/1)
  clade_host_set.tsv         clade_id -> comma-separated host accessions
"""
import csv, json, glob, os
from collections import defaultdict, Counter

RES = "research"
OUT = os.path.join(RES, "cophylogeny", "clade_resolution")
os.makedirs(OUT, exist_ok=True)

def host_of(p): return p.split("_prophage_")[0]

def main():
    # ---- clade memberships ----
    clade_to_pro = {}
    clade_comm   = {}
    for d in sorted(glob.glob(os.path.join(RES, "clades", "*"))):
        if not os.path.isdir(d):
            continue
        comm = os.path.basename(d)
        tj = os.path.join(d, "tight_clades.json")
        if not os.path.exists(tj):
            continue
        for cid, members in json.load(open(tj)).items():
            clade_to_pro[cid] = members
            clade_comm[cid] = comm

    print("total tight clades:", len(clade_to_pro))

    # ---- host phylogroups ----
    phylo = {}
    with open(os.path.join(RES, "phylogroups", "phylogroups.tsv")) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            phylo[r["accession"]] = r["phylogroup"]

    # ---- write clade_memberships + clade_meta ----
    with open(os.path.join(OUT, "clade_memberships.tsv"), "w") as f:
        f.write("prophage_id\tclade_id\tcommunity\thost_accession\tphylogroup\n")
        for cid in sorted(clade_to_pro):
            for p in clade_to_pro[cid]:
                h = host_of(p)
                f.write("%s\t%s\t%s\t%s\t%s\n" % (p, cid, clade_comm[cid], h, phylo.get(h, "Unknown")))

    with open(os.path.join(OUT, "clade_meta.tsv"), "w") as f:
        f.write("clade_id\tcommunity\tn_members\tn_hosts\thost_list\n")
        for cid in sorted(clade_to_pro):
            members = clade_to_pro[cid]
            hosts = sorted({host_of(p) for p in members})
            f.write("%s\t%s\t%d\t%d\t%s\n" % (cid, clade_comm[cid], len(members), len(hosts), ",".join(hosts)))

    with open(os.path.join(OUT, "host_phylogroup_map.tsv"), "w") as f:
        f.write("host_accession\tphylogroup\n")
        for h in sorted(phylo):
            f.write("%s\t%s\n" % (h, phylo[h]))

    # ---- clade x phylogroup counts (association matrix) ----
    # normalize phylogroup set: map the rare/ambiguous to a coarse category
    def coarse(pg):
        if pg in ("A", "B1", "B2", "D", "E", "F", "C", "G"):
            return pg
        if pg in ("cladeI", "cladeII", "cladeIII", "cladeIV", "cladeV"):
            return "clade"
        if pg in ("E or cladeI",):
            return "E"
        return "Other"   # Unknown, Non Escherichia, albertii, H

    counts = Counter()
    rows = []
    with open(os.path.join(OUT, "clade_memberships.tsv")) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            counts[(r["clade_id"], coarse(r["phylogroup"]))] += 1
            rows.append(r)

    phylo_groups = ["A", "B1", "B2", "C", "D", "E", "F", "G", "clade", "Other"]
    clades = sorted(clade_to_pro)
    with open(os.path.join(OUT, "association_matrix.tsv"), "w") as f:
        f.write("clade_id\t" + "\t".join(phylo_groups) + "\ttotal\n")
        for c in clades:
            tot = 0
            cells = []
            for pg in phylo_groups:
                v = counts.get((c, pg), 0)
                tot += v
                cells.append(str(v))
            cells.append(str(tot))
            f.write(c + "\t" + "\t".join(cells) + "\n")

    # ---- host x clade binary matrix (full 26k) ----
    host_clades = defaultdict(set)
    for r in rows:
        host_clades[r["host_accession"]].add(r["clade_id"])
    hosts = sorted(host_clades)
    with open(os.path.join(OUT, "host_clade_matrix.tsv"), "w") as f:
        f.write("host_accession\t" + "\t".join(clades) + "\tn_clades\n")
        for h in hosts:
            cs = host_clades[h]
            f.write(h + "\t" + "\t".join("1" if c in cs else "0" for c in clades) + "\t%d\n" % len(cs))

    # ---- clade host-set (for Jaccard between clades) ----
    with open(os.path.join(OUT, "clade_host_set.tsv"), "w") as f:
        f.write("clade_id\thosts\n")
        for cid in sorted(clade_to_pro):
            hs = sorted({host_of(p) for p in clade_to_pro[cid]})
            f.write("%s\t%s\n" % (cid, ",".join(hs)))

    summary = {
        "n_tight_clades": len(clades),
        "n_prophages_in_clades": len(rows),
        "n_unique_hosts_with_clades": len(hosts),
        "n_phylogroups": len(phylo_groups),
        "phylogroups": phylo_groups,
    }
    with open(os.path.join(OUT, "prep_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
