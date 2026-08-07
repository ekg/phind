#!/usr/bin/env python3
"""Build NTM host clades from whole-genome mash distances.

Inputs (downstream of download-7-352-ntm-genomes):
  * <BASE>/genomes/canonical_objects/<acc>/<acc>.pansn.fa.gz   (7,303 PanSN bgzip)
  * <BASE>/accessions/ntm_accession_manifest.tsv               (species/genus/organism)

Pipeline (upstream)
  1. mash sketch -k 21 -s 10000 (directly on the bgzip PanSN FASTA)
  2. mash dist all-vs-all -> host.dist

This script (from host.dist):
  3. connected components (single-linkage) at dist <= 0.05 (~95% ANI species-level)
  4. each component = a host clade, labeled by dominant NCBI species (majority vote)

Usage:
  1. Build a list of the PanSN bgzip paths and sketch/dist them:
       find <BASE>/genomes/canonical_objects -name '*.pansn.fa.gz' | sort > <BASE>/host_clades/filelist.txt
       cd <BASE>/host_clades && mash sketch -k 21 -s 10000 -l filelist.txt -o host -p 32
       mash dist -p 32 host.msh host.msh > host.dist
  2. Cluster:
       python3 host_clades_mash.py <BASE>

  * host_clades.tsv        accession | host_clade_id | species | genus | organism
  * host_clade_summary.tsv host_clade_id | count | dominant_species | species_distribution

Run (mash must be on PATH):
    python3 host_clades_mash.py <BASE>
"""
import os
import sys
from collections import defaultdict, Counter

THRESHOLD = 0.05  # mash distance ~ 95% ANI (species-level)


def load_manifest(manifest_path):
    ann = {}
    with open(manifest_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            ann[p[idx["primary_acc"]]] = {
                "species": p[idx["species"]],
                "genus": p[idx["genus"]],
                "organism": p[idx["organism"]],
            }
    return ann


def subspecies_of(organism):
    """Split M. abscessus complex into subspecies; else None (use species field)."""
    o = organism.lower()
    for subsp in ("abscessus", "bolletii", "massiliense"):
        if f"abscessus subsp. {subsp}" in o:
            return f"Mycobacteroides abscessus subsp. {subsp}"
    if "mycobacteroides abscessus" in o or "mycobacterium abscessus" in o:
        if "subsp." not in o and "subsp " not in o:
            return "Mycobacteroides abscessus subsp. abscessus"
    return None


def main(base):
    outdir = os.path.join(base, "host_clades")
    os.makedirs(outdir, exist_ok=True)
    ann = load_manifest(os.path.join(base, "accessions", "ntm_accession_manifest.tsv"))

    # 2. clustering from host.dist (sketch+dist done upstream per docstring)
    # union-find over dist <= THRESHOLD
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    genome_accs = set()
    with open(os.path.join(outdir, "host.dist")) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            a = os.path.basename(parts[0]).replace(".pansn.fa.gz", "")
            b = os.path.basename(parts[1]).replace(".pansn.fa.gz", "")
            genome_accs.add(a)
            genome_accs.add(b)
            if float(parts[2]) <= THRESHOLD:
                parent.setdefault(a, a)
                parent.setdefault(b, b)
                union(a, b)

    # 3. connected components -> clades, size-desc order -> ids
    comps = defaultdict(list)
    for acc in genome_accs:
        comps[find(acc) if acc in parent else acc].append(acc)
    clades = sorted(comps.values(), key=len, reverse=True)
    clade_of = {}
    clade_meta = []
    for i, members in enumerate(clades, 1):
        cid = f"host_clade_{i:04d}"
        clade_meta.append((cid, members))
        for m in members:
            clade_of[m] = cid

    # dominant species label per clade (splitting abscessus subspp)
    clade_counter = {}
    for cid, members in clade_meta:
        labels = []
        for m in members:
            a = ann.get(m, {})
            labels.append(subspecies_of(a.get("organism", "")) or a.get("species", "unknown"))
        clade_counter[cid] = Counter(labels)

    # host_clades.tsv
    rows = sorted(
        (m, clade_of[m], ann.get(m, {}).get("species", ""),
         ann.get(m, {}).get("genus", ""), ann.get(m, {}).get("organism", ""))
        for m in genome_accs
    )
    with open(os.path.join(outdir, "host_clades.tsv"), "w") as f:
        f.write("accession\thost_clade_id\tspecies\tgenus\torganism\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    # host_clade_summary.tsv
    with open(os.path.join(outdir, "host_clade_summary.tsv"), "w") as f:
        f.write("host_clade_id\tcount\tdominant_species\tspecies_distribution\n")
        for cid, members in clade_meta:
            dist = ";".join(f"{s}:{n}" for s, n in clade_counter[cid].most_common())
            f.write(f"{cid}\t{len(members)}\t{clade_counter[cid].most_common(1)[0][0]}\t{dist}\n")

    print(f"clades: {len(clade_meta)}; assigned: {len(genome_accs)} accessions")
    print(f"wrote {os.path.join(outdir, 'host_clades.tsv')} and host_clade_summary.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    main(os.path.abspath(sys.argv[1]))