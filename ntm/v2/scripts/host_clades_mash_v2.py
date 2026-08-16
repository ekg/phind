#!/usr/bin/env python3
"""NTM v2 — host clades on the QC-passed cohort (NCBI subset; run-assemblies extendable).

Method (chat-agent update 2026-08-16 — supersedes subsetting v1 host dists):
  full re-sketch over the v2 cohort: mash sketch every genome in
  canonical_objects/ (+ run_assemblies/ when present) -> mash dist all-vs-all
  -> connected components (single-linkage) at dist <= 0.05 (~95% ANI),
  same parameters/logic as v1 ntm/scripts/host_clades_mash.py, fresh clade
  ids for v2.

Annotations: species from the QC-passed accession list
  (NTM_QC_passed_accession_list.tsv, cols genome_id/accession/species/data_source);
  genus = first token of species; organism = species (the v2 accession list
  has no separate organism column). Run-assembly rows (data_source=ASSEMBLY)
  annotate identically once their genomes are ingested.

Pipeline (upstream, must be run first — sketch is the long pole):
  1. find <BASE>/genomes/canonical_objects -name '*.pansn.fa.gz' | sort \\
         > <BASE>/host_clades/filelist.txt
  2. cd <BASE>/host_clades && mash sketch -k 21 -s 10000 -l filelist.txt -o host -p 32
  3. mash dist -p 32 host.msh host.msh > host.dist

This script (from host.dist):
  * clusters at dist <= THRESHOLD (union-find connected components)
  * writes host_clades.tsv (accession | host_clade_id | species | genus | organism)
    and host_clade_summary.tsv (host_clade_id | count | dominant_species | species_distribution)
  * prints summary: clade count, top species, size distribution, TB-complex /
    abscessus-complex composition
  * validates: every filelist accession assigned exactly one clade; output
    rows == filelist accessions; no silent drops

Run-assembly extension: when run_assemblies/ contains *.pansn.fa.gz, append
them to the filelist, re-sketch/dist, and re-run this script unchanged
(schema and clade ids stay fresh). Annotations resolve automatically via the
accession list (run ids included).

Usage:
    python3 host_clades_mash_v2.py <BASE>

(mash must be on PATH for the upstream sketch/dist steps.)
"""
import argparse
import os
import sys
from collections import Counter, defaultdict

THRESHOLD = 0.05  # mash distance ~ 95% ANI (species-level)


def load_accession_annotations(accession_list_path):
    """acc -> dict(species, genus, organism) from the QC-passed accession list.

    First row per accession wins (list is already deduplicated by the
    collaborator; GCA/GCF twins share the same numeric and the canonical
    objects dir name is the preferred accession).
    """
    ann = {}
    with open(accession_list_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            acc = p[idx["accession"]]
            if acc in ann:
                continue
            species = p[idx["species"]].strip()
            genus = species.split(" ", 1)[0] if species else ""
            ann[acc] = {"species": species, "genus": genus, "organism": species}
    return ann


def load_genome_filelist(filelist_path):
    """Set of accession ids from the filelist (basename minus .pansn.fa.gz)."""
    accs = set()
    with open(filelist_path) as f:
        for line in f:
            name = line.strip()
            if not name:
                continue
            accs.add(os.path.basename(name).replace(".pansn.fa.gz", ""))
    return accs


def run_clustering(base, outdir, ann, genome_accs):
    """Union-find connected components over host.dist at <= THRESHOLD (v1 logic)."""
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

    dist_path = os.path.join(outdir, "host.dist")
    n_lines = 0
    n_unions = 0
    with open(dist_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            n_lines += 1
            a = os.path.basename(parts[0]).replace(".pansn.fa.gz", "")
            b = os.path.basename(parts[1]).replace(".pansn.fa.gz", "")
            if a == b or float(parts[2]) > THRESHOLD:
                continue
            parent.setdefault(a, a)
            parent.setdefault(b, b)
            union(a, b)
            n_unions += 1

    # connected components -> clades, size-desc order -> ids
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

    # dominant species label per clade (species strings verbatim from accession list)
    clade_counter = {}
    for cid, members in clade_meta:
        labels = [ann.get(m, {}).get("species") or "unknown" for m in members]
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

    return clade_meta, clade_counter, n_lines, n_unions


def validate(base, outdir, genome_accs, clade_meta, clade_counter, ann):
    """Validation checks; returns (ok, message list)."""
    msgs = []
    ok = True
    _ = base  # annotations already resolved by caller (ann)

    # 1. every filelist accession has exactly one row
    with open(os.path.join(outdir, "host_clades.tsv")) as f:
        f.readline()  # header
        rows = [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]
    out_accs = {r[0] for r in rows}
    if len(rows) != len(genome_accs):
        ok = False
        msgs.append(f"FAIL row count {len(rows)} != filelist {len(genome_accs)}")
    if out_accs != genome_accs:
        ok = False
        msgs.append(f"FAIL accession mismatch: missing={len(genome_accs - out_accs)} extra={len(out_accs - genome_accs)}")
    if len({r[0] for r in rows}) != len(rows):
        ok = False
        msgs.append("FAIL duplicate accession rows")

    # 2. inner-join vs QC-passed accession set (NCBI canonical + runs if present)
    missing_ann = sorted(a for a in genome_accs if a not in ann)
    if missing_ann:
        ok = False
        msgs.append(f"FAIL {len(missing_ann)} accessions missing from accession list: {missing_ann[:5]}")
    else:
        msgs.append("annotation join: all accessions resolve to a species row")

    # 3. clade size distribution sanity: no clade >40% of cohort unless the
    #    dominant species is the abscessus complex (v1 top clade was abscessus)
    #    or the TB complex (kept verbatim per collaborator QC list — ~55% of
    #    the v2 cohort is MTBC-labelled; see chat-agent update 2026-08-16)
    cohort = len(genome_accs)
    top = clade_meta[0]
    frac = len(top[1]) / cohort
    dom = clade_counter[top[0]].most_common(1)[0][0].lower()
    explained = "abscessus" in dom or "tuberculosis" in dom
    msgs.append(f"top clade {top[0]} n={len(top[1])} ({frac:.1%}) dominant={clade_counter[top[0]].most_common(1)[0][0]}")
    if frac > 0.40 and not explained:
        ok = False
        msgs.append(f"FAIL top clade {frac:.1%} > 40% and not abscessus/TB-complex-dominated")
    if len(clade_meta) == 0:
        ok = False
        msgs.append("FAIL no clades produced")
    return ok, msgs


def composition_report(genome_accs, ann):
    """Composition of the cohort (species-level; TB-complex kept verbatim)."""
    counts = Counter(ann.get(a, {}).get("species", "unknown") for a in genome_accs)
    total = len(genome_accs)
    tb = sum(n for s, n in counts.items() if "tuberculosis" in s.lower() and "paratuberculosis" not in s.lower())
    tb_paratuberculosis = sum(n for s, n in counts.items() if "paratuberculosis" in s.lower())
    absc = sum(n for s, n in counts.items() if "abscessus" in s.lower())
    lines = [f"total genomes: {total}"]
    lines.append(f"TB-complex (species containing 'tuberculosis', excl. paratuberculosis): {tb} ({tb/total:.1%}) — kept verbatim per collaborator QC list")
    lines.append(f"paratuberculosis: {tb_paratuberculosis}")
    lines.append(f"abscessus complex (species containing 'abscessus'): {absc} ({absc/total:.1%})")
    lines.append("top 15 species:")
    for s, n in counts.most_common(15):
        lines.append(f"  {n:6d}  {s}")
    return "\n".join(lines)


def main(base):
    outdir = os.path.join(base, "host_clades")
    os.makedirs(outdir, exist_ok=True)

    # genome set: canonical_objects (NCBI) + run_assemblies if present (extension)
    filelist_path = os.path.join(outdir, "filelist.txt")
    if not os.path.exists(filelist_path):
        canon = os.path.join(base, "genomes", "canonical_objects")
        runs = os.path.join(base, "genomes", "run_assemblies")
        paths = []
        for d in (canon, runs):
            if os.path.isdir(d):
                paths.extend(
                    os.path.join(dp, fn)
                    for dp, _, fns in os.walk(d)
                    for fn in fns if fn.endswith(".pansn.fa.gz")
                )
        paths.sort()
        with open(filelist_path, "w") as f:
            for p in paths:
                f.write(p + "\n")
        print(f"wrote {filelist_path}: {len(paths)} genomes")
    genome_accs = load_genome_filelist(filelist_path)
    print(f"filelist genomes: {len(genome_accs)}")

    ann = load_accession_annotations(os.path.join(base, "inputs", "NTM_QC_passed_accession_list.tsv"))

    clade_meta, clade_counter, n_lines, n_unions = run_clustering(base, outdir, ann, genome_accs)

    print(f"host.dist lines: {n_lines}; union ops: {n_unions}")
    print(f"clades: {len(clade_meta)}; assigned: {len(genome_accs)} accessions")
    print(f"wrote {os.path.join(outdir, 'host_clades.tsv')} and host_clade_summary.tsv")

    print("\n=== summary ===")
    size_dist = Counter(len(m) for _, m in clade_meta)
    print("clade size distribution (size:count):")
    for size in sorted(size_dist, reverse=True)[:15]:
        print(f"  n={size}: {size_dist[size]} clades")

    print("\n=== composition ===")
    print(composition_report(genome_accs, ann))

    print("\n=== validation ===")
    ok, msgs = validate(base, outdir, genome_accs, clade_meta, clade_counter, ann)
    for m in msgs:
        print("  " + m)
    if not ok:
        sys.exit(1)
    print("validation: PASS")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("base", help="v2 base dir (e.g. /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2)")
    args = ap.parse_args()
    sys.exit(main(os.path.abspath(args.base)))
