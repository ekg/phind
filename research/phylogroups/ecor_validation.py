#!/usr/bin/env python3
"""
Validate the in-silico ClermonTyping phylogroups against the accepted ECOR
reference phylogroups.

Inputs:
  research/phylogroups/data/phylogroups.tsv   (accession -> phylogroup)
  research/ecor/ecor_phylogroups_known.tsv     (gcf -> known phylogroup)

Prints the agreement table. The "mismatches" are cross-checked against the
discrepancies documented in Waters et al. 2020 (EzClermont paper, Table 2),
where the literature phylogroup differs from the phylogenetic lineage and the
in-silico Clermont tools agree with lineage (or hit documented artifacts).

Usage: python3 ecor_validation.py <phylogroups.tsv>
"""
import csv, os, sys

def main():
    pg_tsv = sys.argv[1]
    ecor_tsv = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ecor", "ecor_phylogroups_known.tsv")

    mine = {}
    for r in csv.DictReader(open(pg_tsv), delimiter="\t"):
        mine[r["accession"]] = r["phylogroup"].strip()

    known = list(csv.DictReader(open(ecor_tsv), delimiter="\t"))
    n = 0; agree = 0
    mismatches = []
    for r in known:
        acc = r["gcf_accession"]; kn = r["known_phylogroup"]
        m = mine.get(acc)
        n += 1
        if m == kn:
            agree += 1
        else:
            mismatches.append((r["ecor_strain"], acc, kn, r["reported_phylogroup"], m))

    print(f"ECOR strains compared: {n}")
    print(f"Agreement: {agree}/{n} = {100.0*agree/n:.1f}%  (requirement: >=90%)")
    print()
    print(f"Mismatches ({len(mismatches)}):")
    for st, acc, kn, rep, m in mismatches:
        print(f"  {st} {acc}: known={kn} (reported={rep}) vs ClermonTyping={m}")

if __name__ == "__main__":
    main()