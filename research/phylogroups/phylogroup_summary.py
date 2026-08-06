#!/usr/bin/env python3
"""
Summarise the ClermonTyping phylogroup calls for the 26k cohort.

Reads research/phylogroups/data/phylogroups.tsv (accession, phylogroup, ...)
Writes:
  research/phylogroups/data/summary_counts.tsv   one row per phylogroup: count
  research/phylogroups/data/assignment_stats.tsv coverage / assigned-rate stats

Usage: python3 phylogroup_summary.py <data/phylogroups.tsv>
"""
import csv, os, sys
from collections import Counter

def main():
    tsv = sys.argv[1]
    outdir = os.path.dirname(tsv)
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))

    n_total = len(rows)
    known_groups = {"A","B1","B2","C","D","E","F","G","U"}
    # assigned = got a concrete single phylogroup (not unknown/non-escherichia/empty)
    skip = {"", "Non Escherichia", "Unknown", "E or cladeI", "NonEscherichia"}
    counter = Counter()
    assigned = 0
    per_note = Counter()
    for r in rows:
        pg = r["phylogroup"].strip()
        counter[pg] += 1
        if pg not in skip:
            assigned += 1
        else:
            per_note[r.get("note","")] += 1

    with open(os.path.join(outdir, "summary_counts.tsv"), "w", newline="") as f:
        f.write("phylogroup\tcount\n")
        for pg, c in sorted(counter.items(), key=lambda x: (-x[1], x[0])):
            f.write(f"{pg}\t{c}\n")

    assigned_pct = 100.0*assigned/n_total if n_total else 0
    with open(os.path.join(outdir, "assignment_stats.tsv"), "w", newline="") as f:
        f.write("metric\tvalue\n")
        f.write(f"total_genomes\t{n_total}\n")
        f.write(f"assigned_phylogroup\t{assigned}\n")
        f.write(f"assigned_pct\t{assigned_pct:.2f}\n")
        f.write(f"unassigned\t{n_total-assigned}\n")

    print(f"total={n_total} assigned={assigned} ({assigned_pct:.2f}%)")
    print("distribution:", dict(counter.most_common()))

if __name__ == "__main__":
    main()