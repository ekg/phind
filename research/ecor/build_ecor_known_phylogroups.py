#!/usr/bin/env python3
"""
Build research/ecor/ecor_phylogroups_known.tsv -- the accepted (ground-truth)
Clermont phylogroup per ECOR strain, used to validate the in-silico ClermonTyping
assignment (task in-silico-e).

Ground-truth source: the ECOR phylogroups recorded in the EzClermont validation
metadata (Waters et al. 2020, https://github.com/nickp60/EzClermont,
docs/analysis/validate/validation_metadata.csv), whose "reported_phylogroup"
column is taken from the Clermont 2013/2015 phylotyping literature. All 72 ECOR
strains are present.

The recorded values carry Clermont subtype suffixes (e.g. "B2 II", "D (CGA +)").
For validation we compare only the broad phylogroup (A, B1, B2, C, D, E, F), the
category granularity produced by ClermonTyping for E. coli.

Output columns: ecor_strain, gcf_accession (RefSeq paired accession used in the
26k cohort), known_phylogroup, reported_phylogroup.
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
STRAIN_RECON = os.path.join(HERE, "ecor_strain_reconciliation.tsv")
KNOWN_JSON = "/tmp/ecor_pg.json"          # strain -> {reported, broad}
OUT = os.path.join(HERE, "ecor_phylogroups_known.tsv")

def broad(pg: str) -> str:
    # "B2 II" -> B2 ; "D (CGA +)" -> D ; "B2 UA" -> B2 ; "G" -> G
    return pg.split(" ")[0].split("(")[0].strip()

def main():
    # strain -> gcf from the reconciliation table
    strain_to_gcf = {}
    with open(STRAIN_RECON) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            strain_to_gcf[r["ecor_strain"]] = r["gcf"]

    known = json.load(open(KNOWN_JSON))

    missing_gcf = [s for s in known if s not in strain_to_gcf]
    if missing_gcf:
        sys.exit(f"ERROR: no GCF for strains: {missing_gcf}")

    rows = []
    for strain in sorted(known, key=lambda s: int(s.split("-")[1])):
        rows.append({
            "ecor_strain": strain,
            "gcf_accession": strain_to_gcf[strain],
            "known_phylogroup": known[strain]["broad"],
            "reported_phylogroup": known[strain]["reported"],
        })

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # report distribution
    from collections import Counter
    dist = Counter(r["known_phylogroup"] for r in rows)
    print(f"Wrote {OUT} with {len(rows)} strains")
    print("Broad phylogroup distribution:", dict(sorted(dist.items())))

if __name__ == "__main__":
    main()