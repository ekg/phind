#!/usr/bin/env python3
"""
Build the authoritative ECOR manifest for the all-prophage collection.

ECOR = the historical E. coli Reference collection (Ochman & Selander 1984),
72 strains (ECOR-1 .. ECOR-72). Their 2018 reference draft genomes
(Patel et al., "Draft Genome Sequences of the Escherichia coli Reference
(ECOR) Collection", Microbiol Resour Announc 7:e01133-18, 2018; BioProject
PRJNA230969) were deposited with GenBank WGS master accessions QOWM00000000.1
(ECOR-1) .. QOZE00000000.1 (ECOR-72); ECOR-59 was deposited out-of-order as
QOZF00000000.1.

Primary strain number source: NCBI assembly_summary `infraspecific_name`
`strain=MOD1-ECOR<n>` labels for BioProject PRJNA230969 (authoritative, exact
for all 72). The paper Table 1 is embedded below as a cross-check that every
strain number is accounted for.

This script crosswalks each ECOR strain to its GenBank assembly (GCA_...) and
its identical RefSeq assembly (GCF_...) from the NCBI assembly summaries, then
joins the GCF accessions to the 26k prophage collection keyed by the `genome`
column of 26k_prophage1.csv and to the extracted full_prophages.fa headers.

Coordinate policy: 26k_prophage1.csv begin/end are 1-based inclusive
(prophage-semantics v2, C1_RAW_1_BASED_CLOSED); length = end - begin + 1.

Outputs (under research/ecor/):
  ecor_manifest.csv               - one row per ECOR prophage element
  ecor_strain_reconciliation.tsv  - per-strain status/counts
  ecor_reconciliation_summary.json - exact counts + tag merge rate
"""
from __future__ import annotations

import csv
import json
import re
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ECOR_DIR = ROOT / "research" / "ecor"

GENBANK_SUMMARY = ECOR_DIR / "ecori_genbank.tsv"       # taxid-562 filtered
REFSEC_SUMMARY = ECOR_DIR / "ecori_refseq.tsv"
HEADER = ECOR_DIR / "header.txt"
COHORT_ACC = ROOT / "26k_ecoli_accession.txt"
PROPHAGE_CSV = ROOT / "26k_prophage1.csv"
FULL_PROPHAGES_FA = ROOT / "prophage_homology_survey" / "full_prophages.fa"

OUT_MANIFEST = ECOR_DIR / "ecor_manifest.csv"
OUT_RECONC = ECOR_DIR / "ecor_strain_reconciliation.tsv"
OUT_SUMMARY = ECOR_DIR / "ecor_reconciliation_summary.json"
OUT_ALTERNATES = ECOR_DIR / "ecor_alternate_assemblies.tsv"
OUT_LEAF_TAGS = ECOR_DIR / "ecor_leaf_tags.tsv"

# --- Paper Table 1 (Patel et al. 2018, MRA 7:e01133-18) ---
# ECOR strain -> GenBank WGS master accession (for cross-check only).
ECOR_WGS_PAPER: dict[int, str] = {
    1: "QOWM00000000.1", 2: "QOWN00000000.1", 3: "QOWO00000000.1",
    4: "QOWP00000000.1", 5: "QOWQ00000000.1", 6: "QOWR00000000.1",
    7: "QOWS00000000.1", 8: "QOWT00000000.1", 9: "QOWU00000000.1",
    10: "QOWV00000000.1", 11: "QOWW00000000.1", 12: "QOWX00000000.1",
    13: "QOWY00000000.1", 14: "QOWZ00000000.1", 15: "QOXA00000000.1",
    16: "QOXB00000000.1", 17: "QOXC00000000.1", 18: "QOXD00000000.1",
    19: "QOXE00000000.1", 20: "QOXF00000000.1", 21: "QOXG00000000.1",
    22: "QOXH00000000.1", 23: "QOXI00000000.1", 24: "QOXJ00000000.1",
    25: "QOXK00000000.1", 26: "QOXL00000000.1", 27: "QOXM00000000.1",
    28: "QOXN00000000.1", 29: "QOXO00000000.1", 30: "QOXP00000000.1",
    31: "QOXQ00000000.1", 32: "QOXR00000000.1", 33: "QOXS00000000.1",
    34: "QOXT00000000.1", 35: "QOXU00000000.1", 36: "QOXV00000000.1",
    37: "QOXW00000000.1", 38: "QOXX00000000.1", 39: "QOXY00000000.1",
    40: "QOXZ00000000.1", 41: "QOYA00000000.1", 42: "QOYB00000000.1",
    43: "QOYC00000000.1", 44: "QOYD00000000.1", 45: "QOYE00000000.1",
    46: "QOYF00000000.1", 47: "QOYG00000000.1", 48: "QOYH00000000.1",
    49: "QOYI00000000.1", 50: "QOYJ00000000.1", 51: "QOYK00000000.1",
    52: "QOYL00000000.1", 53: "QOYM00000000.1", 54: "QOYN00000000.1",
    55: "QOYO00000000.1", 56: "QOYP00000000.1", 57: "QOYQ00000000.1",
    58: "QOYR00000000.1", 59: "QOZF00000000.1", 60: "QOYS00000000.1",
    61: "QOYT00000000.1", 62: "QOYU00000000.1", 63: "QOYV00000000.1",
    64: "QOYW00000000.1", 65: "QOYX00000000.1", 66: "QOYY00000000.1",
    67: "QOYZ00000000.1", 68: "QOZA00000000.1", 69: "QOZB00000000.1",
    70: "QOZC00000000.1", 71: "QOZD00000000.1", 72: "QOZE00000000.1",
}


# Paper footnote (a): possible contamination in the Messerer et al. 2017
# (PRJNA224116) assembly of that strain (multiple O/H serotyping loci).
# Marks apply to the Messerer column, NOT to the 2018 reference assemblies.
ECOR_PAPER_CONTAM_FLAG: set[int] = {1, 2, 4, 5, 6, 8, 10, 12, 14, 30, 33, 34, 38, 39, 40,
                                     42, 43, 44, 45, 46, 47, 48, 51, 57, 58, 59, 60, 61, 62, 66, 70, 71}


def load_summary(path: Path) -> list[dict[str, str]]:
    with open(HEADER) as f:
        header = f.readline().strip().lstrip("#").split("\t")
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            rows.append(dict(zip(header, parts)))
    return rows


def main() -> int:
    genbank = load_summary(GENBANK_SUMMARY)
    refseq = load_summary(REFSEC_SUMMARY)
    refseq_accs = {r.get("assembly_accession") for r in refseq}

    # 1) Derive ECOR strains from NCBI strain labels in PRJNA230969
    strain_map: OrderedDict[int, dict] = OrderedDict()
    by_wgs: dict[str, dict] = {}
    for r in genbank:
        by_wgs.setdefault(r.get("wgs_master", ""), r)
        m = re.search(r"MOD1[-_ ]?ECOR[-_ ]?(\d{1,2})\b",
                      r.get("infraspecific_name", ""), re.IGNORECASE)
        if not m:
            continue
        n = int(m.group(1))
        if n in strain_map:
            # duplicate strain number across rows: keep latest status, flag
            strain_map[n]["duplicate_rows"] = strain_map[n].get("duplicate_rows", 0) + 1
            continue
        strain_map[n] = {
            "strain": n,
            "wgs": r.get("wgs_master", ""),
            "gca": r.get("assembly_accession", ""),
            "gcf": r.get("gbrs_paired_asm", "") or "",
            "paired_comp": r.get("paired_asm_comp", ""),
            "excluded": r.get("excluded_from_refseq", ""),
            "bioproject": r.get("bioproject", ""),
            "version_status": r.get("version_status", ""),
            "assembly_level": r.get("assembly_level", ""),
            "duplicate_rows": 0,
        }

    missing_nums = [n for n in range(1, 73) if n not in strain_map]
    problems: list[str] = []
    if missing_nums:
        problems.append(f"no NCBI label row for ECOR strains: {missing_nums}")
    if len(strain_map) != 72:
        problems.append(f"expected 72 ECOR strains from labels, found {len(strain_map)}")

    # 2) Cross-check wgs against paper Table 1
    paper_by_wgs = {v: k for k, v in ECOR_WGS_PAPER.items()}
    paper_discrepancies = []
    for n, info in strain_map.items():
        if info["wgs"] not in paper_by_wgs:
            paper_discrepancies.append(f"ECOR-{n}: NCBI wgs {info['wgs']} not in paper table")
        elif paper_by_wgs[info["wgs"]] != n:
            paper_discrepancies.append(
                f"ECOR-{n}: NCBI wgs {info['wgs']} is paper ECOR-{paper_by_wgs[info['wgs']]}")
    # reverse: paper strains missing from NCBI labels
    for n, wgs in ECOR_WGS_PAPER.items():
        if n not in strain_map:
            paper_discrepancies.append(f"paper ECOR-{n} ({wgs}) has no NCBI label row")
    if paper_discrepancies:
        problems.append("paper-table discrepancies: " + "; ".join(paper_discrepancies))

    # 3) Validate GCF presence in refseq summary + cohort
    cohort = set(l.strip() for l in COHORT_ACC.open() if l.strip())
    for n, info in strain_map.items():
        gcf = info["gcf"]
        if gcf:
            info["gcf_in_refseq_summary"] = gcf in refseq_accs
            if not info["gcf_in_refseq_summary"]:
                problems.append(f"ECOR-{n}: GCF {gcf} not in refseq summary")
        else:
            info["gcf_in_refseq_summary"] = False
        info["in_cohort"] = gcf in cohort if gcf else False

    # 4) Load prophage records keyed by genome accession
    proph: dict[str, list[dict]] = defaultdict(list)
    with PROPHAGE_CSV.open(newline="") as f:
        for row in csv.DictReader(f):
            proph[row["genome"]].append(row)

    # 5) Load full_prophages.fa headers
    fa_headers: set[str] = set()
    with FULL_PROPHAGES_FA.open("rb") as f:
        for line in f:
            if line.startswith(b">"):
                fa_headers.add(line[1:].strip().decode())
    print(f"full_prophages.fa headers: {len(fa_headers)}", file=sys.stderr)

    # 6) Build manifest rows + per-strain reconciliation
    manifest_rows: list[dict] = []
    per_strain: OrderedDict[int, dict] = OrderedDict()
    for n in range(1, 73):
        info = strain_map.get(n)
        if info is None:
            per_strain[n] = {
                "ecor_strain": f"ECOR-{n}",
                "wgs_master": ECOR_WGS_PAPER.get(n, ""),
                "gca": "", "gcf": "",
                "in_cohort": False,
                "n_prophage_records": 0, "n_in_full_prophages_fa": 0,
                "tag_merge_rate": 0.0,
                "status": "not-found-by-accession",
                "notes": "no NCBI strain-label row for ECOR number",
            }
            continue
        gcf = info["gcf"] or ""
        records = proph.get(gcf, []) if gcf else []
        n_el = len(records)
        n_fa = sum(1 for r in records if r["prophage_id"] in fa_headers)
        if info["in_cohort"]:
            status = "matched"
        elif gcf and info["gcf_in_refseq_summary"]:
            status = "gcf-resolved-not-in-cohort"
        elif gcf:
            status = "gcf-not-in-refseq-summary"
        else:
            status = "not-found-by-accession"
        per_strain[n] = {
            "ecor_strain": f"ECOR-{n}",
            "wgs_master": info["wgs"],
            "gca": info["gca"],
            "gcf": gcf,
            "in_cohort": info["in_cohort"],
            "n_prophage_records": n_el,
            "n_in_full_prophages_fa": n_fa,
            "tag_merge_rate": round(n_fa / n_el, 4) if n_el else 0.0,
            "status": status,
            "notes": f"paired_comp={info['paired_comp']}; level={info['assembly_level']}",
        }
        for r in sorted(records, key=lambda x: (x["scaffold"], int(float(x["begin"])))):
            b = int(float(r["begin"]))
            e = int(float(r["end"]))
            manifest_rows.append({
                "ecor_strain": f"ECOR-{n}",
                "assembly_accession": gcf,
                "gca_accession": info["gca"],
                "wgs_master": info["wgs"],
                "prophage_id": r["prophage_id"],
                "source_contig": r["scaffold"],
                "start": b,
                "end": e,
                "length": e - b + 1,
                "transposable": r["transposable"],
                "taxonomy": r["taxonomy"],
            })

    # 7) Ambiguity analysis: other ECOR-labeled GCF assemblies in the cohort
    #    (e.g. BioProject PRJNA224116, Messerer et al. 2017, flagged as
    #    potentially contaminated in the 2018 reference paper).
    alternate_rows: list[dict] = []
    strain_has_alternate: dict[int, int] = defaultdict(int)
    for r in refseq:
        if r.get("assembly_accession") not in cohort:
            continue
        blob = "|".join([r.get("infraspecific_name", ""), r.get("isolate", ""),
                         r.get("asm_name", ""), r.get("organism_name", "")])
        m = re.search(r"ECOR[-_ ]?(\d{1,2})\b", blob, re.IGNORECASE)
        if not m:
            continue
        n = int(m.group(1))
        gb = by_wgs.get(r.get("wgs_master", ""), {})
        if r["assembly_accession"] != per_strain[n]["gcf"]:
            strain_has_alternate[n] += 1
            alternate_rows.append({
                "ecor_strain": f"ECOR-{n}",
                "assembly_accession": r["assembly_accession"],
                "bioproject": r.get("bioproject", ""),
                "wgs_master": r.get("wgs_master", ""),
                "asm_submitter": r.get("asm_submitter", ""),
                "seq_rel_date": r.get("seq_rel_date", ""),
                "assembly_level": r.get("assembly_level", ""),
                "strain_label": blob.replace("|", "; "),
                "paper_contamination_flag": "yes" if n in ECOR_PAPER_CONTAM_FLAG else "no",
                "note": "alternate ECOR-labeled assembly in cohort; canonical is the 2018 PRJNA230969 reference",
            })
    ambiguous_strains = sorted(strain_has_alternate)

    # 8) Reconciliation summary
    status_counts: dict[str, int] = defaultdict(int)
    for p in per_strain.values():
        status_counts[p["status"]] += 1
    gcf_counts = {"strains_with_gcf": sum(1 for p in per_strain.values() if p["gcf"]),
                  "strains_gcf_in_cohort": sum(1 for p in per_strain.values() if p["in_cohort"]),
                  "strains_gcf_not_in_cohort": sum(1 for p in per_strain.values() if p["gcf"] and not p["in_cohort"])}
    n_el = sum(p["n_prophage_records"] for p in per_strain.values())
    n_fa = sum(p["n_in_full_prophages_fa"] for p in per_strain.values())

    summary = {
        "schema": "ecor-manifest-reconciliation-v1",
        "canonical_ecor_strains": 72,
        "reference": {
            "paper": "Patel et al. 2018, Microbiol Resour Announc 7:e01133-18 (Draft Genome Sequences of the Escherichia coli Reference (ECOR) Collection)",
            "bioproject": "PRJNA230969",
            "wgs_range": "QOWM00000000.1..QOZE00000000.1 (ECOR-59 = QOZF00000000.1)",
            "strain_number_source": "NCBI assembly_summary infraspecific_name 'strain=MOD1-ECOR<n>'",
        },
        "strain_status_counts": dict(status_counts),
        "coverage": {**gcf_counts,
                     "fraction_gcf_in_cohort": round(gcf_counts["strains_gcf_in_cohort"] / 72, 4)},
        "ambiguity": {
            "strains_with_alternate_ecor_labeled_assemblies_in_cohort": len(ambiguous_strains),
            "alternate_assembly_count": len(alternate_rows),
            "note": "alternates are BioProject PRJNA224116 (Messerer et al. 2017) and a few others; the 2018 reference paper flags ~half of PRJNA224116 assemblies as potentially contaminated; canonical mapping uses PRJNA230969.",
        },
        "elements": {
            "total_ecor_prophage_records": n_el,
            "records_in_full_prophages_fa": n_fa,
            "tag_merge_rate": round(n_fa / n_el, 6) if n_el else 0.0,
            "total_manifest_rows": len(manifest_rows),
            "sum_element_bp": sum(r["length"] for r in manifest_rows),
        },
        "problems": problems,
    }

    # 9) Leaf tags for the MASH tree (one row per full_prophages.fa header)
    id_to_strain = {r["prophage_id"]: r["ecor_strain"] for r in manifest_rows}
    with OUT_LEAF_TAGS.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["prophage_id", "is_ecor", "ecor_strain"])
        for hid in sorted(fa_headers):
            strain = id_to_strain.get(hid, "")
            w.writerow([hid, "TRUE" if strain else "FALSE", strain])

    # 8) Write outputs
    ECOR_DIR.mkdir(parents=True, exist_ok=True)
    with OUT_MANIFEST.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "ecor_strain", "assembly_accession", "gca_accession", "wgs_master",
            "prophage_id", "source_contig", "start", "end", "length",
            "transposable", "taxonomy"])
        w.writeheader()
        w.writerows(manifest_rows)
    with OUT_RECONC.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "ecor_strain", "wgs_master", "gca", "gcf", "in_cohort",
            "n_prophage_records", "n_in_full_prophages_fa", "tag_merge_rate",
            "status", "notes"], delimiter="\t")
        w.writeheader()
        for p in per_strain.values():
            w.writerow({**p, "in_cohort": str(p["in_cohort"]),
                        "n_prophage_records": str(p["n_prophage_records"]),
                        "n_in_full_prophages_fa": str(p["n_in_full_prophages_fa"])})
    with OUT_ALTERNATES.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "ecor_strain", "assembly_accession", "bioproject", "wgs_master",
            "asm_submitter", "seq_rel_date", "assembly_level", "strain_label",
            "paper_contamination_flag", "note"], delimiter="\t")
        w.writeheader()
        for a in sorted(alternate_rows, key=lambda x: (int(x["ecor_strain"].split("-")[1]), x["assembly_accession"])):
            w.writerow(a)
    with OUT_SUMMARY.open("w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
