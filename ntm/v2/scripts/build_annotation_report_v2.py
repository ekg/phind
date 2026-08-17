#!/usr/bin/env python3
"""
build_annotation_report_v2.py — per-genome functional-QC report for the NTM v2
generated phage genomes (2,388 ML + 1,251 ancestral = 3,639), from Pharokka +
CheckV outputs produced by run_annotation_v2.py.

Adapted from the v1 builder (ntm/scripts/build_annotation_report.py) with two
changes requested for v2:

  1. ML genomes are the PRIMARY cohort and ancestral genomes an explicit
     COMPARATOR. Summaries are broken out by cohort:
       ntm2_ml_reconstructed  (release status=ml,   >=2 clade members)
       ntm2_ml_singleton      (release status=singleton, pass-through fragment)
       ntm2_anc               (ancestral-state genomes)
     and by source mode (ntm2_ml / ntm2_anc).

  2. A graded `candidate_functionality` classification replaces any binary
     functional/non-functional claim, and never equates annotation failure
     with biological absence (see CLASSIFICATION below).

Inputs (all under the v2 analysis root, default
/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/annotation):
  input/genome_index.tsv                                 (this pipeline)
  pharokka_out/pharokka_cds_final_merged_output.tsv      (per-CDS)
  pharokka_out/pharokka_length_gc_cds_density.tsv        (per-contig)
  checkv_out/quality_summary.tsv                         (per-contig)

Outputs (ROOT/report/):
  per_genome_functional_qc.tsv          one row per genome, sorted by genome_id
  summary_by_cohort.tsv                 cohort-level aggregates
  summary_by_source.tsv                 source-mode aggregates
  candidate_functionality_by_cohort.tsv tier x cohort matrix

Determinism: rows sorted by genome_id, fixed 3-decimal rounding, LF endings,
no timestamps — reruns are byte-identical (SHA-256 stable).

CLASSIFICATION (graded candidate tiers; precedence top-down)
------------------------------------------------------------
A_strong_candidate        CheckV Complete/High/Medium AND coherent head-packaging
                          (>=2 of terminase/portal/capsid) AND tail AND lysis
                          AND truncated_frac <= 0.2 AND contamination == 0.
B_moderate_candidate      CheckV Complete/High/Medium AND >=1 head-packaging
                          module AND (tail OR lysis) AND truncated_frac <= 0.35
                          AND contamination == 0 — core present but not the full
                          coherent set of tier A.
C_partial_module_evidence >=1 key module detected, but either the module set or
                          the CheckV grade is insufficient for A/B.
D_modules_not_detected    No key-module annotation detected. This is a DETECTION
                          result under mmseqs2-only PHROG matching, NOT proof of
                          biological absence: novel/divergent ORFans (common in
                          mycobacteriophages), short fragments, and annotation
                          depth all reduce detection.
no_annotation             Pharokka called 0 CDS for the record — an identifier /
                          tool failure to investigate, never a biological claim.

Integrase is recorded (temperate prophage support) but deliberately NOT part of
the tier ladder: it is expected for prophage-derived genomes, not universally
required of a functional phage genome.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import statistics
import sys
from collections import defaultdict

DEFAULT_ROOT = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/annotation"

# --- key-module detection (same matching as v1, kept for comparability) ------
KEY = {
    "terminase": re.compile(r"terminase", re.I),
    "portal": re.compile(r"portal", re.I),
    "capsid": re.compile(r"capsid|major head|head protein", re.I),
    "tail": re.compile(r"tail|sheath|base ?plate|tail ?fiber|tail ?tube", re.I),
    "integrase": re.compile(r"integrase", re.I),
    # v1 regex was r"holin|lysin|endolysin|spanin|R|lysis" — the bare "R"
    # alternative matched ANY capital R (e.g. "RNA polymerase"), over-calling
    # lysis. "R"/"Rz" are spanin subunit gene names: require word boundaries.
    "lysis": re.compile(r"holin|lysin|endolysin|spanin|\bRz?\b|lysis", re.I),
}
# PHROG category -> inferred module (category evidence supplements annot text)
CAT2KEY = {
    "tail": "tail",
    "integration and excision": "integrase",
    "lysis": "lysis",
}
HEAD_MODULES = ("terminase", "portal", "capsid")
KEYLIST = ["terminase", "portal", "capsid", "tail", "integrase", "lysis"]

CHECKV_GRADE_PASS = {"complete", "high-quality", "medium-quality"}

# cohort derivation from genome_id: ntm2_<clade>_ML / ntm2_<clade>_ANCESTRAL
RE_ML = re.compile(r"^ntm2_.+_ML$")
RE_ANC = re.compile(r"^ntm2_.+_ANCESTRAL$")


# ----------------------------------------------------------------------------
# input parsing
# ----------------------------------------------------------------------------

def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_genome_index(path):
    """genome_index.tsv -> {genome_id: row}. Validates ID uniqueness."""
    rows = read_tsv(path)
    meta = {}
    for r in rows:
        gid = r["genome_id"]
        if gid in meta:
            raise ValueError(f"duplicate genome_id in index: {gid}")
        meta[gid] = r
    return meta


def cohort_of(genome_id: str, index_row: dict) -> str:
    """Primary/comparator cohort label.

    status=ml -> ntm2_ml_reconstructed; status=singleton -> ntm2_ml_singleton;
    ancestral FASTA -> ntm2_anc. Falls back to parsing the genome_id suffix.
    """
    status = (index_row.get("status") or "").strip()
    if status == "ml":
        return "ntm2_ml_reconstructed"
    if status == "singleton":
        return "ntm2_ml_singleton"
    if status == "ancestral":
        return "ntm2_anc"
    if RE_ML.match(genome_id):
        return "ntm2_ml_reconstructed"
    if RE_ANC.match(genome_id):
        return "ntm2_anc"
    raise ValueError(f"cannot derive cohort for genome_id={genome_id!r} "
                     f"(status={status!r})")


def parse_pharokka_cds(path, keep_ids=None):
    """-> per-genome {gene_count, truncated, modules:set, cats:set, hypothetical}"""
    genes = defaultdict(int)
    partial = defaultdict(int)
    keys = defaultdict(set)
    cats = defaultdict(set)
    hypoth = defaultdict(int)
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            gid = r["contig"]
            if keep_ids is not None and gid not in keep_ids:
                continue
            genes[gid] += 1
            if (r.get("partial") or "00").strip() != "00":
                partial[gid] += 1
            annot = r.get("annot", "") or ""
            for k, rx in KEY.items():
                if annot and rx.search(annot):
                    keys[gid].add(k)
            cat = (r.get("category", "") or "").strip().lower()
            if cat:
                cats[gid].add(cat)
                if cat in CAT2KEY:
                    keys[gid].add(CAT2KEY[cat])
            if cat == "unknown function" or (annot and "hypothetical" in annot.lower()):
                hypoth[gid] += 1
    return genes, partial, keys, cats, hypoth


def normalize_checkv_grade(raw: str) -> str:
    g = (raw or "").strip()
    low = g.lower()
    if low in ("complete", "high-quality", "medium-quality",
               "low-quality", "not-determined", "genome fragment",
               "provirus", "proviral contig"):
        return g if g else "Not-determined"
    if "not determined" in low:
        return "Not-determined"
    if "low quality" in low:
        return "Low-quality"
    return g or "Not-determined"


def parse_checkv(path):
    cv = {}
    for r in read_tsv(path):
        cv[r["contig_id"]] = r
    return cv


# ----------------------------------------------------------------------------
# classification
# ----------------------------------------------------------------------------

def _to_float(x, default=None):
    try:
        return float(str(x).strip())
    except (TypeError, ValueError):
        return default


def classify(gene_count, modules, checkv_grade, truncated_frac,
             contamination, contamination_known=True):
    """Graded candidate-functionality tier.

    Returns (tier, rationale). Missing-module outcomes are 'not detected'
    statements about the annotation evidence, never claims of biological
    absence.
    """
    head_n = len([m for m in HEAD_MODULES if m in modules])
    head_ok = head_n >= 2
    tail_ok = "tail" in modules
    lysis_ok = "lysis" in modules
    grade_pass = checkv_grade.lower() in CHECKV_GRADE_PASS
    contam = contamination if contamination_known else None

    bits = []
    bits.append(f"checkv={checkv_grade or 'missing'}")
    bits.append(f"head={head_n}/3[{'+'.join(m for m in HEAD_MODULES if m in modules) or '-'}]")
    bits.append(f"tail={int(tail_ok)}")
    bits.append(f"lysis={int(lysis_ok)}")
    bits.append(f"integrase={int('integrase' in modules)}(support-only)")
    bits.append(f"trunc_frac={truncated_frac:.3f}" if truncated_frac is not None else "trunc_frac=NA")
    bits.append(f"contam={contamination:g}" if contam is not None else "contam=NA")

    if gene_count == 0:
        tier = "no_annotation"
        bits.append("0 CDS from Pharokka — identifier/tool failure to investigate, not a biological claim")
    elif contam is not None and contam > 5.0:
        tier = "C_partial_module_evidence"
        bits.append("contamination>5% makes module evidence unreliable")
    elif (grade_pass and head_ok and tail_ok and lysis_ok
          and truncated_frac is not None and truncated_frac <= 0.20
          and contam is not None and contam == 0.0):
        tier = "A_strong_candidate"
        bits.append("all tier-A criteria met")
    elif (grade_pass and head_n >= 1 and (tail_ok or lysis_ok)
          and truncated_frac is not None and truncated_frac <= 0.35
          and (contam is None or contam == 0.0)):
        tier = "B_moderate_candidate"
        missing = [n for n, ok in (("head>=2", head_ok), ("tail", tail_ok), ("lysis", lysis_ok)) if not ok]
        bits.append("core modules present, incomplete vs tier A: " + (",".join(missing) or "hygiene"))
    elif len(modules & set(KEYLIST)) >= 1:
        tier = "C_partial_module_evidence"
        why = []
        if not grade_pass:
            why.append(f"checkv grade '{checkv_grade}' below Medium")
        if head_n == 0:
            why.append("no head-packaging module detected")
        elif not head_ok:
            why.append(f"only {head_n}/3 head-packaging modules detected")
        if not tail_ok and head_ok:
            why.append("tail not detected")
        if not lysis_ok and head_ok:
            why.append("lysis not detected")
        if contam is not None and contam > 0:
            why.append(f"contamination {contam:g}%")
        if truncated_frac is not None and truncated_frac > 0.35:
            why.append(f"truncated_frac {truncated_frac:.2f} > 0.35")
        bits.append("partial module evidence: " + "; ".join(why))
    else:
        tier = "D_modules_not_detected"
        bits.append("no key-module PHROG/mmseqs2 hits — detection gap (novel/divergent ORFan, "
                    "fragment, or annotation depth), NOT proof of biological absence")
    return tier, "; ".join(bits)


# ----------------------------------------------------------------------------
# report building
# ----------------------------------------------------------------------------

PER_GENOME_COLS = [
    "genome_id", "source", "cohort", "clade_id", "length", "gc_perc",
    "cds_density", "gene_count", "truncated_genes", "truncated_frac",
    "terminase", "portal", "capsid", "tail", "integrase", "lysis",
    "n_key_proteins", "head_packaging_modules", "has_head_packaging_cat",
    "has_tail_cat", "has_lysis_cat", "has_integration_cat", "temperate_signal",
    "hypothetical", "hypothetical_frac", "checkv_miuvig", "checkv_grade",
    "completeness_pct", "completeness_method", "contamination", "viral_genes",
    "host_genes", "candidate_functionality", "tier_rationale", "flag",
]


def build_rows(index_path, cds_path, lgd_path, checkv_path,
               missing_checkv_ok=False):
    meta = load_genome_index(index_path)
    genes, partial, keys, cats, hypoth = parse_pharokka_cds(cds_path, keep_ids=set(meta))
    lgd = {r["contig"]: r for r in read_tsv(lgd_path)}
    cv = parse_checkv(checkv_path) if checkv_path and os.path.exists(checkv_path) else {}

    missing_checkv = sorted(set(meta) - set(cv))
    if missing_checkv and not missing_checkv_ok:
        raise RuntimeError(
            f"{len(missing_checkv)} genomes missing from CheckV quality_summary "
            f"(first: {missing_checkv[:3]}) — run CheckV before building the report "
            f"or pass --missing-checkv-ok")

    rows = []
    for gid in sorted(meta):
        m = meta[gid]
        g = genes.get(gid, 0)
        ks = keys.get(gid, set())
        lg = lgd.get(gid, {})
        c = cv.get(gid, {})
        grade = normalize_checkv_grade(c.get("checkv_quality"))
        contam = _to_float(c.get("contamination"))
        trunc = partial.get(gid, 0)
        trunc_frac = round(trunc / g, 3) if g else None
        tier, rationale = classify(g, ks, grade, trunc_frac, contam)

        flags = []
        if g == 0:
            flags.append("no_genes")
        if g > 0 and not ({"terminase", "portal", "capsid"} & ks):
            flags.append("missing_core_structural")
        if grade.lower() in ("low-quality", "not-determined", "genome fragment"):
            flags.append("low_completeness")
        if contam is not None and contam > 0:
            flags.append("contamination")
        hgenes = c.get("host_genes") or ""
        try:
            if hgenes and int(hgenes) > 0:
                flags.append("host_genes")
        except ValueError:
            pass
        if trunc_frac is not None and trunc_frac > 0.2:
            flags.append("many_truncated")

        rows.append({
            "genome_id": gid,
            "source": m["source"],
            "cohort": cohort_of(gid, m),
            "clade_id": m.get("clade_id", ""),
            "length": m.get("length", ""),
            "gc_perc": lg.get("gc_perc", ""),
            "cds_density": lg.get("cds_coding_density", ""),
            "gene_count": g,
            "truncated_genes": trunc,
            "truncated_frac": trunc_frac if trunc_frac is not None else "",
            "terminase": int("terminase" in ks), "portal": int("portal" in ks),
            "capsid": int("capsid" in ks), "tail": int("tail" in ks),
            "integrase": int("integrase" in ks), "lysis": int("lysis" in ks),
            "n_key_proteins": len(ks & set(KEYLIST)),
            "head_packaging_modules": len([k for k in HEAD_MODULES if k in ks]),
            "has_head_packaging_cat": int("head and packaging" in cats.get(gid, set())),
            "has_tail_cat": int("tail" in cats.get(gid, set())),
            "has_lysis_cat": int("lysis" in cats.get(gid, set())),
            "has_integration_cat": int("integration and excision" in cats.get(gid, set())),
            "temperate_signal": int("integrase" in ks),
            "hypothetical": hypoth.get(gid, 0),
            "hypothetical_frac": round(hypoth.get(gid, 0) / g, 3) if g else "",
            "checkv_miuvig": (c.get("miuvig_quality") or "").strip(),
            "checkv_grade": grade,
            "completeness_pct": (c.get("completeness") or "").strip(),
            "completeness_method": (c.get("completeness_method") or "").strip(),
            "contamination": (c.get("contamination") or "").strip(),
            "viral_genes": (c.get("viral_genes") or "").strip(),
            "host_genes": hgenes.strip(),
            "candidate_functionality": tier,
            "tier_rationale": rationale,
            "flag": ";".join(flags) or "ok",
        })
    return rows, missing_checkv


SUMMARY_COLS = [
    "group", "n_genomes", "median_genes", "median_key_proteins",
    "pct_with_terminase", "pct_with_capsid", "pct_with_tail", "pct_with_lysis",
    "pct_temperate", "pct_missing_core", "pct_low_completeness",
    "pct_contaminated", "median_truncated_frac", "median_hypothetical_frac",
    "pct_tier_A", "pct_tier_B", "pct_tier_C", "pct_tier_D",
    "pct_no_annotation",
]

TIERS = ["A_strong_candidate", "B_moderate_candidate",
         "C_partial_module_evidence", "D_modules_not_detected", "no_annotation"]


def summarize(rows, key):
    groups = defaultdict(list)
    for r in rows:
        groups[r[key]].append(r)
    out = []
    for grp, rs in sorted(groups.items()):
        n = len(rs)

        def med(k, cast=float):
            vals = [cast(r[k]) for r in rs if r[k] != "" and r[k] is not None]
            return round(statistics.median(vals), 3) if vals else ""

        def pct(f):
            return round(100 * sum(1 for r in rs if f(r)) / n, 1)

        out.append({
            "group": grp, "n_genomes": n,
            "median_genes": int(med("gene_count")) if med("gene_count") != "" else "",
            "median_key_proteins": med("n_key_proteins", int),
            "pct_with_terminase": pct(lambda r: r["terminase"]),
            "pct_with_capsid": pct(lambda r: r["capsid"]),
            "pct_with_tail": pct(lambda r: r["tail"]),
            "pct_with_lysis": pct(lambda r: r["lysis"]),
            "pct_temperate": pct(lambda r: r["temperate_signal"]),
            "pct_missing_core": pct(lambda r: "missing_core_structural" in r["flag"]),
            "pct_low_completeness": pct(lambda r: "low_completeness" in r["flag"]),
            "pct_contaminated": pct(lambda r: "contamination" in r["flag"]),
            "median_truncated_frac": med("truncated_frac"),
            "median_hypothetical_frac": med("hypothetical_frac"),
            "pct_tier_A": pct(lambda r: r["candidate_functionality"] == TIERS[0]),
            "pct_tier_B": pct(lambda r: r["candidate_functionality"] == TIERS[1]),
            "pct_tier_C": pct(lambda r: r["candidate_functionality"] == TIERS[2]),
            "pct_tier_D": pct(lambda r: r["candidate_functionality"] == TIERS[3]),
            "pct_no_annotation": pct(lambda r: r["candidate_functionality"] == TIERS[4]),
        })
    return out


def tier_matrix(rows):
    by = defaultdict(lambda: defaultdict(int))
    for r in rows:
        by[r["cohort"]][r["candidate_functionality"]] += 1
    cohorts = sorted(by)
    out = [["cohort"] + TIERS + ["total"]]
    for c in cohorts:
        out.append([c] + [by[c].get(t, 0) for t in TIERS] + [sum(by[c].values())])
    return out


def write_tsv(path, cols, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--root", default=DEFAULT_ROOT,
                    help=f"v2 analysis root (default {DEFAULT_ROOT})")
    ap.add_argument("--outdir", default=None,
                    help="report output dir (default ROOT/report)")
    ap.add_argument("--missing-checkv-ok", action="store_true",
                    help="do not fail when genomes are missing from CheckV "
                         "output (they are reported with empty CheckV fields)")
    args = ap.parse_args(argv)
    root = args.root
    outdir = args.outdir or os.path.join(root, "report")

    rows, missing_checkv = build_rows(
        os.path.join(root, "input", "genome_index.tsv"),
        os.path.join(root, "pharokka_out", "pharokka_cds_final_merged_output.tsv"),
        os.path.join(root, "pharokka_out", "pharokka_length_gc_cds_density.tsv"),
        os.path.join(root, "checkv_out", "quality_summary.tsv"),
        missing_checkv_ok=args.missing_checkv_ok,
    )

    write_tsv(os.path.join(outdir, "per_genome_functional_qc.tsv"),
              PER_GENOME_COLS, rows)
    write_tsv(os.path.join(outdir, "summary_by_cohort.tsv"),
              SUMMARY_COLS, summarize(rows, "cohort"))
    write_tsv(os.path.join(outdir, "summary_by_source.tsv"),
              SUMMARY_COLS, summarize(rows, "source"))
    mat = tier_matrix(rows)
    with open(os.path.join(outdir, "candidate_functionality_by_cohort.tsv"),
              "w", newline="") as f:
        csv.writer(f, delimiter="\t", lineterminator="\n").writerows(mat)

    print(f"per-genome -> {outdir}/per_genome_functional_qc.tsv ({len(rows)} genomes)")
    print(f"summaries  -> {outdir}/summary_by_cohort.tsv, summary_by_source.tsv, "
          f"candidate_functionality_by_cohort.tsv")
    for r in sorted(summarize(rows, "cohort"), key=lambda x: x["group"]):
        print(f"  {r['group']}: n={r['n_genomes']}  A={r['pct_tier_A']}%  "
              f"D={r['pct_tier_D']}%  %terminase={r['pct_with_terminase']}")
    if missing_checkv:
        print(f"  WARNING: {len(missing_checkv)} genomes missing CheckV results",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
