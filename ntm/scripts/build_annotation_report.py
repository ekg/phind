#!/usr/bin/env python3
"""
build_annotation_report.py — per-genome gene-annotation + QC report for all
phage genomes (E. coli ML, NTM ML, NTM ancestral), from Pharokka + CheckV.

Inputs:
  pharokka_out/pharokka_cds_final_merged_output.tsv  (per-CDS: contig, partial, phrog, annot, category)
  pharokka_out/pharokka_length_gc_cds_density.tsv     (per-contig: length, gc, cds density)
  checkv_out/quality_summary.tsv                      (completeness, contamination, viral/host genes)
  input/genome_index.tsv                              (source, clade_id)

Per-genome report -> annotation/report/per_genome_annotation_qc.tsv:
  gene_count, truncated(partial) gene count, key proteins present
  (terminase/portal/capsid/tail/integrase/holin-lysin), hypothetical fraction,
  PHROG category coverage, CheckV completeness/contamination, and a 'flag'
  column (missing core structural = terminase+portal+capsid all absent;
  low completeness; host contamination; many truncated).
"""
import csv
import os
import re
import statistics
from collections import defaultdict

BASE = "/mnt/nvme3n1/erikg/phind-genome-work/annotation"
CDS = f"{BASE}/pharokka_out/pharokka_cds_final_merged_output.tsv"
LGD = f"{BASE}/pharokka_out/pharokka_length_gc_cds_density.tsv"
CHECKV = f"{BASE}/checkv_out/quality_summary.tsv"
INDEX = f"{BASE}/input/genome_index.tsv"
OUT = f"{BASE}/report"
os.makedirs(OUT, exist_ok=True)

KEY = {
    "terminase": re.compile(r"terminase", re.I),
    "portal": re.compile(r"portal", re.I),
    "capsid": re.compile(r"capsid|major head|head protein", re.I),
    "tail": re.compile(r"tail|sheath|base ?plate|tail ?fiber|tail ?tube", re.I),
    "integrase": re.compile(r"integrase", re.I),
    "lysis": re.compile(r"holin|lysin|endolysin|spanin|R|lysis", re.I),
}
CAT2KEY = {"tail": "tail", "integration and excision": "integrase", "lysis": "lysis"}


def main():
    meta = {}
    with open(INDEX) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            meta[r["genome_id"]] = r

    lgd = {}
    with open(LGD) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            lgd[r["contig"]] = r

    genes = defaultdict(int)
    partial = defaultdict(int)
    keys = defaultdict(set)
    cats = defaultdict(set)
    hypoth = defaultdict(int)
    with open(CDS) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            gid = r["contig"]
            if gid not in meta:
                continue
            genes[gid] += 1
            if (r.get("partial") or "00").strip() != "00":
                partial[gid] += 1
            annot = r.get("annot", "") or ""
            for k, rx in KEY.items():
                if rx.search(annot):
                    keys[gid].add(k)
            cat = (r.get("category", "") or "").strip().lower()
            if cat:
                cats[gid].add(cat)
                if cat in CAT2KEY:
                    keys[gid].add(CAT2KEY[cat])
            if cat == "unknown function" or "hypothetical" in annot.lower():
                hypoth[gid] += 1

    cv = {}
    if os.path.exists(CHECKV):
        with open(CHECKV) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                cv[r["contig_id"]] = r

    KEYLIST = ["terminase", "portal", "capsid", "tail", "integrase", "lysis"]
    rows = []
    for gid, m in meta.items():
        g = genes.get(gid, 0)
        ks = keys.get(gid, set())
        lg = lgd.get(gid, {})
        c = cv.get(gid, {})
        comp = (c.get("checkv_quality") or "").strip()
        comp_pct = (c.get("completeness") or "").strip()
        cont = (c.get("contamination") or "").strip()
        vgenes = (c.get("viral_genes") or "").strip()
        hgenes = (c.get("host_genes") or "").strip()
        flags = []
        if g == 0:
            flags.append("no_genes")
        if g > 0 and not ({"terminase", "portal", "capsid"} & ks):
            flags.append("missing_core_structural")
        if any(x in comp.lower() for x in ("low quality", "low-quality", "not determined", "genome fragment")):
            flags.append("low_completeness")
        try:
            if cont and float(cont) > 0:
                flags.append("contamination")
        except ValueError:
            pass
        try:
            if hgenes and int(hgenes) > 0:
                flags.append("host_genes")
        except ValueError:
            pass
        if g > 0 and partial.get(gid, 0) / g > 0.2:
            flags.append("many_truncated")
        rows.append({
            "genome_id": gid, "source": m["source"], "clade_id": m["clade_id"],
            "length": m["length"], "gc_perc": lg.get("gc_perc", ""),
            "cds_density": lg.get("cds_coding_density", ""),
            "gene_count": g, "truncated_genes": partial.get(gid, 0),
            "terminase": int("terminase" in ks), "portal": int("portal" in ks),
            "capsid": int("capsid" in ks), "tail": int("tail" in ks),
            "integrase": int("integrase" in ks), "lysis": int("lysis" in ks),
            "n_key_proteins": len(ks & set(KEYLIST)),
            "has_head_packaging": int("head and packaging" in cats.get(gid, set())),
            "has_tail_cat": int("tail" in cats.get(gid, set())),
            "has_lysis_cat": int("lysis" in cats.get(gid, set())),
            "has_integration_cat": int("integration and excision" in cats.get(gid, set())),
            "hypothetical": hypoth.get(gid, 0),
            "hypothetical_frac": round(hypoth.get(gid, 0) / g, 3) if g else "",
            "checkv_miuvig": comp, "completeness_pct": comp_pct,
            "contamination": cont, "viral_genes": vgenes, "host_genes": hgenes,
            "flag": ";".join(flags) or "ok",
        })

    cols = list(rows[0].keys())
    with open(f"{OUT}/per_genome_annotation_qc.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    by_src = defaultdict(list)
    for r in rows:
        by_src[r["source"]].append(r)
    with open(f"{OUT}/summary_by_source.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["source", "n_genomes", "median_genes", "median_key_proteins",
                    "pct_with_terminase", "pct_with_capsid", "pct_missing_core",
                    "pct_low_completeness", "pct_contaminated", "median_truncated_frac",
                    "median_hypothetical_frac"])
        for src, rs in sorted(by_src.items()):
            n = len(rs)
            def med(key):
                return statistics.median(r[key] for r in rs)
            def pct(f):
                return round(100 * sum(1 for r in rs if f(r)) / n, 1)
            w.writerow([src, n, int(med("gene_count")), med("n_key_proteins"),
                        pct(lambda r: r["terminase"]), pct(lambda r: r["capsid"]),
                        pct(lambda r: "missing_core_structural" in r["flag"]),
                        pct(lambda r: "low_completeness" in r["flag"]),
                        pct(lambda r: "contamination" in r["flag"]),
                        round(statistics.median(r["truncated_genes"] / r["gene_count"] for r in rs if r["gene_count"]), 3),
                        round(statistics.median(r["hypothetical_frac"] for r in rs if r["hypothetical_frac"] != ""), 3)])

    print(f"report -> {OUT}/per_genome_annotation_qc.tsv ({len(rows)} genomes)")
    print(f"summary -> {OUT}/summary_by_source.tsv")
    for src, rs in sorted(by_src.items()):
        n = len(rs)
        t = sum(r["terminase"] for r in rs)
        mc = sum("missing_core_structural" in r["flag"] for r in rs)
        print(f"  {src}: n={n}  %terminase={100*t//n}%  %missing_core={100*mc//n}%")


if __name__ == "__main__":
    main()
