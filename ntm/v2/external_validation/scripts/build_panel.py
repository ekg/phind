#!/usr/bin/env python3
"""
build_panel.py — build the public mycobacteriophage reference panel from
cached source responses (offline, deterministic).

Reads   <data_dir>/cache/**  + fetch_manifest.json + taxonomy_hosts.json
Writes  <data_dir>/panel/:
          panel.full.fasta          (all retained records; unique panel IDs)
          manifest_panel.tsv        (one row per retained record)
          crosswalk.tsv             (aliases, duplicates, variants)
          exclusions.tsv            (explicit reasons)
          source_summary.tsv        (count reconciliation)
          host_distribution.tsv     (host-label / genus distribution)
          REPORT.md                 (coverage, caveats, licensing)
          PROVENANCE.md             (source URLs, timestamps, checksums)
        plus per-scope fasta: panel.mycobacteriophage.fasta,
          panel.other_actinobacteriophage.fasta

Re-running from the same cache reproduces every file byte-for-byte.

Usage:
  python3 build_panel.py [--data-dir DIR] [--out-dir DIR]
"""

from __future__ import annotations

import argparse
import collections
import glob
import gzip
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import panel_lib as pl  # noqa: E402

DEFAULT_DATA_DIR = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/external_validation"


def load_all(cache_dir: str):
    pdb_tsv = os.path.join(cache_dir, "phagesdb", "phagesdb_sequenced_phages_full.tsv")
    pdb_fa = os.path.join(cache_dir, "phagesdb", "Actinobacteriophages-All.fasta")
    inp_tsv = os.path.join(cache_dir, "inphared",
                           "7Apr2026_millardlab_website_table.txt.gz")
    inp_fa = os.path.join(cache_dir, "inphared", "7Apr2026_genomes.fa.gz")
    ncbi_dir = os.path.join(cache_dir, "ncbi")

    man = pl.read_fetch_manifest(cache_dir)

    print("[build] parsing phagesdb metadata + fasta …", flush=True)
    pdb_meta = pl.parse_phagesdb_tsv(pdb_tsv)
    pdb_fasta, pdb_extras = pl.parse_phagesdb_fasta(pdb_fa)
    print(f"[build] phagesdb: {len(pdb_meta)} metadata rows, "
          f"{len(pdb_fasta)} fasta records (+{len(pdb_extras)} dup/malformed)",
          flush=True)

    print("[build] parsing inphared table …", flush=True)
    inp_rows = pl.parse_inphared_table(inp_tsv)
    wanted = {r["accession"] for r in inp_rows if pl.inphared_scope(r)}
    print(f"[build] inphared: {len(inp_rows)} rows, {len(wanted)} in scope", flush=True)
    print(f"[build] streaming inphared genomes fasta ({len(wanted)} wanted) …",
          flush=True)
    inp_seqs = pl.load_inphared_seqs(inp_fa, wanted)
    print(f"[build] inphared sequences matched: {len(inp_seqs)}", flush=True)

    print("[build] parsing ncbi summaries + fasta …", flush=True)
    ss_files = sorted(glob.glob(os.path.join(ncbi_dir, "esummary_*.json")))
    fa_files = sorted(glob.glob(os.path.join(ncbi_dir, "efetch_fasta_*.txt")))
    bf_ss = sorted(glob.glob(os.path.join(ncbi_dir, "backfill_esummary_*.json")))
    bf_fa = sorted(glob.glob(os.path.join(ncbi_dir, "backfill_efetch_*.txt")))
    summ = pl.parse_ncbi_esummary(ss_files + bf_ss)
    ncbi_fa = pl.parse_ncbi_fasta(fa_files + bf_fa)
    print(f"[build] ncbi: {len(summ)} uids (incl {len(bf_ss)} backfill batches), "
          f"{len(ncbi_fa)} fasta records", flush=True)

    taxonomy = {}
    tpath = os.path.join(cache_dir, "taxonomy_hosts.json")
    if os.path.exists(tpath):
        taxonomy = json.load(open(tpath))
    gb = pl.parse_gb_remainder(os.path.join(ncbi_dir, "gb_remainder.txt"))

    return pl.merge_records(pdb_meta, pdb_fasta, pdb_extras, inp_rows, inp_seqs,
                            summ, ncbi_fa, taxonomy, gb, man), man


def build(data_dir: str, out_dir: str):
    cache_dir = os.path.join(data_dir, "cache")
    merged, man = load_all(cache_dir)
    records = merged["records"]

    if merged["needs_gb_remainder"]:
        sys.stderr.write(
            f"[build] {len(merged['needs_gb_remainder'])} NCBI-only records "
            "still need GenBank host metadata; run "
            "`fetch_sources.py --stage gb` first.\n"
            f"examples: {merged['needs_gb_remainder'][:5]}\n")
        raise SystemExit(2)

    os.makedirs(out_dir, exist_ok=True)
    pl.write_panel_fasta(os.path.join(out_dir, "panel.full.fasta"), records)
    by_scope = collections.defaultdict(list)
    for r in records:
        by_scope[r["scope"]].append(r)
    for scope, rs in by_scope.items():
        pl.write_panel_fasta(os.path.join(out_dir, f"panel.{scope}.fasta"), rs)

    pl.write_tsv(os.path.join(out_dir, "manifest_panel.tsv"),
                 pl.MANIFEST_COLUMNS, pl.manifest_rows(records))
    pl.write_tsv(os.path.join(out_dir, "crosswalk.tsv"),
                 pl.CROSSWALK_COLUMNS,
                 sorted(merged["crosswalk"], key=lambda r: (r["panel_id"], r["alias"])))
    pl.write_tsv(os.path.join(out_dir, "exclusions.tsv"),
                 pl.EXCLUSION_COLUMNS,
                 sorted(merged["exclusions"],
                        key=lambda r: (r["source"], r["reason"], r["identifier"])))

    # ---- count reconciliation ------------------------------------------
    sc = merged["source_counts"]
    bf_path = os.path.join(cache_dir, "ncbi", "backfill_found.json")
    n_backfill = len(json.load(open(bf_path))) if os.path.exists(bf_path) else 0
    rows = [
        {"metric": "phagesdb.metadata_rows", "count": sc["phagesdb"]["metadata_rows"]},
        {"metric": "phagesdb.fasta_records", "count": sc["phagesdb"]["fasta_records"]},
        {"metric": "phagesdb.retained_entries", "count": sc["phagesdb"]["retained"]},
        {"metric": "inphared.table_rows", "count": sc["inphared"]["table_rows"]},
        {"metric": "inphared.in_scope_rows", "count": sc["inphared"]["in_scope_rows"]},
        {"metric": "inphared.out_of_scope_rows", "count": sc["inphared"]["out_of_scope"]},
        {"metric": "inphared.with_sequence", "count": sc["inphared"]["with_sequence"]},
        {"metric": "ncbi.uids", "count": sc["ncbi"]["uids"]},
        {"metric": "ncbi.esearch_uids", "count": sc["ncbi"]["uids"] - n_backfill},
        {"metric": "ncbi.backfill_accessions_resolved", "count": n_backfill},
        {"metric": "ncbi.with_fasta", "count": sc["ncbi"]["with_fasta"]},
        {"metric": "ncbi.ncbi_only_retained", "count": sc["ncbi"]["ncbi_only_retained"]},
        {"metric": "entries.total", "count": sc["entries"]["total"]},
        {"metric": "entries.merged_duplicates", "count": sc["entries"]["merged_duplicates"]},
        {"metric": "entries.retained", "count": sc["entries"]["survivors"]},
        {"metric": "panel.records", "count": len(records)},
        {"metric": "panel.mycobacteriophage", "count": len(by_scope.get("mycobacteriophage", []))},
        {"metric": "panel.other_actinobacteriophage",
         "count": len(by_scope.get("other_actinobacteriophage", []))},
        {"metric": "panel.phage_host_unresolved",
         "count": len(by_scope.get("phage_host_unresolved", []))},
        {"metric": "exclusions.total", "count": len(merged["exclusions"])},
    ]
    for (src, reason), n in sorted(collections.Counter(
            (e["source"], e["reason"]) for e in merged["exclusions"]).items()):
        rows.append({"metric": f"exclusions.{src}.{reason}", "count": n})
    pl.write_tsv(os.path.join(out_dir, "source_summary.tsv"),
                 ["metric", "count"], rows)

    # ---- host distribution ----------------------------------------------
    host_counts = collections.Counter(
        (r["host_genus"] or "(unresolved)", r["scope"]) for r in records)
    host_rows = [{"host_genus": g, "scope": s, "records": n}
                 for (g, s), n in sorted(host_counts.items())]
    pl.write_tsv(os.path.join(out_dir, "host_distribution.tsv"),
                 ["host_genus", "scope", "records"], host_rows)

    write_report(out_dir, records, merged, man, by_scope, n_backfill)
    write_provenance(out_dir, man)
    print(f"[build] panel: {len(records)} records "
          f"({len(by_scope.get('mycobacteriophage', []))} mycobacteriophage, "
          f"{len(by_scope.get('other_actinobacteriophage', []))} other actino) "
          f"-> {out_dir}", flush=True)
    return records, merged


def write_report(out_dir, records, merged, man, by_scope, n_backfill=0):
    sc = merged["source_counts"]
    ev = collections.Counter(r["evidence_class"] for r in records)
    comp = collections.Counter(r["completeness"] or "(none)" for r in records)
    top = collections.Counter(r["host_genus"] or "(unresolved)" for r in records)
    src_combo = collections.Counter(r["sources"] for r in records)
    ncbi_primary = sum(1 for r in records if r["primary_source"] == "ncbi")
    versioned = sum(1 for r in records
                    if pl.is_versioned_accession(r["primary_accession"]))
    host_labels = sum(1 for r in records if r["host_reported"])
    host_unresolved = len(records) - host_labels

    lines = []
    a = lines.append
    a("# Public mycobacteriophage reference panel — build report")
    a("")
    a(f"Panel size: **{len(records)}** retained records "
      f"({len(by_scope.get('mycobacteriophage', []))} mycobacteriophages / "
      f"{len(by_scope.get('other_actinobacteriophage', []))} other "
      "actinobacteriophages from the PhagesDB leg).")
    a("")
    a("## Source coverage and reconciliation")
    a("")
    a("| source | metric | count |")
    a("|---|---|---:|")
    for m in [
        ("PhagesDB", "metadata rows (sequenced phages)", sc["phagesdb"]["metadata_rows"]),
        ("PhagesDB", "bulk FASTA records", sc["phagesdb"]["fasta_records"]),
        ("PhagesDB", "retained entries", sc["phagesdb"]["retained"]),
        ("INPHARED", "table rows (release 7Apr2026)", sc["inphared"]["table_rows"]),
        ("INPHARED", "in-scope mycobacteriophage rows", sc["inphared"]["in_scope_rows"]),
        ("INPHARED", "out-of-scope rows (non-myco hosts; by design)", sc["inphared"]["out_of_scope"]),
        ("INPHARED", "in-scope with sequence", sc["inphared"]["with_sequence"]),
        ("NCBI", "esearch uids", sc["ncbi"]["uids"] - n_backfill),
        ("NCBI", "backfill accessions resolved by id (esearch index lag)",
         n_backfill),
        ("NCBI", "records with FASTA", sc["ncbi"]["with_fasta"]),
        ("NCBI", "NCBI-only records retained after screening",
         sc["ncbi"]["ncbi_only_retained"]),
        ("merge", "accession-union entries", sc["entries"]["total"]),
        ("merge", "exact-sequence duplicates merged", sc["entries"]["merged_duplicates"]),
        ("merge", "retained panel records", sc["entries"]["survivors"]),
    ]:
        a(f"| {m[0]} | {m[1]} | {m[2]} |")
    a("")
    a(f"Exclusions with explicit reasons: {len(merged['exclusions'])} "
      "(see `exclusions.tsv`).")
    a("")
    a("NCBI versioned accessions used as primary identifiers for "
      f"**{ncbi_primary}** records; **{versioned}** of all primary accessions "
      "carry an explicit `.N` version. Records whose primary accession is a "
      "PhagesDB-only phage (no GenBank deposit) keep the PhagesDB name as "
      "identifier.")
    a("")
    a("## Source-combination distribution")
    a("")
    a("| sources | records |")
    a("|---|---:|")
    for k, n in src_combo.most_common():
        a(f"| `{k}` | {n} |")
    a("")
    a("## Host-label distribution (top genera)")
    a("")
    a("| host genus (reported) | records |")
    a("|---|---:|")
    for g, n in top.most_common(15):
        a(f"| {g} | {n} |")
    a("")
    a(f"Host labels resolved: {host_labels}; unresolved: {host_unresolved} "
      "(see completeness caveats).")
    a("")
    a("## Evidence class (cultured vs MAG/provirus)")
    a("")
    a("| evidence_class | records |")
    a("|---|---:|")
    for k, n in ev.most_common():
        a(f"| {k} | {n} |")
    a("")
    a("`isolated_sequenced` = PhagesDB plaque-purified, sequenced isolate. "
      "`deposited_isolate_sequence` = GenBank/INPHARED deposit without "
      "culture status in machine-readable form. `predicted_prophage` / "
      "`metagenome_assembled` are keyword-derived from source descriptions "
      "(heuristic — see caveats).")
    a("")
    a("## Completeness / topology")
    a("")
    a("| completeness label | records |")
    a("|---|---:|")
    for k, n in comp.most_common():
        a(f"| {k} | {n} |")
    a("")
    a("## Caveats")
    a("")
    a("1. **Reported host ≠ verified host range.** Every host field is a "
      "*reported label* (`host_evidence` column records its provenance: "
      "PhagesDB metadata, INPHARED table, or a GenBank `/host=` qualifier). "
      "None of these represent experimentally verified host range.")
    a("2. **RefSeq/GenBank dual membership.** Many phages carry both a "
      "GenBank and a RefSeq accession for the same sequence; these are "
      "merged as exact duplicates and both accessions remain in the "
      "crosswalk (`relation=same_sequence`).")
    a("3. **PhagesDB terminal overhangs.** PhagesDB bulk FASTA includes "
      "terminal overhang / extended ends; the GenBank (NCBI) version of the "
      "same phage can be shorter. Where both exist and bytes differ, the NCBI "
      "record is the canonical carrier and the PhagesDB variant is recorded "
      "in the crosswalk (`relation=sequence_variant`) with its own sha256.")
    a("4. **Exact-duplicate dedup only.** Dedup is by accession alias and "
      "by exact (uppercase) sequence hash; it is not reverse-complement or "
      "near-duplicate aware.")
    a("5. **Evidence-class heuristics.** MAG/provirus classification is "
      "keyword-derived from source descriptions where PhagesDB culture "
      "metadata is absent; treat as provisional.")
    a("6. **INPHARED scope.** Only Mycobacteriaceae-host (or "
      "mycobacteriophage-named) INPHARED rows are included; other actino "
      "hosts come from the PhagesDB leg only (task scope).")
    a("7. **Completeness labels are heterogeneous** (PhagesDB 'finished', "
      "NCBI title 'complete genome', INPHARED length-only) and are not "
      "comparable across sources without the `completeness_source` column.")
    a("")
    a("## Licensing / redistribution")
    a("")
    a("- **PhagesDB**: public research database; bulk downloads provided for "
      "research use (no explicit dataset license statement; site terms at "
      "phagesdb.org). Redistribute the cached raw FASTA with attribution.")
    a("- **INPHARED (Millard Lab, 7Apr2026 release)**: aggregates GenBank/"
      "RefSeq/ENA records plus PhagesDB; redistribution constraints follow "
      "the underlying repositories (mostly public-domain/CC0 for RefSeq; "
      "GenBank submitters retain rights but records are published openly).")
    a("- **NCBI GenBank/RefSeq**: U.S. Government public domain for RefSeq; "
      "GenBank records are openly available; submitters retain copyright on "
      "submitted sequences.")
    a("- The **panel FASTA is stored externally** at "
      "`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/external_validation/panel/`; "
      "the repository commits code, manifests, summaries and provenance only.")
    a("")
    a("## Determinism")
    a("")
    a("All outputs are pure functions of `cache/` (see PROVENANCE.md). "
      "Rebuilds from the same cache are byte-identical; "
      "`verify_panel.py --determinism` enforces this.")
    a("")
    open(os.path.join(out_dir, "REPORT.md"), "w").write("\n".join(lines))


def write_provenance(out_dir, man):
    lines = ["# Provenance — raw source artifacts", ""]
    a = lines.append
    a("| key | retrieved_at (UTC) | bytes | sha256 | url |")
    a("|---|---|---:|---|---|")
    for k in sorted(man):
        v = man[k]
        a(f"| `{k}` | {v.get('retrieved_at','')} | {v.get('bytes','')} | "
          f"`{v.get('sha256','')}` | {v.get('url','')} |")
    a("")
    a("Cached raw responses live under `cache/` (external storage). "
      "Per-batch NCBI efetch/esummary files are siblings of the recorded "
      "anchor file.")
    open(os.path.join(out_dir, "PROVENANCE.md"), "w").write("\n".join(lines))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=DEFAULT_DATA_DIR)
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args(argv)
    out = args.out_dir or os.path.join(args.data_dir, "panel")
    build(args.data_dir, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
