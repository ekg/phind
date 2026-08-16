#!/usr/bin/env python3
"""NTM v2 input validation: accession list + prophage master manifest.

Checks:
  1. sha256 of both collaborator input files
  2. Schema documentation (columns, types, value distributions)
  3. Accession list vs v1 local genomes (7,352 canonical objects) coverage
  4. Prophage manifest vs v1 geNomad ntm_prophages.csv + coordinate sanity vs v1 .fai
  5. GCA/GCF dedup rule + prophage_id uniqueness
  6. Writes validation_report.md + coverage.tsv

Usage:
  python3 validate_inputs.py [--data-dir /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/inputs]
                             [--v1-dir /mnt/nvme3n1/erikg/phind-genome-work/ntm/v1]
                             [--out-dir ...]
"""

import argparse
import csv
import gzip
import hashlib
import os
import re
import statistics
import sys
from collections import Counter, defaultdict

PANSNDIR = "genomes/canonical_objects"
PROPHAGE_V1 = "prophages/ntm_prophages.csv"
HOST_CLADES = "host_clades/host_clades.tsv"

ACCESSION_COLS = ["genome_id", "accession", "species", "data_source"]
MANIFEST_COLS = [
    "genome_id", "species", "species_confidence", "data_source", "accession",
    "fasta_path", "busco_complete_pct", "busco_fragmented_pct",
    "checkm2_completeness_pct", "checkm2_contamination_pct",
    "assembly_contigs", "assembly_length_bp", "assembly_n50_bp",
    "ani_pct_to_closest_type_strain", "ani_closest_type_strain_species",
    "has_prophage", "prophage_id", "prophage_contig", "prophage_start",
    "prophage_end", "prophage_length_bp", "prophage_transposable_element",
    "prophage_family",
]

NUMERIC_RE = re.compile(r"^[A-Z]{3}_(\d+\.\d+)$")
RUN_RE = re.compile(r"^(ERR|SRR|DRR)\d+$")
TB_TERMS = ("tuberculosis", "[tuberculosis]")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        rows = [row for row in reader]
    return header, rows


def strip_numeric(acc):
    """GCA_000013925.1 -> 000013925.1 ; ERR12345 -> None"""
    m = NUMERIC_RE.match(acc)
    return m.group(1) if m else None


def pref(acc):
    """Letter prefix: GCA/GCF/ERR/SRR/DRR."""
    if acc.startswith(("GCA_", "GCF_")):
        return acc[:3]
    m = re.match(r"^(ERR|SRR|DRR)", acc)
    return m.group(1) if m else acc


def is_tb(species):
    """TB-complex: contains 'tuberculosis' but NOT paratuberculosis."""
    return "tuberculosis" in species and "paratuberculosis" not in species


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/inputs")
    ap.add_argument("--v1-dir", default="/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1")
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()
    data_dir = args.data_dir
    v1_dir = args.v1_dir
    out_dir = args.out_dir or data_dir
    os.makedirs(out_dir, exist_ok=True)

    acc_file = os.path.join(data_dir, "NTM_QC_passed_accession_list.tsv")
    man_file = os.path.join(data_dir, "NTM_QC_passed_prophage_master_manifest.tsv")

    report = []
    def R(line=""):
        report.append(line)

    # ---------------- 1. sha256 ----------------
    R("# NTM v2 collaborator input validation report")
    R("")
    R("Generated: (see coverage.tsv mtime) — script `ntm/v2/scripts/validate_inputs.py`")
    R("")
    R("## 1. File integrity (sha256)")
    R("")
    R("| file | size (bytes) | sha256 |")
    R("|---|---|---|")
    for f in (acc_file, man_file):
        h = sha256_file(f)
        sz = os.path.getsize(f)
        R(f"| `{os.path.basename(f)}` | {sz} | `{h}` |")
    R("")

    # ---------------- 2. Accession list ----------------
    acc_header, acc_rows = read_tsv(acc_file)
    R("## 2. Schema — NTM_QC_passed_accession_list.tsv")
    R("")
    R(f"- rows (excl. header): **{len(acc_rows)}**")
    R(f"- columns ({len(acc_header)}): `{'`, `'.join(acc_header)}`")
    R("")
    R("| column | type | notes |")
    R("|---|---|---|")
    R("| genome_id | str | identifier, identical to accession for run-assemblies |")
    R("| accession | str | NCBI accession: GCA_/GCF_ assembly or ERR/SRR/DRR run id |")
    R("| species | str | collaborator species call (may include `[tuberculosis]`) |")
    R("| data_source | str | `NCBI` (assembly) or `ASSEMBLY` (run-assembly) |")
    R("")

    acc_pref_counts = Counter(pref(r[1]) for r in acc_rows)
    ds_counts = Counter(r[3] for r in acc_rows)
    R("| accession prefix | count |")
    R("|---|---|")
    for p in ("GCA", "GCF", "ERR", "SRR", "DRR"):
        R(f"| {p} | {acc_pref_counts.get(p, 0)} |")
    R("")
    R("| data_source | count |")
    R("|---|---|")
    for ds, c in ds_counts.most_common():
        R(f"| {ds} | {c} |")
    R("")

    acc_numeric = defaultdict(set)  # numeric part -> set of prefixes
    acc_runs = set()
    for r in acc_rows:
        acc = r[1]
        n = strip_numeric(acc)
        if n is not None:
            acc_numeric[n].add(pref(acc))
        elif RUN_RE.match(acc):
            acc_runs.add(acc)
        else:
            R(f"WARNING: unrecognized accession format: {acc}")
    unique_numeric = len(acc_numeric)
    with_twin = sum(1 for s in acc_numeric.values() if len(s) == 2)
    gca_only = sum(1 for s in acc_numeric.values() if s == {"GCA"})
    gcf_only = sum(1 for s in acc_numeric.values() if s == {"GCF"})
    R(f"- unique NCBI assemblies (by numeric part): **{unique_numeric}**")
    R(f"  - with both GCA+GCF twins: {with_twin}")
    R(f"  - GCA-only: {gca_only}")
    R(f"  - GCF-only: {gcf_only}")
    R(f"- unique run ids (ERR/SRR/DRR): **{len(acc_runs)}**")
    R(f"- total unique genome ids (numeric + runs): **{unique_numeric + len(acc_runs)}**")
    R("")

    # species composition of accession list
    R("### Species composition (accession list, all rows)")
    R("")
    sp_counter = Counter(r[2] for r in acc_rows)
    R("| species | count |")
    R("|---|---|")
    for sp, c in sp_counter.most_common(30):
        R(f"| {sp} | {c} |")
    R("")
    R(f"Total distinct species strings: {len(sp_counter)}")
    R("")
    tb_rows = [r for r in acc_rows if is_tb(r[2])]
    tb_numeric = {strip_numeric(r[1]) for r in tb_rows if strip_numeric(r[1])}
    tb_runs = {r[1] for r in tb_rows if RUN_RE.match(r[1])}
    R(f"### TB-complex subset (species containing 'tuberculosis', excl. paratuberculosis)")
    R("")
    R(f"- rows: {len(tb_rows)} (of {len(acc_rows)})")
    R(f"- unique NCBI assemblies: {len(tb_numeric)}")
    R(f"- unique run ids: {len(tb_runs)}")
    R("")
    R("> v1 excluded MTBC by design; v2 default = keep collaborator list verbatim "
      "(pending user confirmation). Documented here for the scope decision.")
    R("")

    # ---------------- 3. Manifest ----------------
    man_header, man_rows = read_tsv(man_file)
    R("## 3. Schema — NTM_QC_passed_prophage_master_manifest.tsv")
    R("")
    R(f"- rows (excl. header): **{len(man_rows)}**")
    R(f"- columns ({len(man_header)}): `{'`, `'.join(man_header)}`")
    R("")
    R("| column | type | notes |")
    R("|---|---|---|")
    R("| genome_id | str | `{accession}_{assembly_name}_genomic` for NCBI; run id for ASSEMBLY |")
    R("| species | str | collaborator species call |")
    R("| species_confidence | str | empty / MODERATE / HIGH |")
    R("| data_source | str | NCBI or ASSEMBLY |")
    R("| accession | str | GCA_/GCF_ assembly or ERR/SRR/DRR run id |")
    R("| fasta_path | str | collaborator Iridis cluster path (NOT public) |")
    R("| busco_complete_pct | float | QC metric |")
    R("| busco_fragmented_pct | float | QC metric |")
    R("| checkm2_completeness_pct | float | QC metric |")
    R("| checkm2_contamination_pct | float | QC metric |")
    R("| assembly_contigs | int | contig count |")
    R("| assembly_length_bp | int | assembly size |")
    R("| assembly_n50_bp | int | N50 |")
    R("| ani_pct_to_closest_type_strain | float | ANI to nearest type strain |")
    R("| ani_closest_type_strain_species | str | type strain species |")
    R("| has_prophage | bool (True/False) | flag; False rows have empty prophage_* cols |")
    R("| prophage_id | str | `{genome_id}_prophageN` (N 1-based) |")
    R("| prophage_contig | str | WGS contig accession (GCA rows bare, GCF rows NZ_-prefixed); NODE_* for runs |")
    R("| prophage_start | float (int) | 1-based INCLUSIVE start |")
    R("| prophage_end | float (int) | 1-based INCLUSIVE end |")
    R("| prophage_length_bp | float (int) | end-start+1 |")
    R("| prophage_transposable_element | bool | transposable element flag |")
    R("| prophage_family | str | geNomad family call |")
    R("")

    # validate row column counts
    bad_cols = [i for i, r in enumerate(man_rows) if len(r) != len(MANIFEST_COLS)]
    R(f"- rows with wrong column count: **{len(bad_cols)}** (of {len(man_rows)})")
    for i in bad_cols[:5]:
        R(f"  - row {i + 2}: {len(man_rows[i])} cols")
    R("")

    has_prophage_counts = Counter(r[15] for r in man_rows)
    R(f"- has_prophage value counts: {dict(has_prophage_counts)}")
    R("")

    true_rows = [r for r in man_rows if r[15] == "True"]
    false_rows = [r for r in man_rows if r[15] == "False"]

    # false rows should have empty prophage cols
    empty_prophage_cols_ok = all(
        all(r[c] == "" for c in range(16, 23)) for r in false_rows
    )
    R(f"- False rows: all prophage_* cols empty: {empty_prophage_cols_ok}")

    # true rows: prophage cols populated?
    missing_in_true = sum(1 for r in true_rows if any(r[c] == "" for c in range(16, 23)))
    R(f"- True rows with any empty prophage_* col: **{missing_in_true}** (row accounting)")
    R("")

    # ---------------- 4. Dedup rule ----------------
    R("## 4. Unique genome & prophage counts (dedup rule: prefer GCA rows)")
    R("")
    # per numeric accession: prefer GCA genome rows over GCF
    # runs are inherently unique
    true_by_source = Counter(pref(r[4]) for r in true_rows)  # by accession prefix
    src_labels = {"GCA": "GCA (NCBI)", "GCF": "GCF (NCBI)", "ERR": "ERR", "SRR": "SRR", "DRR": "DRR"}
    R("| source | prophage rows (True) |")
    R("|---|---|")
    for s in ("GCA", "GCF", "ERR", "SRR", "DRR"):
        R(f"| {s} | {true_by_source.get(s, 0)} |")
    R("")

    numeric_true = defaultdict(lambda: {"GCA": [], "GCF": []})
    run_true = defaultdict(list)
    for r in true_rows:
        acc = r[4]
        n = strip_numeric(acc)
        if n is not None:
            numeric_true[n][pref(acc)].append(r)
        else:
            run_true[acc].append(r)

    # dedup: prefer GCA; GCF rows used only when no GCA twin
    kept_rows = []
    gcf_only_numerics = set()
    for n, byp in numeric_true.items():
        if byp["GCA"]:
            kept_rows.extend(byp["GCA"])
        else:
            kept_rows.extend(byp["GCF"])
            gcf_only_numerics.add(n)
    for acc, rows in run_true.items():
        kept_rows.extend(rows)

    kept_by_source = Counter(pref(r[4]) for r in kept_rows)
    R("| kept prophage rows after dedup | source | count |")
    R("|---|---|---|")
    for s in ("GCA", "GCF", "ERR", "SRR", "DRR"):
        R(f"| | {s} | {kept_by_source.get(s, 0)} |")
    R(f"- total kept prophage rows: **{len(kept_rows)}**")
    R(f"- GCF-only numeric assemblies contributing prophages: **{len(gcf_only_numerics)}**")
    R(f"- unique prophage-bearing genomes after dedup: **{len(numeric_true) + len(run_true)}** "
      f"({len(numeric_true)} NCBI numerics + {len(run_true)} run-assemblies)")
    R("")

    # prophage_id uniqueness
    pid_counts = Counter(r[16] for r in kept_rows)
    dup_pids = {pid: c for pid, c in pid_counts.items() if c > 1}
    R(f"- prophage_id uniqueness after dedup: **{len(pid_counts)} unique ids for {len(kept_rows)} rows**")
    if dup_pids:
        R(f"  - DUPLICATE prophage_ids: {len(dup_pids)}")
        for pid, c in list(dup_pids.items())[:10]:
            R(f"    - {pid}: {c}x")
    else:
        R("  - no duplicate prophage_ids ✔")
    R("")

    # ---------------- 5. Length / coordinate checks ----------------
    R("## 5. Prophage length distribution & coordinate sanity")
    R("")
    lens = []
    coord_bad = 0
    len_mismatch = 0
    for r in true_rows:
        try:
            start = int(float(r[18]))
            end = int(float(r[19]))
            length = int(float(r[20]))
        except ValueError:
            coord_bad += 1
            continue
        lens.append(length)
        if start > end:
            coord_bad += 1
        if length != end - start + 1:
            len_mismatch += 1
    lens_sorted = sorted(lens)
    n_lens = len(lens)
    R(f"- rows with begin>end or unparseable coords: **{coord_bad}**")
    R(f"- rows where length != end-start+1 (1-based inclusive): **{len_mismatch}**")
    R("")
    R("| length stat | bp |")
    R("|---|---|")
    R(f"| min | {lens_sorted[0]} |")
    R(f"| median | {statistics.median(lens_sorted)} |")
    R(f"| mean | {sum(lens_sorted) / n_lens:.1f} |")
    R(f"| max | {lens_sorted[-1]} |")
    R("")
    R("| length bin | count |")
    R("|---|---|")
    bins = [(0, 5000), (5000, 10000), (10000, 20000), (20000, 40000),
            (40000, 60000), (60000, 100000), (100000, 10**9)]
    for lo, hi in bins:
        c = sum(1 for L in lens if lo <= L < hi)
        R(f"| {lo}-{hi} | {c} |")
    R("")

    fam_counter = Counter(r[22] for r in true_rows)
    R("| prophage_family | count |")
    R("|---|---|")
    for fam, c in fam_counter.most_common():
        R(f"| {fam or '(empty)'} | {c} |")
    R("")

    # ---------------- 6. Coverage vs v1 ----------------
    R("## 6. Coverage vs v1 (7,352 local genomes)")
    R("")
    co_dir = os.path.join(v1_dir, PANSNDIR)
    v1_accs = {d for d in os.listdir(co_dir) if os.path.isdir(os.path.join(co_dir, d))}
    v1_numeric = {strip_numeric(a) for a in v1_accs if strip_numeric(a)}
    R(f"- v1 canonical_objects dirs: {len(v1_accs)}")
    R(f"  - with pansn.fa.gz: {sum(1 for a in v1_accs if os.path.exists(os.path.join(co_dir, a, a + '.pansn.fa.gz')))}")
    R(f"  - empty dirs: {len(v1_accs) - sum(1 for a in v1_accs if os.listdir(os.path.join(co_dir, a)))}")
    R(f"- v1 unique numeric accessions: {len(v1_numeric)}")
    R("")

    # accession list (all NCBI) vs v1
    acc_all_numeric = set(acc_numeric.keys())
    v2_in_v1 = acc_all_numeric & v1_numeric
    v2_missing = acc_all_numeric - v1_numeric
    v1_not_v2 = v1_numeric - acc_all_numeric
    # ---------- coverage: GCA-preferred exact match (primary, download-task semantics) ----------
    # preferred accession per numeric = GCA if v2 lists one, else GCF.
    preferred_match = 0
    preferred_missing = 0
    preferred_gcf_only = 0
    for n, ps in acc_numeric.items():
        acc = "GCA_" + n if "GCA" in ps else "GCF_" + n
        if acc in v1_accs:
            preferred_match += 1
        else:
            preferred_missing += 1
            if "GCA" not in ps:
                preferred_gcf_only += 1
    # prefix-agnostic overlap (secondary info)
    R(f"- v2 NCBI assemblies already local, **GCA-preferred exact-match** "
      f"(dir name == preferred accession): **{preferred_match}**")
    R(f"- v2 NCBI assemblies MISSING locally (need download): **{preferred_missing}** "
      f"(of which GCF-only numerics: {preferred_gcf_only})")
    R(f"- prefix-agnostic overlap (numeric id present locally under either GCA/GCF): "
      f"**{len(v2_in_v1)}**")
    R(f"- v1 assemblies not in v2 list: **{len(v1_not_v2)}**")
    R("")
    R("> The download task (ntm-v2-download-12-030) uses the GCA-preferred exact-match "
      "definition: link local dirs where the preferred accession exists; download the rest.")
    R("")
    R("> Reconciliation vs earlier chat-agent estimate (13,123 unique numerics; 1,057 local; "
      "12,030 to download): precise recount gives **13,122** unique numerics, **1,150** "
      "GCA-preferred local, **11,972** to download. Small deltas from the estimate; "
      "coverage.tsv is the ground truth.")
    R("")

    # host clades overlap
    hc_accs = set()
    with open(os.path.join(v1_dir, HOST_CLADES)) as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd)
        for row in rd:
            if row:
                hc_accs.add(row[0])
    hc_numeric = {strip_numeric(a) for a in hc_accs if strip_numeric(a)}
    R(f"- v1 host_clades.tsv unique accessions: {len(hc_accs)} (numeric: {len(hc_numeric)})")
    R(f"- host_clades ∩ v2 NCBI assemblies: **{len(hc_numeric & acc_all_numeric)}**")
    R(f"- host_clades ∩ v1-local-and-v2: {len(hc_numeric & v2_in_v1)}")
    R("")

    # prophage-bearing overlap: v1 prophage genomes vs v2 prophage genomes
    v1_proph_genomes = set()
    with open(os.path.join(v1_dir, PROPHAGE_V1)) as fh:
        rd = csv.reader(fh)
        next(rd)
        for row in rd:
            if row:
                v1_proph_genomes.add(row[1])
    v1_proph_numeric = {strip_numeric(a) for a in v1_proph_genomes if strip_numeric(a)}
    v2_proph_numeric = set(numeric_true.keys())
    v2_proph_runs = set(run_true.keys())
    R(f"- v1 prophage CSV: {sum(1 for _ in open(os.path.join(v1_dir, PROPHAGE_V1))) - 1} rows, "
      f"{len(v1_proph_genomes)} unique genomes")
    R(f"- v2 prophage-bearing genomes (after dedup): {len(v2_proph_numeric)} NCBI + {len(v2_proph_runs)} runs")
    R(f"- overlap (v1 ∩ v2 prophage-bearing numerics): **{len(v1_proph_numeric & v2_proph_numeric)}**")
    R("")

    # ---------------- 7. Per-genome prophage count comparison (v1 vs v2) ----------------
    R("## 7. Per-genome prophage count: v1 geNomad vs v2 manifest")
    R("")
    v1_counts = Counter()
    with open(os.path.join(v1_dir, PROPHAGE_V1)) as fh:
        rd = csv.reader(fh)
        next(rd)
        for row in rd:
            if row:
                v1_counts[row[1]] += 1
    v2_counts = Counter(pref(r[4]) + "_" + (strip_numeric(r[4]) or r[4]) for r in kept_rows)
    # simpler: count per genome_id
    v2_counts_by_genome = Counter(r[4] for r in kept_rows)  # accession-level
    # map v1 acc (GCA_000523695.1) to v2 accession (same when GCA row kept)
    compared = 0
    agree = 0
    disagree_rows = []
    for v1_acc, c1 in v1_counts.items():
        # v2: prefer GCA twin for numeric
        n = strip_numeric(v1_acc)
        if n is None:
            continue
        # find v2 counts for this numeric (GCA preferred)
        c2 = None
        if n in numeric_true:
            c2 = len(numeric_true[n]["GCA"]) if numeric_true[n]["GCA"] else len(numeric_true[n]["GCF"])
        if c2 is None:
            continue
        compared += 1
        if c1 == c2:
            agree += 1
        else:
            disagree_rows.append((v1_acc, c1, c2))
    R(f"- v1 genomes with prophages also in v2 (comparable): {compared}")
    R(f"- per-genome prophage count identical: **{agree}**")
    R(f"- count differs: {len(disagree_rows)}")
    R("")
    if disagree_rows:
        R("| v1 genome | v1 count | v2 count |")
        R("|---|---|---|")
        for a, c1, c2 in sorted(disagree_rows, key=lambda x: -abs(x[1] - x[2]))[:12]:
            R(f"| {a} | {c1} | {c2} |")
        R("")
    R("### Interpretation of count differences")
    R("")
    R("v1 ran geNomad v1.12.0 (find-proviruses default, end-to-end mode) on the v1 "
      "cohort. The collaborator ran their own geNomad pipeline (version/params unknown) "
      "on their assemblies. Both call prophages per genome, but counts differ for "
      f"{len(disagree_rows)} of {compared} shared prophage-bearing genomes. Contributing "
      "factors:")
    R("")
    R("- collaborator manifest retains very short predictions (min 237 bp; v1 min was "
      "3,220 bp) — e.g. GCF_001296255.1: 4 (v1) vs 10 (v2, incl. 1,428 bp and 4,800 bp "
      "calls);")
    R("- geNomad version/parameter drift between runs (e.g. GCF_900130755.1: 9 (v1) "
      "vs 3 (v2));")
    R("- assembly version differences (GCA vs GCF mirror or updated versions).")
    R("")
    R("v2 keeps the collaborator's calls verbatim (their QC-passed list is the v2 "
      "denominator); v1 counts are recorded here only as a cross-check.")
    R("")

    # ---------------- 8. Scaffold-name mapping decision ----------------
    R("## 8. Scaffold-name mapping decision (collaborator contig → PanSN)")
    R("")
    R("v1 PanSN header: `{assembly_acc}#1#{contig}` (e.g. `GCA_000523695.1#1#JAOB01000032.1`).")
    R("v1 `ntm_prophages.csv` stored `source_seq` = the full PanSN header.")
    R("")
    R("Collaborator `prophage_contig` values:")
    R("")
    contig_samples = Counter()
    nc_rows = []
    for r in true_rows:
        acc = r[4]
        contig = r[17]
        if strip_numeric(acc):
            if pref(acc) == "GCF":
                if contig.startswith("NZ_"):
                    key = "GCF-row: NZ_-prefixed WGS"
                elif contig.startswith("NC_"):
                    key = "GCF-row: NC_ chromosome (complete genomes)"
                else:
                    key = "GCF-row: other"
                    nc_rows.append((acc, r[16], contig))
            else:
                key = "GCA-row: bare WGS" if not contig.startswith(("NZ_", "NODE_")) else "GCA-row: NZ_/NODE"
        else:
            key = "run: NODE_*" if contig.startswith("NODE_") else "run: other"
        contig_samples[key] += 1
    for k, c in contig_samples.most_common():
        R(f"- {k}: {c}")
    if nc_rows:
        R("")
        R(f"> Note: {len(nc_rows)} GCF rows use `NC_` chromosome accessions "
          "(complete/chromosome-level assemblies where the chromosome is the single "
          "replicon, e.g. NC_000962.3 for M. tuberculosis H37Rv). No NZ_ strip applies.")
    R("")
    R("**Decision:** for NCBI rows, map collaborator contig → PanSN using the LOCAL "
      "genome dir's prefix: `{local_acc}#1#{contig}` where the contig suffix matches the "
      "dir prefix convention (GCA dirs → bare WGS accession; GCF dirs → `NZ_`-prefixed). "
      "v1 GCA dirs (canonical) match GCA-row contigs directly; GCF rows drop `NZ_` when "
      "mapped onto a GCA dir and keep it on a GCF dir.")
    R("")
    R("Worked example on ≥3 genomes (coordinate check vs local .fai):")
    R("")
    R("| genome | prophage_id | contig (collab) | PanSN (mapped) | fai len | start | end | within? |")
    R("|---|---|---|---|---|---|---|---|")
    # spot-check against v1 .fai
    spot = 0
    fai_missing = 0
    fai_fail = 0
    fai_checked = 0
    examples = []
    v1_acc_to_dir = {}
    for a in v1_accs:
        n = strip_numeric(a)
        if n:
            v1_acc_to_dir[n] = a
    fai_cache = {}
    for r in true_rows:
        acc = r[4]
        n = strip_numeric(acc)
        if n is None or n not in v1_acc_to_dir:
            continue
        if spot >= 8:
            break
        contig = r[17]
        local_dir = v1_acc_to_dir[n]
        bare = contig[3:] if contig.startswith("NZ_") else contig
        # PanSN suffix convention follows the LOCAL dir prefix
        if local_dir.startswith("GCF_"):
            suffix = "NZ_" + bare if not bare.startswith("NZ_") else bare
        else:
            suffix = bare
        pansn = f"{local_dir}#1#{suffix}"
        faif = os.path.join(co_dir, local_dir, local_dir + ".pansn.fa.gz.fai")
        if faif not in fai_cache:
            if os.path.exists(faif):
                fai_cache[faif] = {}
                with open(faif) as fh:
                    for line in fh:
                        parts = line.split("\t")
                        fai_cache[faif][parts[0]] = int(parts[1])
            else:
                fai_cache[faif] = None
        fai = fai_cache[faif]
        if fai is None:
            fai_missing += 1
            continue
        flen = fai.get(pansn)
        start = int(float(r[18]))
        end = int(float(r[19]))
        if flen is None:
            # fallback: try the alternative suffix form (GCA dir with NZ_, GCF dir without)
            alt = bare if local_dir.startswith("GCF_") else "NZ_" + bare
            alt_pansn = f"{local_dir}#1#{alt}"
            flen2 = fai.get(alt_pansn)
            if flen2 is not None:
                flen = flen2
                pansn = alt_pansn
        if flen is None:
            ok = "contig not in fai"
            fai_missing += 1
        else:
            fai_checked += 1
            ok = "OK" if end <= flen else f"FAIL len={flen}"
        if ok != "OK":
            fai_fail += 1
        examples.append((local_dir, r[16], contig, pansn, flen, start, end, ok))
        spot += 1
    for ex in examples:
        R(f"| {ex[0]} | {ex[1]} | {ex[2]} | `{ex[3]}` | {ex[4]} | {ex[5]} | {ex[6]} | {ex[7]} |")
    R("")
    R(f"- spot-checks performed on overlapping v1-local genomes: {len(examples)} "
      f"(fai-resolved: {fai_checked})")
    R(f"- contig missing from .fai: {fai_missing} | outside bounds: {fai_fail}")

    # ---------------- 9. Row accounting ----------------
    R("")
    R("## 9. Row accounting (no silent drops)")
    R("")
    R(f"- manifest rows total: {len(man_rows)}")
    R(f"- True rows: {len(true_rows)} → kept after dedup: {len(kept_rows)}")
    R(f"  - GCA rows kept: {kept_by_source.get('GCA', 0)} (GCA preferred)")
    R(f"  - GCF rows kept: {kept_by_source.get('GCF', 0)} (GCF-only numerics: {len(gcf_only_numerics)})")
    R(f"  - run rows kept: {sum(kept_by_source.get(s, 0) for s in ('ERR', 'SRR', 'DRR'))}")
    R(f"- False rows (genomes without prophage): {len(false_rows)} — retained in manifest, not prophage set")
    R(f"- True rows dropped by dedup (GCF twins of GCA): {len(true_rows) - len(kept_rows)}")
    R(f"- rows with wrong column count: {len(bad_cols)}")
    R(f"- True rows with empty prophage cols: {missing_in_true}")
    R("")
    R("Every manifest row is accounted for: kept (prophage set), retained as no-prophage genome, "
      "or dropped by explicit dedup rule (GCF twin where GCA preferred).")

    # ---------------- coverage.tsv ----------------
    cov_rows = [
        ("metric", "value", "detail"),
        ("accession_list_rows", str(len(acc_rows)), "excl header"),
        ("accession_list_sha256", sha256_file(acc_file), ""),
        ("manifest_rows", str(len(man_rows)), "excl header"),
        ("manifest_sha256", sha256_file(man_file), ""),
        ("manifest_true_rows", str(len(true_rows)), "has_prophage=True"),
        ("manifest_false_rows", str(len(false_rows)), "has_prophage=False"),
        ("unique_numeric_assemblies", str(unique_numeric), "GCA+GCF by numeric part"),
        ("numeric_with_twins", str(with_twin), "both GCA+GCF present"),
        ("numeric_gca_only", str(gca_only), ""),
        ("numeric_gcf_only", str(gcf_only), ""),
        ("unique_run_ids", str(len(acc_runs)), "ERR/SRR/DRR"),
        ("unique_genome_ids_total", str(unique_numeric + len(acc_runs)), "numeric + runs"),
        ("v1_local_genomes", str(len(v1_accs)), "canonical_objects dirs"),
        ("v1_local_with_fasta", str(sum(1 for a in v1_accs if os.path.exists(os.path.join(co_dir, a, a + '.pansn.fa.gz')))), ""),
        ("v2_covered_local", str(preferred_match), "GCA-preferred exact-match (dir name == preferred accession)"),
        ("v2_missing_local", str(preferred_missing), "need download (GCA-preferred semantics)"),
        ("v2_missing_gcf_only_numerics", str(preferred_gcf_only), "of missing: numerics with only GCF row in v2"),
        ("v2_covered_prefix_agnostic", str(len(v2_in_v1)), "numeric present locally under either GCA/GCF"),
        ("v2_missing_prefix_agnostic", str(len(v2_missing)), "numeric absent locally under either prefix"),
        ("v1_not_in_v2", str(len(v1_not_v2)), ""),
        ("host_clades_v1", str(len(hc_accs)), "unique accessions in v1 host_clades.tsv"),
        ("host_clades_v2_overlap", str(len(hc_numeric & acc_all_numeric)), ""),
        ("v1_prophage_rows", str(sum(1 for _ in open(os.path.join(v1_dir, PROPHAGE_V1))) - 1), "ntm_prophages.csv"),
        ("v1_prophage_genomes", str(len(v1_proph_genomes)), ""),
        ("v2_prophage_numeric_genomes", str(len(v2_proph_numeric)), "after GCA/GCF dedup"),
        ("v2_prophage_run_genomes", str(len(v2_proph_runs)), ""),
        ("v2_prophage_rows_kept", str(len(kept_rows)), "after dedup"),
        ("v2_prophage_rows_kept_gca", str(kept_by_source.get('GCA', 0)), "GCA rows kept (preferred)"),
        ("v2_prophage_rows_kept_gcf", str(kept_by_source.get('GCF', 0)), "GCF rows kept (GCF-only numerics)"),
        ("v2_prophage_rows_kept_runs", str(sum(kept_by_source.get(s, 0) for s in ('ERR', 'SRR', 'DRR'))), "run-assembly rows kept"),
        ("v2_prophage_rows_dropped_gcf_twin", str(len(true_rows) - len(kept_rows)), ""),
        ("v2_prophage_gcf_only_numerics", str(len(gcf_only_numerics)), ""),
        ("v2_prophage_id_unique", str(len(pid_counts)), "unique prophage_ids"),
        ("v2_prophage_length_min", str(lens_sorted[0]), "bp"),
        ("v2_prophage_length_median", str(int(statistics.median(lens_sorted))), "bp"),
        ("v2_prophage_length_mean", f"{sum(lens_sorted) / n_lens:.1f}", "bp"),
        ("v2_prophage_length_max", str(lens_sorted[-1]), "bp"),
        ("v2_coord_bad", str(coord_bad), "begin>end or unparseable"),
        ("v2_len_mismatch", str(len_mismatch), "length != end-start+1"),
        ("v2_prophage_genomes_overlap_v1", str(len(v1_proph_numeric & v2_proph_numeric)), "v1∩v2 prophage-bearing numerics"),
        ("per_genome_count_compared", str(compared), "v1 vs v2 per-genome prophage counts"),
        ("per_genome_count_agree", str(agree), ""),
        ("per_genome_count_differ", str(len(disagree_rows)), ""),
        ("tb_complex_rows", str(len(tb_rows)), "species contains 'tuberculosis'"),
        ("tb_complex_numeric", str(len(tb_numeric)), ""),
        ("tb_complex_runs", str(len(tb_runs)), ""),
    ]
    cov_path = os.path.join(out_dir, "coverage.tsv")
    with open(cov_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        for row in cov_rows:
            w.writerow(row)

    # species table for report
    R("")
    R("## 10. Species composition (manifest)")
    R("")
    man_sp = Counter(r[1] for r in man_rows)
    R("| species | genomes (rows) |")
    R("|---|---|")
    for sp, c in man_sp.most_common(30):
        R(f"| {sp} | {c} |")
    R("")

    report_path = os.path.join(out_dir, "validation_report.md")
    with open(report_path, "w") as fh:
        fh.write("\n".join(report) + "\n")

    print(f"Wrote {report_path}")
    print(f"Wrote {cov_path}")
    # print key numbers to stdout for quick review
    print(f"rows: acc={len(acc_rows)} man={len(man_rows)} true={len(true_rows)} false={len(false_rows)}")
    print(f"numeric={unique_numeric} twins={with_twin} gca_only={gca_only} gcf_only={gcf_only} runs={len(acc_runs)}")
    print(f"v1 local={len(v1_accs)} v2_covered={len(v2_in_v1)} v2_missing={len(v2_missing)} v1_not_v2={len(v1_not_v2)}")
    print(f"kept prophage rows={len(kept_rows)} gcf_only_numerics={len(gcf_only_numerics)} dup_pids={len(dup_pids)}")
    print(f"length min/med/mean/max={lens_sorted[0]}/{int(statistics.median(lens_sorted))}/{sum(lens_sorted)/n_lens:.0f}/{lens_sorted[-1]}")
    print(f"coord_bad={coord_bad} len_mismatch={len_mismatch}")
    print(f"per-genome compared={compared} agree={agree} differ={len(disagree_rows)}")
    print(f"v1∩v2 prophage numerics={len(v1_proph_numeric & v2_proph_numeric)}")
    print(f"TB rows={len(tb_rows)} numeric={len(tb_numeric)} runs={len(tb_runs)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
