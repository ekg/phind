#!/usr/bin/env python3
"""
NTM v3 — build the unified all-inclusive prophage manifest (BV-BRC phigaro
primary, v2 collaborator calls for genomes absent from the BV-BRC export).

Convention (v2 extract_report.md, carried into v3):
  prophage coordinates come from the source manifests — we do NOT re-run the
  caller. v3 reconciles the two callers per genome:
    * genomes present in the BV-BRC phigaro export  -> phigaro coordinates
      (source=BV-BRC, caller=phigaro), regardless of export row namespace
      (ASSEMBLY run / NCBI assembly / BV-BRC taxid.version);
    * genomes absent from the export (v2-only NCBI assemblies and the 135
      v2-only runs) -> v2 collaborator coordinates (source=V2, caller=v2-collab).
  A genome NEVER carries calls from both callers: on overlap the v2 set is
  recorded as superseded (counted here, not emitted).

Export-internal dedup (measured 2026-09-13, see extract_report.md):
  1. stub rows — the coordinates CSV carries 22,225 rows with an empty
     prophage_id (float coords, empty transposable). 5,927 of them duplicate a
     named row of the same genome (same scaffold/begin/end) in 2,549 genomes
     whose summary prophage_count double-counts; the other 16,298 are the only
     representation of their prophage and keep a synthesized prophage_id
     `<scaffold>_prophage<N>` (next free integer per scaffold, coordinate order).
  2. twin assemblies — 1,870 (numeric, version) assemblies are held as two
     canonical objects (GCA_ and GCF_ twins of the same sequences, both
     linked from their own export rows); and 364 further objects receive
     calls from more than one export genome (BV-BRC taxid.version variants,
     GCA/GCF rows sharing one object). Exactly one export genome is kept per
     canonical assembly, GCA-preferred (v2 convention): the row whose
     canonical object is the GCA twin first, then the NCBI row whose
     accession IS the canonical object, then any NCBI row (GCA preferred),
     then the BV-BRC row with the most prophages (highest taxid.version
     tie-break). Version-distinct objects (16 numerics with >1 version in
     the export) stay separate — they are distinct assemblies the
     acquisition downloaded at exact versions.

v2 side dedup replicates ntm/v2/scripts/extract_full_ntm_prophages.py exactly
(GCA preferred per numeric.version key; 23,243 kept rows). Supersession is
numeric-level for assemblies (GCA/GCF twins and version skew collapse to one
genome) and accession-level for runs.

Inputs (repo):
  ntm/v3/inputs/ntm_qc_passed_phigaro_summary_20260909.csv
  ntm/v3/inputs/ntm_qc_passed_phigaro_coordinates_20260909.csv
  ntm/v3/inputs/v3_acquisition_manifest.tsv.gz
Inputs (NVMe, frozen):
  {work}/ntm/v2/inputs/NTM_QC_passed_prophage_master_manifest.tsv
Outputs:
  repo  ntm/v3/inputs/v3_prophage_manifest.tsv.gz   (gzip -n, deterministic)
  NVMe  {work}/ntm/v3/inputs/v3_prophage_manifest.tsv  (plain copy)
        {work}/ntm/v3/scratch/v3_prophage_reconciliation.json (stats)

Usage:
  python3 ntm/v3/scripts/build_v3_prophage_manifest.py \
      [--work-dir /mnt/nvme3n1/erikg/phind-genome-work]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SUMMARY_CSV = REPO / "ntm/v3/inputs/ntm_qc_passed_phigaro_summary_20260909.csv"
COORDS_CSV = REPO / "ntm/v3/inputs/ntm_qc_passed_phigaro_coordinates_20260909.csv"
ACQ_TSVGZ = REPO / "ntm/v3/inputs/v3_acquisition_manifest.tsv.gz"

GC_RE = re.compile(r"^GC([AF])_(\d{9})\.(\d+)$")
RUN_RE = re.compile(r"^(ERR|SRR|DRR)\d+$")
V2_NUMERIC_RE = re.compile(r"^[A-Z]{3}_(\d+\.\d+)$")   # v2 script's exact key
PID_N_RE = re.compile(r"_prophage(\d+)$")

FIELDS = ["genome_id", "source", "accession", "species", "scaffold",
          "prophage_id", "begin", "end", "length", "transposable", "caller",
          "taxonomy", "canonical_acc", "genome_status", "notes"]


def numeric_of(acc: str) -> str | None:
    m = GC_RE.match(acc)
    return m.group(2) if m else None


def v2_preferred(rows: list[dict]) -> list[dict]:
    """v2 GCA/GCF dedup, replicating the v2 extract script exactly.

    Key = numeric.version of the accession (v2's NUMERIC_RE group(1)); run
    accessions fall back to the accession itself.
    """
    by_key: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        m = V2_NUMERIC_RE.match(r["accession"])
        by_key[m.group(1) if m else r["accession"]].append(r)
    kept: list[dict] = []
    for _key, rl in by_key.items():
        accs = {r["accession"] for r in rl}
        if len(accs) == 1:
            pa = next(iter(accs))
        else:
            pa = next((a for a in accs if a.startswith("GCA_")), None)
            if pa is None:
                pa = next((a for a in accs if a.startswith("GCF_")), None)
            if pa is None:
                pa = next(iter(accs))
        kept.extend(r for r in rl if r["accession"] == pa)
    return kept


def icoord(v: str) -> int:
    return int(float(v))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", default="/mnt/nvme3n1/erikg/phind-genome-work")
    args = ap.parse_args()
    work = Path(args.work_dir)
    v2_manifest = work / "ntm/v2/inputs/NTM_QC_passed_prophage_master_manifest.tsv"
    out_gz = REPO / "ntm/v3/inputs/v3_prophage_manifest.tsv.gz"
    out_nvme = work / "ntm/v3/inputs/v3_prophage_manifest.tsv"
    stats_path = work / "ntm/v3/scratch/v3_prophage_reconciliation.json"
    out_nvme.parent.mkdir(parents=True, exist_ok=True)
    stats_path.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------- load
    summary: dict[str, dict] = {}
    with open(SUMMARY_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            summary[r["genome_id"]] = r
    assert len(summary) == 26499, len(summary)

    acq: dict[str, dict] = {}
    with gzip.open(ACQ_TSVGZ, "rt", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            acq[r["genome_key"]] = r
    assert len(acq) == 34846, len(acq)

    coords: list[dict] = []
    with open(COORDS_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            coords.append(r)
    assert len(coords) == 47994, len(coords)

    # ------------------------------------------- stub dedup (export rows)
    # unique key (genome_id, scaffold, begin, end); a named row wins over a
    # stub row with the same key.
    by_uniq: dict[tuple, dict] = {}
    n_stub_rows = n_stub_dup = n_stub_kept = 0
    for r in coords:
        is_stub = not r["prophage_id"].strip()
        if is_stub:
            n_stub_rows += 1
        key = (r["genome_id"], r["scaffold"], icoord(r["begin"]), icoord(r["end"]))
        cur = by_uniq.get(key)
        if cur is None:
            by_uniq[key] = r
            if is_stub:
                n_stub_kept += 1
        elif is_stub and cur["prophage_id"].strip():
            n_stub_dup += 1          # stub duplicate of a named row: drop stub
        elif not is_stub and not cur["prophage_id"].strip():
            by_uniq[key] = r          # named row replaces stub: un-count stub
            n_stub_kept -= 1
            n_stub_dup += 1
    unique_rows = list(by_uniq.values())
    assert len(unique_rows) == 42067, len(unique_rows)
    assert n_stub_rows == n_stub_dup + n_stub_kept, (n_stub_rows, n_stub_dup, n_stub_kept)

    # synthesize prophage_ids for kept stubs: next free integer per
    # (genome, scaffold), assigned in begin order
    rows_by_gid: dict[str, list[dict]] = defaultdict(list)
    for r in unique_rows:
        rows_by_gid[r["genome_id"]].append(r)
    for gid, rs in rows_by_gid.items():
        max_n: dict[str, int] = defaultdict(int)
        for r in rs:
            m = PID_N_RE.search(r["prophage_id"]) if r["prophage_id"].strip() else None
            if m:
                max_n[r["scaffold"]] = max(max_n[r["scaffold"]], int(m.group(1)))
        stubs = sorted((r for r in rs if not r["prophage_id"].strip()),
                       key=lambda r: (icoord(r["begin"]), icoord(r["end"])))
        for r in stubs:
            max_n[r["scaffold"]] += 1
            r["prophage_id"] = f"{r['scaffold']}_prophage{max_n[r['scaffold']]}"

    # --------------------------- export-internal assembly dedup
    # genome key: run accession for ASSEMBLY rows; (numeric, version) for
    # rows resolved to a GC object (GCA/GCF twins of one assembly collapse,
    # v2 convention; version skew stays distinct); other canonical_acc
    # (BV-BRC downloads) as-is.
    def export_key(gid: str) -> tuple:
        a = acq[gid]
        if a["source_ns"] == "ASSEMBLY":
            return ("run", gid)
        m = GC_RE.match(a["canonical_acc"])
        if m:
            return ("asm", m.group(2), m.group(3))
        return ("obj", a["canonical_acc"] or gid)

    def export_rank(gid: str):
        a, s = acq[gid], summary[gid]
        n = len(rows_by_gid[gid])
        can_gca = 0 if a["canonical_acc"].startswith("GCA_") else 1
        if s["source"] == "NCBI":
            exact = 0 if a["accession"] == a["canonical_acc"] else 1
            pref = 0 if a["accession"].startswith("GCA_") else 1
            return (0, can_gca, exact, pref, a["accession"], -n)
        if s["source"] == "BV-BRC":
            return (1, can_gca, 1, 1, gid, -n)
        return (2, 1, 1, 1, gid, -n)

    groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for gid in rows_by_gid:
        groups[export_key(gid)].append(gid)
    kept_export_gid: dict[str, list[str]] = {}
    export_internal_superseded_rows = 0
    export_internal_superseded_genomes = 0
    for key, gids in groups.items():
        pick = sorted(gids, key=export_rank)[0]
        kept_export_gid[pick] = gids
        export_internal_superseded_rows += sum(
            len(rows_by_gid[g]) for g in gids if g != pick)
        export_internal_superseded_genomes += len(gids) - 1
    multi_groups = sum(1 for g in groups.values() if len(g) > 1)
    twin_asm_groups = sum(
        1 for k, g in groups.items()
        if k[0] == "asm" and len({acq[x]["canonical_acc"] for x in g}) > 1)

    # export presence sets (whole export, not just prophage-bearing genomes)
    export_numerics: set[str] = set()
    export_runs: set[str] = set()
    for gid, s in summary.items():
        a = acq[gid]
        if a["source_ns"] == "ASSEMBLY":
            export_runs.add(gid)
        else:
            m = GC_RE.match(a["accession"]) or GC_RE.match(a["canonical_acc"])
            if m:
                export_numerics.add(m.group(2))

    # ------------------------------------------------------- export rows
    out_rows: list[dict] = []
    for pick in sorted(kept_export_gid, key=lambda g: (export_key(g)[0], export_key(g)[1])):
        a, s = acq[pick], summary[pick]
        status = "blocked_run" if a["plan"] == "blocked_run" else "extractable"
        can = a["canonical_acc"] if status != "blocked_run" else a["accession"]
        notes = ""
        if len(kept_export_gid[pick]) > 1:
            sup = sorted(g for g in kept_export_gid[pick] if g != pick)
            notes = (f"supersedes export twin(s): {','.join(sup)}"
                     if len(sup) <= 4 else
                     f"supersedes {len(sup)} export twin rows")
        for r in sorted(rows_by_gid[pick],
                        key=lambda r: (r["scaffold"], icoord(r["begin"]),
                                       icoord(r["end"]), r["prophage_id"])):
            begin, end = icoord(r["begin"]), icoord(r["end"])
            out_rows.append({
                "genome_id": pick,
                "source": "BV-BRC",
                "accession": s["accession"],
                "species": s["species"],
                "scaffold": r["scaffold"],
                "prophage_id": r["prophage_id"],
                "begin": str(begin),
                "end": str(end),
                "length": str(end - begin + 1),
                "transposable": r["transposable"],
                "caller": "phigaro",
                "taxonomy": r["taxonomy"],
                "canonical_acc": can,
                "genome_status": status,
                "notes": notes,
            })

    n_export_rows = len(out_rows)

    # ----------------------------------------------------------- v2 rows
    with open(v2_manifest, newline="") as fh:
        v2_all = list(csv.DictReader(fh, delimiter="\t"))
    v2_true = [r for r in v2_all if r["has_prophage"] == "True"]
    v2_kept = v2_preferred(v2_true)
    assert len(v2_kept) == 23243, len(v2_kept)

    # acquisition lookup for v2-only genomes: numeric -> V2_NUMERIC row,
    # run -> V2_RUN row
    acq_v2num: dict[str, dict] = {}
    acq_v2run: dict[str, dict] = {}
    for a in acq.values():
        if a["origin"] == "v2_only" and a["source_ns"] == "V2_NUMERIC":
            acq_v2num[a["numeric"]] = a
        elif a["origin"] == "v2_only" and a["source_ns"] == "V2_RUN":
            acq_v2run[a["accession"]] = a

    n_v2_sup_num_rows = n_v2_sup_run_rows = n_v2_kept_rows = 0
    v2_sup_numerics: set[str] = set()
    v2_sup_runs = 0
    for r in v2_kept:
        acc = r["accession"]
        if RUN_RE.match(acc):
            if acc in export_runs:
                n_v2_sup_run_rows += 1
                v2_sup_runs += 1
                continue
            a = acq_v2run[acc]
            can, status = acc, "blocked_run"
        else:
            m = V2_NUMERIC_RE.match(acc)
            assert m, acc
            num = m.group(1).split(".")[0]
            if num in export_numerics:
                n_v2_sup_num_rows += 1
                v2_sup_numerics.add(num)
                continue
            a = acq_v2num[num]
            can, status = a["canonical_acc"], "extractable"
        n_v2_kept_rows += 1
        begin, end = icoord(r["prophage_start"]), icoord(r["prophage_end"])
        out_rows.append({
            "genome_id": r["genome_id"],
            "source": "V2",
            "accession": acc,
            "species": r["species"],
            "scaffold": r["prophage_contig"],
            "prophage_id": r["prophage_id"],
            "begin": str(begin),
            "end": str(end),
            "length": str(end - begin + 1),
            "transposable": r["prophage_transposable_element"],
            "caller": "v2-collab",
            "taxonomy": r["prophage_family"],
            "canonical_acc": can,
            "genome_status": status,
            "notes": "v2 genome absent from BV-BRC export",
        })
    assert n_v2_kept_rows == 619, n_v2_kept_rows

    # ------------------------------------------------------ sanity gates
    keys = [(r["canonical_acc"], r["prophage_id"]) for r in out_rows]
    dup = [k for k, c in Counter(keys).items() if c > 1]
    assert not dup, f"duplicate (canonical_acc, prophage_id): {dup[:5]}"
    genomes = {(r["canonical_acc"], r["source"]) for r in out_rows}
    per_genome_src: dict[str, set[str]] = defaultdict(set)
    for _, src in genomes:
        pass
    for r in out_rows:
        per_genome_src[r["canonical_acc"]].add(r["source"])
    mixed = [g for g, s in per_genome_src.items() if len(s) > 1]
    assert not mixed, f"genome with two caller sources: {mixed[:5]}"
    for r in out_rows:
        assert int(r["length"]) == int(r["end"]) - int(r["begin"]) + 1, r

    # ---------------------------------------------------------- outputs
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=FIELDS, delimiter="\t",
                       lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
    w.writeheader()
    for r in out_rows:
        w.writerow(r)
    payload = buf.getvalue().encode()
    with open(out_gz, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9, mtime=0) as fh:
            fh.write(payload)
    with open(out_nvme, "wb") as fh:
        fh.write(payload)

    extractable = sum(1 for r in out_rows if r["genome_status"] == "extractable")
    blocked = len(out_rows) - extractable
    stats = {
        "coordinates_csv_rows": len(coords),
        "stub_rows": n_stub_rows,
        "stub_rows_dropped_as_duplicates": n_stub_dup,
        "stub_rows_kept_synthesized_pid": n_stub_kept,
        "unique_prophages_after_stub_dedup": len(unique_rows),
        "genomes_double_counted_in_summary": 2549,
        "export_internal_multi_genome_groups": multi_groups,
        "export_internal_twin_assembly_groups": twin_asm_groups,
        "export_internal_superseded_genomes": export_internal_superseded_genomes,
        "export_internal_superseded_rows": export_internal_superseded_rows,
        "export_genomes_kept": len(kept_export_gid),
        "export_manifest_rows": n_export_rows,
        "export_rows_extractable": sum(
            1 for r in out_rows
            if r["source"] == "BV-BRC" and r["genome_status"] == "extractable"),
        "export_rows_blocked_run": sum(
            1 for r in out_rows
            if r["source"] == "BV-BRC" and r["genome_status"] == "blocked_run"),
        "v2_manifest_rows": len(v2_all),
        "v2_has_prophage_rows": len(v2_true),
        "v2_kept_after_gca_dedup": len(v2_kept),
        "v2_superseded_numeric_rows": n_v2_sup_num_rows,
        "v2_superseded_numeric_genomes": len(v2_sup_numerics),
        "v2_superseded_run_rows": n_v2_sup_run_rows,
        "v2_superseded_run_genomes": v2_sup_runs,
        "v2_kept_manifest_rows": n_v2_kept_rows,
        "v2_kept_numeric_rows": sum(
            1 for r in out_rows
            if r["source"] == "V2" and not RUN_RE.match(r["accession"])),
        "v2_kept_run_rows": sum(
            1 for r in out_rows
            if r["source"] == "V2" and RUN_RE.match(r["accession"])),
        "unified_manifest_rows": len(out_rows),
        "unified_extractable_rows": extractable,
        "unified_blocked_run_rows": blocked,
    }
    stats_path.write_text(json.dumps(stats, indent=1) + "\n")
    print(json.dumps(stats, indent=1))
    print(f"wrote {out_gz} ({len(out_rows)} rows) + plain copy {out_nvme}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
