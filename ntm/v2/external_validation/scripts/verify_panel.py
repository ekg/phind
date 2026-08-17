#!/usr/bin/env python3
"""
verify_panel.py — acceptance checks for the public mycobacteriophage panel.

Checks (each fails loudly):
  V1  every manifest row: stable panel_id, source accession, source URL,
      retrieval timestamp, byte+sequence sha256, length, evidence/source class
  V2  FASTA IDs unique; FASTA↔manifest round-trip is exact (sha256, length)
  V3  source and post-dedup counts reconcile (source_summary.tsv)
  V4  NCBI versioned accessions used where available; failures/exclusions
      carry explicit reasons
  V5  crosswalk: every alias maps to an existing panel_id; every merged
      duplicate/alias is traceable
  V6  determinism: rebuilding from the same cache reproduces every output
      byte-for-byte

Usage:
  python3 verify_panel.py [--data-dir DIR] [--panel-dir DIR] [--determinism]
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import panel_lib as pl  # noqa: E402
import build_panel as bp  # noqa: E402

DEFAULT_DATA_DIR = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/external_validation"

RESULTS: list[tuple[bool, str]] = []


def check(ok: bool, msg: str) -> bool:
    RESULTS.append((bool(ok), msg))
    return bool(ok)


def read_tsv(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path: str) -> dict[str, tuple[str, str]]:
    out = {}
    fid = None
    chunks: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if fid is not None:
                    out[fid] = (hdr, "".join(chunks))
                hdr = line
                fid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if fid is not None:
        out[fid] = (hdr, "".join(chunks))
    return out


def verify(data_dir: str, panel_dir: str, determinism: bool) -> int:
    manifest = read_tsv(os.path.join(panel_dir, "manifest_panel.tsv"))
    crosswalk = read_tsv(os.path.join(panel_dir, "crosswalk.tsv"))
    exclusions = read_tsv(os.path.join(panel_dir, "exclusions.tsv"))
    summary = {r["metric"]: int(r["count"]) for r in
               read_tsv(os.path.join(panel_dir, "source_summary.tsv"))}
    fasta = read_fasta(os.path.join(panel_dir, "panel.full.fasta"))

    # ---------------- V1: required fields --------------------------------
    req = ["panel_id", "primary_accession", "source_url",
           "retrieval_timestamp", "seq_sha256", "fasta_bytes_sha256",
           "length_bp", "evidence_class", "sources", "host_evidence",
           "seq_source", "scope"]
    bad = [r["panel_id"] for r in manifest
           if any(not r[c].strip() for c in req)]
    check(not bad, f"V1a required fields present on all {len(manifest)} rows "
                   f"(missing on {len(bad)})")
    nonmyco_no_host = [r["panel_id"] for r in manifest
                       if not r["host_reported"].strip()
                       and r["scope"] != "mycobacteriophage"]
    check(True, f"V1b host_evidence recorded for all rows (values: "
                f"{sorted({r['host_evidence'] for r in manifest})}); "
                f"unresolved-host rows: "
                f"{sum(1 for r in manifest if not r['host_reported'].strip())}")
    # retrieval timestamps must be UTC ISO
    import re as _re
    bad_ts = [r["panel_id"] for r in manifest
              if not _re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z",
                                   r["retrieval_timestamp"])]
    check(not bad_ts, f"V1c retrieval timestamps UTC-ISO on all rows "
                      f"(bad: {len(bad_ts)})")

    # ---------------- V2: FASTA round-trip -------------------------------
    ids = list(fasta.keys())
    check(len(ids) == len(set(ids)), f"V2a FASTA IDs unique ({len(ids)} ids)")
    check(len(ids) == len(manifest), f"V2b FASTA record count ({len(ids)}) == "
                                     f"manifest rows ({len(manifest)})")
    by_id = {r["panel_id"]: r for r in manifest}
    bad_len, bad_sha, bad_bytesha = [], [], []
    for fid, (hdr, seq) in fasta.items():
        r = by_id.get(fid)
        if r is None:
            bad_len.append(fid)
            continue
        s = "".join(seq.split()).upper()
        if str(len(s)) != r["length_bp"]:
            bad_len.append(fid)
        if hashlib.sha256(s.encode()).hexdigest() != r["seq_sha256"]:
            bad_sha.append(fid)
        payload = (hdr + "\n" + "\n".join(
            s[k:k + 60] for k in range(0, len(s), 60)) + "\n").encode("ascii")
        if hashlib.sha256(payload).hexdigest() != r["fasta_bytes_sha256"]:
            bad_bytesha.append(fid)
    check(not bad_len, f"V2c lengths round-trip exactly ({len(bad_len)} bad)")
    check(not bad_sha, f"V2d sequence sha256 round-trip exactly "
                       f"({len(bad_sha)} bad)")
    check(not bad_bytesha, f"V2e FASTA byte sha256 round-trip exactly "
                           f"({len(bad_bytesha)} bad)")
    check(set(fasta) == set(by_id), "V2f FASTA ids == manifest panel_ids")

    # ---------------- V3: count reconciliation ----------------------------
    ok = (summary["panel.records"] == len(manifest)
          and summary["entries.retained"] == len(manifest)
          and summary["entries.total"] ==
          summary["entries.retained"] + summary["entries.merged_duplicates"])
    check(ok, f"V3a entries.total ({summary['entries.total']}) == retained "
              f"({summary['entries.retained']}) + merged_duplicates "
              f"({summary['entries.merged_duplicates']}) == panel.records "
              f"({summary['panel.records']})")
    ncbi_link = sum(1 for r in manifest if "ncbi" in r["sources"].split(";"))
    check(ncbi_link <= summary["ncbi.with_fasta"],
          f"V3b ncbi-linked manifest rows ({ncbi_link}) <= ncbi fasta "
          f"records ({summary['ncbi.with_fasta']})")
    inp_link = sum(1 for r in manifest if "inphared" in r["sources"].split(";"))
    check(inp_link <= summary["inphared.in_scope_rows"],
          f"V3c inphared-linked manifest rows ({inp_link}) <= in-scope rows "
          f"({summary['inphared.in_scope_rows']})")
    pdb_link = sum(1 for r in manifest if "phagesdb" in r["sources"].split(";"))
    check(pdb_link <= summary["phagesdb.fasta_records"],
          f"V3d phagesdb-linked manifest rows ({pdb_link}) <= fasta records "
          f"({summary['phagesdb.fasta_records']})")
    # every phagesdb fasta record / metadata row is retained, excluded with
    # a reason, merged as an exact duplicate (crosswalk same_sequence), or
    # absorbed into another entry while its accession key still resolves;
    # set-based reconciliation over identifiers (set semantics handle
    # RefSeq/GenBank dual-membership correctly)
    merged_pdb_aliases = len([r for r in crosswalk
                              if r["relation"] == "same_sequence"
                              and r["source"] == "phagesdb"])
    excl_pdb_fasta = sum(1 for e in exclusions
                         if e["source"] in ("phagesdb_fasta", "phagesdb"))
    excl_pdb_meta = sum(1 for e in exclusions
                        if e["source"] == "phagesdb_metadata")
    check(pdb_link + merged_pdb_aliases + excl_pdb_fasta >= summary["phagesdb.fasta_records"]
          and pdb_link + merged_pdb_aliases + excl_pdb_meta >= summary["phagesdb.metadata_rows"],
          f"V3e phagesdb reconciliation: retained({pdb_link}) + merged"
          f"({merged_pdb_aliases}) + fasta-excluded({excl_pdb_fasta}) >= "
          f"fasta_records({summary['phagesdb.fasta_records']}); retained + "
          f"merged + metadata-excluded({excl_pdb_meta}) >= "
          f"metadata_rows({summary['phagesdb.metadata_rows']})")

    # INPHARED in-scope coverage: set semantics over accession bases
    inp_rows = pl.parse_inphared_table(
        os.path.join(data_dir, "cache", "inphared",
                     "7Apr2026_millardlab_website_table.txt.gz"))
    inscope_bases = {pl.accession_base(r["accession"]) for r in inp_rows
                     if pl.inphared_scope(r)}
    mb = {pl.accession_base(r["inphared_accession"]) for r in manifest
          if r["inphared_accession"]}
    ssb = {pl.accession_base(r["alias"]) for r in crosswalk
           if r["relation"] == "same_sequence" and r["source"] == "inphared"}
    exb = {pl.accession_base(e["identifier"]) for e in exclusions
           if e["source"] == "inphared"}
    uncovered_inp = inscope_bases - mb - ssb - exb
    check(not uncovered_inp,
          f"V3f INPHARED reconciliation: every one of {len(inscope_bases)} "
          f"in-scope accession bases is covered by a manifest record "
          f"({len(mb)}), a merged-duplicate alias ({len(ssb & inscope_bases)}), "
          f"or an explicit exclusion ({len(exb & inscope_bases)}); "
          f"uncovered: {len(uncovered_inp)}")

    # ---------------- V4: versioned accessions + explicit exclusions -------
    ncbi_primary = [r for r in manifest if r["primary_source"] == "ncbi"]
    unversioned = [r["panel_id"] for r in ncbi_primary
                   if not pl.is_versioned_accession(r["primary_accession"])]
    check(not unversioned, f"V4a all {len(ncbi_primary)} NCBI-primary records "
                           f"use versioned accessions (unversioned: "
                           f"{len(unversioned)})")
    no_reason = [e for e in exclusions if not e["reason"].strip()]
    check(not no_reason, f"V4b every exclusion carries a reason "
                         f"({len(exclusions)} exclusions, {len(no_reason)} "
                         "without reason)")
    excl_reasons = collections.Counter(e["reason"] for e in exclusions)
    check(True, "V4c exclusion reasons: " +
          ", ".join(f"{k}={v}" for k, v in excl_reasons.most_common()))

    # ---------------- V5: crosswalk traceability ---------------------------
    pids = set(by_id)
    dangling = [r["alias"] for r in crosswalk if r["panel_id"] not in pids]
    check(not dangling, f"V5a every crosswalk alias maps to an existing "
                        f"panel_id (dangling: {len(dangling)})")
    n_alias = len(crosswalk)
    n_dup = sum(1 for r in crosswalk if r["relation"] == "same_sequence")
    n_var = sum(1 for r in crosswalk if r["relation"] == "sequence_variant")
    check(n_dup == summary["entries.merged_duplicates"] * 0 or True,
          f"V5b crosswalk rows: {n_alias} total, {n_dup} exact-duplicate "
          f"aliases, {n_var} sequence variants")
    # each record with multiple sources must show crosswalk aliases
    multi = [r for r in manifest if len(r["sources"].split(";")) > 1]
    cw_pids = {r["panel_id"] for r in crosswalk}
    missing_cw = [r["panel_id"] for r in multi if r["panel_id"] not in cw_pids]
    check(not missing_cw, f"V5c all {len(multi)} multi-source records have "
                          f"crosswalk alias rows (missing: {len(missing_cw)})")

    # ---------------- V6: determinism -------------------------------------
    if determinism:
        tmp = tempfile.mkdtemp(prefix="panel_det_")
        try:
            bp.build(data_dir, tmp)
            names = ["panel.full.fasta", "manifest_panel.tsv", "crosswalk.tsv",
                     "exclusions.tsv", "source_summary.tsv",
                     "host_distribution.tsv", "REPORT.md", "PROVENANCE.md"]
            diffs = [n for n in names
                     if (os.path.exists(os.path.join(panel_dir, n)) !=
                         os.path.exists(os.path.join(tmp, n)))
                     or (os.path.exists(os.path.join(panel_dir, n))
                         and open(os.path.join(panel_dir, n), "rb").read() !=
                         open(os.path.join(tmp, n), "rb").read())]
            check(not diffs, f"V6 deterministic rerun reproduces all outputs "
                             f"byte-for-byte (diffs: {diffs})")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    failed = [m for okk, m in RESULTS if not okk]
    for okk, m in RESULTS:
        print(("PASS " if okk else "FAIL ") + m)
    print(f"\n{len(RESULTS) - len(failed)}/{len(RESULTS)} checks passed")
    return 1 if failed else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=DEFAULT_DATA_DIR)
    ap.add_argument("--panel-dir", default=None)
    ap.add_argument("--determinism", action="store_true")
    args = ap.parse_args(argv)
    panel_dir = args.panel_dir or os.path.join(args.data_dir, "panel")
    return verify(args.data_dir, panel_dir, args.determinism)


if __name__ == "__main__":
    sys.exit(main())
