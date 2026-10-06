#!/usr/bin/env python3
"""
NTM v3 — derive the delivery-updated acquisition manifest from the frozen
acquisition manifest (build_v3_acquisition_manifest.py output).

The 2026-09 collaborator delivery landed the previously-pending run assemblies
locally as `{work}/ntm/v3/genomes/delivery_raw/sra_assembled/<run>.fasta`
(17,920 files). This script flips the acquisition manifest's `blocked_run`
rows for which the delivered FASTA now exists:

  plan              blocked_run -> delivered_run
  resolution_method *           -> collaborator_delivery
  canonical_acc     ""          -> accession

Every other row and every other column is left byte-identical. The frozen
repo inputs are never modified in place: the output is a new manifest.

A blocked_run row with no delivered FASTA (the 135 V2_RUN rows) is left
untouched (still plan=blocked_run, still an empty canonical_acc) — the
collaborator delivery covers only the export ASSEMBLY runs, so those runs
stay genuinely blocked. The script fails loudly if any non-blocked_run row
would be flipped, or if a delivered FASTA matches no blocked_run accession.

Output (repo, committed, gzip -n deterministic):
  ntm/v3/inputs/v3_acquisition_manifest_delivery.tsv.gz

Usage:
  python3 ntm/v3/scripts/build_v3_acquisition_manifest_delivery.py \
      [--work-dir /mnt/nvme3n1/erikg/phind-genome-work] \
      [--acquisition ntm/v3/inputs/v3_acquisition_manifest.tsv.gz] \
      [--raw-dir {work}/ntm/v3/genomes/delivery_raw/sra_assembled] \
      [--out-manifest ntm/v3/inputs/v3_acquisition_manifest_delivery.tsv.gz]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
ACQ_TSVGZ = REPO / "ntm/v3/inputs/v3_acquisition_manifest.tsv.gz"

FIELDS = ["genome_key", "origin", "source_ns", "accession", "numeric",
          "species", "prophage_count", "canonical_acc", "plan",
          "resolution_method", "link_gz", "bvbrc_contigs", "notes"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", default="/mnt/nvme3n1/erikg/phind-genome-work")
    ap.add_argument("--acquisition", default=str(ACQ_TSVGZ),
                    help="frozen acquisition manifest (gzipped TSV) to read")
    ap.add_argument("--raw-dir", default=None,
                    help="delivered run FASTA dir "
                         "(default {work}/ntm/v3/genomes/delivery_raw/sra_assembled)")
    ap.add_argument("--out-manifest", default=str(
        REPO / "ntm/v3/inputs/v3_acquisition_manifest_delivery.tsv.gz"),
        help="delivery-updated acquisition manifest output path (gzip -n)")
    args = ap.parse_args()
    work = Path(args.work_dir)
    acquisition_path = Path(args.acquisition)
    raw_dir = Path(args.raw_dir) if args.raw_dir else \
        work / "ntm/v3/genomes/delivery_raw/sra_assembled"
    out_manifest = Path(args.out_manifest)
    out_manifest.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------- load
    with gzip.open(acquisition_path, "rt", newline="") as fh:
        raw = fh.read()
    rows = list(csv.DictReader(io.StringIO(raw), delimiter="\t"))
    assert list(rows[0].keys()) == FIELDS, list(rows[0].keys())
    assert len(rows) == 34846, len(rows)

    delivered = {p.name[:-len(".fasta")]
                 for p in raw_dir.iterdir() if p.name.endswith(".fasta")}
    print(f"delivered run FASTAs in {raw_dir}: {len(delivered)}")

    # a delivered FASTA must only ever match a blocked_run row (captured
    # before any flip, since the flip itself clears blocked_run)
    pre_non_blocked = {r["accession"] for r in rows if r["plan"] != "blocked_run"}
    assert not (delivered & pre_non_blocked), \
        f"delivered FASTA matching a non-blocked row: {sorted(delivered & pre_non_blocked)[:5]}"

    # ------------------------------------------------------- flip rows
    n_flip = 0
    flipped_ns: Counter = Counter()
    for r in rows:
        if r["plan"] != "blocked_run":
            continue
        if r["accession"] not in delivered:
            continue
        r["plan"] = "delivered_run"
        r["resolution_method"] = "collaborator_delivery"
        r["canonical_acc"] = r["accession"]
        n_flip += 1
        flipped_ns[r["source_ns"]] += 1

    flipped_accs = {r["accession"] for r in rows if r["plan"] == "delivered_run"}
    assert flipped_accs <= delivered, "flipped accession without delivered FASTA"
    leftover = delivered - flipped_accs
    assert not leftover, f"delivered FASTA matching no blocked_run row: {sorted(leftover)[:5]}"

    # ------------------------------------------ not-found (still blocked)
    still_blocked = [r for r in rows if r["plan"] == "blocked_run"]
    still_blocked_ns: Counter = Counter(r["source_ns"] for r in still_blocked)

    # ---------------------------------------------------------- output
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=FIELDS, delimiter="\t",
                       lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
    w.writeheader()
    for r in rows:
        w.writerow(r)
    with open(out_manifest, "wb") as raw_out:
        with gzip.GzipFile(fileobj=raw_out, mode="wb", compresslevel=9,
                           mtime=0) as fh:
            fh.write(buf.getvalue().encode("utf-8"))

    print(f"flipped blocked_run -> delivered_run: {n_flip} "
          f"({dict(flipped_ns)})")
    print(f"still blocked_run: {len(still_blocked)} ({dict(still_blocked_ns)})")
    print(f"wrote {out_manifest} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
