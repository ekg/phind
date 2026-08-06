#!/usr/bin/env python3
"""
QC for the ML phage genome release (research/ml_phage_genomes/).

For every clade in the release FASTA:
  * sequence length (already in manifest, re-verified from the FASTA)
  * GC content of the ML genome
  * longest run of ambiguous Ns
  * completeness: fraction of clade members whose sequence is
    covered >= 50% by alignment to the ML genome (mappy/minimap2)
  * identity distribution: per-member identity of the best
    ML-genome-to-member alignment (median / min / max over members)

Emits:
  research/ml_phage_genomes/qc_table.tsv          per-clade QC rows
  research/ml_phage_genomes/qc_flagged.tsv        clades that fail QC

Flags (implausible):
  length < 10 kb                 too_short
  length > 200 kb                too_long
  max_N_run >= 100               long_n_run
  completeness < 0.5             low_completeness
  median_identity < 0.8          low_identity
  member_count == 1 and status==singleton are NOT flagged (expected)

Usage:
  python3 scripts/qc_ml_phage_genomes.py \
      --fasta research/ml_phage_genomes/all_ml_phage_genomes.fa \
      --manifest research/ml_phage_genomes/release_manifest.tsv \
      --clades-root /home/erikg/phind/research/clades \
      --out research/ml_phage_genomes/qc_table.tsv \
      [--threads 8]
"""
import argparse
import csv
import json
import os
import re
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor

import mappy


MIN_LEN = 10_000
MAX_LEN = 200_000
N_RUN_FLAG = 100
COMPLETENESS_THRESHOLD = 0.5
COMPLETE_MEMBER_COV = 0.5  # member considered aligned if >= this fraction covered
IDENTITY_FLAG = 0.80


def load_fasta(path):
    """Return {header_id: sequence}."""
    out = {}
    cur = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(seq)
                cur = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if cur is not None:
        out[cur] = "".join(seq)
    return out


def load_members(clades_root, community, clade_id):
    """Return {prophage_element_id: sequence} for a clade."""
    path = os.path.join(clades_root, community, clade_id, "sequences.fa")
    if not os.path.exists(path):
        return None
    return load_fasta(path)


def gc_and_nrun(seq):
    gc = (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0
    runs = re.findall(r"N+", seq.upper())
    max_run = max((len(r) for r in runs), default=0)
    return gc, max_run


def map_stats(query, members):
    """Map ML genome to members; return (complete_frac, identities list)."""
    if not members:
        return 0.0, []
    names = list(members.keys())
    seqs = list(members.values())
    # build index over members via a temp FASTA (mappy file-based index)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
        for name, mseq in zip(names, seqs):
            tf.write(f">{name}\n{mseq}\n")
        tmp_path = tf.name
    try:
        a = mappy.Aligner(fn_idx_in=tmp_path, preset="asm5", best_n=1, n_threads=1)
    finally:
        os.unlink(tmp_path)
    if a is None:
        return 0.0, []
    member_covered = 0
    identities = []
    for name, mseq in zip(names, seqs):
        best = None  # best (fraction_covered, identity)
        for hit in a.map(query):
            if hit.ctg == name:
                qlen = hit.q_en - hit.q_st
                frac = qlen / len(mseq) if len(mseq) else 0.0
                if frac < COMPLETE_MEMBER_COV:
                    continue
                ident = 1.0 - hit.mlen / qlen if qlen else 0.0
                if best is None or frac > best[0]:
                    best = (frac, ident)
        if best is not None:
            member_covered += 1
            identities.append(best[1])
    complete_frac = member_covered / len(members)
    return complete_frac, identities


def fasta_id_to_clade(fasta_id):
    """clade_0_0000_ML -> 0_0000"""
    fid = fasta_id
    if fid.startswith("clade_"):
        fid = fid[len("clade_"):]
    if fid.endswith("_ML"):
        fid = fid[: -len("_ML")]
    return fid


def qc_clade(args):
    fasta, manifest_row, clades_root, fasta_id = args
    seq = fasta.get(fasta_id)
    if seq is None:
        return None
    clade_id = fasta_id_to_clade(fasta_id)
    length = len(seq)
    gc, max_n_run = gc_and_nrun(seq)
    community = manifest_row.get("community", clade_id.split("_")[0])
    members = load_members(clades_root, community, clade_id)
    complete_frac, identities = map_stats(seq, members) if members else (None, [])
    ident_med = sorted(identities)[len(identities) // 2] if identities else None
    ident_min = min(identities) if identities else None
    ident_max = max(identities) if identities else None
    n_members = len(members) if members is not None else int(manifest_row.get("n_members", 0))
    status = manifest_row.get("status", "ml")
    return {
        "clade_id": clade_id,
        "community": community,
        "n_members": n_members,
        "status": status,
        "length_bp": length,
        "gc_content": round(gc, 4),
        "max_n_run": max_n_run,
        "completeness": round(complete_frac, 4) if complete_frac is not None else "",
        "identity_median": round(ident_med, 4) if ident_med is not None else "",
        "identity_min": round(ident_min, 4) if ident_min is not None else "",
        "identity_max": round(ident_max, 4) if ident_max is not None else "",
        "n_members_aligned": len(identities),
    }


def flag_row(row):
    flags = []
    try:
        length = int(row["length_bp"])
        if length < MIN_LEN:
            flags.append("too_short")
        if length > MAX_LEN:
            flags.append("too_long")
    except (ValueError, TypeError):
        pass
    try:
        if int(row["max_n_run"]) >= N_RUN_FLAG:
            flags.append("long_n_run")
    except (ValueError, TypeError):
        pass
    try:
        if float(row["completeness"]) < COMPLETENESS_THRESHOLD:
            flags.append("low_completeness")
    except (ValueError, TypeError):
        pass
    try:
        if float(row["identity_median"]) < IDENTITY_FLAG:
            flags.append("low_identity")
    except (ValueError, TypeError):
        pass
    return ",".join(flags)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--clades-root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    fasta = load_fasta(args.fasta)
    manifest = {}
    with open(args.manifest) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            manifest[row["clade_id"]] = row

    clade_ids = [cid for cid in fasta.keys()]
    tasks = [
        (fasta, manifest.get(fasta_id_to_clade(cid), {}), args.clades_root, cid)
        for cid in clade_ids
    ]
    rows = []
    with ProcessPoolExecutor(max_workers=args.threads) as ex:
        for res in ex.map(qc_clade, tasks):
            if res is not None:
                rows.append(res)

    # flag
    for r in rows:
        r["flags"] = flag_row(r)

    # write QC table
    cols = [
        "clade_id", "community", "n_members", "status", "length_bp",
        "gc_content", "max_n_run", "completeness", "n_members_aligned",
        "identity_median", "identity_min", "identity_max", "flags",
    ]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["clade_id"]):
            w.writerow({k: r.get(k, "") for k in cols})

    # flagged clades
    flagged = [r for r in rows if r.get("flags")]
    outdir = os.path.dirname(args.out)
    flag_path = os.path.join(outdir, "qc_flagged.tsv")
    with open(flag_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in sorted(flagged, key=lambda x: x["clade_id"]):
            w.writerow({k: r.get(k, "") for k in cols})

    n_flagged = len(flagged)
    n_short = sum(1 for r in flagged if "too_short" in r["flags"])
    n_long = sum(1 for r in flagged if "too_long" in r["flags"])
    n_nrun = sum(1 for r in flagged if "long_n_run" in r["flags"])
    n_comp = sum(1 for r in flagged if "low_completeness" in r["flags"])
    n_id = sum(1 for r in flagged if "low_identity" in r["flags"])
    print(f"QC done: {len(rows)} clades, {n_flagged} flagged "
          f"(too_short={n_short} too_long={n_long} n_run={n_nrun} "
          f"low_completeness={n_comp} low_identity={n_id})")
    print(f"wrote {args.out} and {flag_path}")


if __name__ == "__main__":
    sys.exit(main())
