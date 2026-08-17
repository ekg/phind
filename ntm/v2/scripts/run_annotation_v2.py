#!/usr/bin/env python3
"""
run_annotation_v2.py — reproducible, restart-safe driver for the NTM v2
functional-QC pipeline (Pharokka + CheckV) on the current v2 generated
genomes: 2,388 ML + 1,251 ancestral = 3,639 records from
/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/release/.

Adapted from the established v1 E. coli/NTM workflow (research/phage_annotation)
with the same pinned tools, databases and fast configuration:

  pharokka v1.10.1  /home/erikg/micromamba/envs/pharokka/bin/pharokka
    pharokka run -m --mmseqs2_only --skip_extra_annotations --skip_mash
                 -g prodigal-gv -t <=64
  checkv  v1.1.1    /home/erikg/micromamba/envs/phage-annot/bin/checkv
    checkv end_to_end -t <=64
  PHAROKKA_DB  /mnt/nvme3n1/erikg/phind-genome-work/annotation/pharokka_db
               (PHROG DB, 9 Sep 2025 release — READ-ONLY shared copy)
  CHECKV_DB    /mnt/nvme3n1/erikg/phind-genome-work/annotation/checkv_db/checkv-db-v1.5
               (READ-ONLY shared copy)

The v1 evidence tree (/mnt/nvme3n1/erikg/phind-genome-work/annotation) is only
ever READ (databases); every write goes to the separate v2 analysis root
(default /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/annotation). The driver
refuses to run if the v2 root is inside the v1 tree.

Layout (v2 analysis root):
  input/all_v2_phage_genomes.fa   3,639 records, header = exact genome ID only
  input/genome_index.tsv          genome_id, source (ntm2_ml|ntm2_anc), cohort,
                                  clade_id, length, status, host_clades
  pharokka_out/  checkv_out/  report/
  run_state/<stage>.done.json     restart markers: command, versions, DB paths,
                                  thread count, timestamp, exit code
  run_state/tool_provenance.json  pinned versions + DB identity

Stages (in order; each skipped if its marker exists and outputs are present):
  prepare  build combined input FASTA + index; FASTA/index round-trip check
           (exact, duplicate-free, 2,388 ML + 1,251 ancestral)
  pharokka annotate all records (single multi-FASTA run, meta mode)
  checkv   end_to_end on the same input FASTA
  report   build_annotation_report_v2.py merge + graded classification

Usage:
  run_annotation_v2.py --dry-run          # prepare + validate + print commands
  run_annotation_v2.py                    # run everything (restart-safe)
  run_annotation_v2.py --stages pharokka  # run one stage
  run_annotation_v2.py --force --stages checkv   # redo a stage
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys

# --- pinned provenance -------------------------------------------------------
RELEASE_DIR = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/release"
ML_FA = os.path.join(RELEASE_DIR, "all_ntm_ml_phage_genomes.fa")
ANC_FA = os.path.join(RELEASE_DIR, "all_ntm_ancestral_phage_genomes.fa")
V1_ANNOTATION_DIR = "/mnt/nvme3n1/erikg/phind-genome-work/annotation"  # read-only
PHAROKKA_BIN = "/home/erikg/micromamba/envs/pharokka/bin/pharokka"
CHECKV_BIN = "/home/erikg/micromamba/envs/phage-annot/bin/checkv"
PHAROKKA_DB = os.path.join(V1_ANNOTATION_DIR, "pharokka_db")
CHECKV_DB = os.path.join(V1_ANNOTATION_DIR, "checkv_db", "checkv-db-v1.5")
DEFAULT_ROOT = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/annotation"
MAX_THREADS = 64
EXPECT_ML = 2388
EXPECT_ANC = 1251
EXPECT_TOTAL = EXPECT_ML + EXPECT_ANC  # 3,639

STAGES = ["prepare", "pharokka", "checkv", "report"]

HERE = os.path.dirname(os.path.abspath(__file__))


# --- small helpers -----------------------------------------------------------

def log(msg):
    print(f"[run_annotation_v2] {msg}", flush=True)


def sha256_file(path, _buf=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(_buf):
            h.update(chunk)
    return h.hexdigest()


def now_iso():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def tool_version(binpath, args):
    try:
        out = subprocess.run([binpath] + args, capture_output=True, text=True,
                             timeout=120)
        return (out.stdout + out.stderr).strip().splitlines()[0] if (out.stdout or out.stderr) else "unknown"
    except Exception as e:  # noqa: BLE001
        return f"unavailable: {e}"


def db_identity(db_dir):
    """Cheap, stable identity of a tool DB: path + per-file (name,size,mtime)."""
    info = {"path": db_dir, "files": []}
    if os.path.isdir(db_dir):
        for name in sorted(os.listdir(db_dir)):
            p = os.path.join(db_dir, name)
            if os.path.isfile(p):
                st = os.stat(p)
                info["files"].append({"name": name, "size": st.st_size,
                                      "mtime": int(st.st_mtime)})
    return info


# --- FASTA / index handling --------------------------------------------------

def iter_fasta(path):
    """Yield (header_first_token, description, [sequence_lines])."""
    with open(path) as f:
        gid = desc = None
        seq = []
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if gid is not None:
                    yield gid, desc, seq
                parts = line[1:].split(None, 1)
                gid = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                seq = []
            elif line:
                seq.append(line.strip())
        if gid is not None:
            yield gid, desc, seq


def parse_kv_desc(desc):
    out = {}
    for tok in desc.split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
    return out


def stage_prepare(root, force=False):
    """Build input/all_v2_phage_genomes.fa + input/genome_index.tsv.

    Genome IDs are preserved EXACTLY (first token of each release header).
    The combined FASTA is written with bare IDs so Pharokka/CheckV contig IDs
    round-trip 1:1 with the index. Duplicate IDs abort.
    """
    indir = os.path.join(root, "input")
    os.makedirs(indir, exist_ok=True)
    fa_out = os.path.join(indir, "all_v2_phage_genomes.fa")
    idx_out = os.path.join(indir, "genome_index.tsv")

    index_rows = []      # dicts in file order
    seen = {}
    dupes = []
    # single streaming pass: ML first, then ancestral; write bare-ID FASTA
    with open(fa_out, "w", newline="\n") as out:
        for src, path, status_default in (
            ("ntm2_ml", ML_FA, "ml"),
            ("ntm2_anc", ANC_FA, "ancestral"),
        ):
            if not os.path.exists(path):
                raise FileNotFoundError(path)
            for gid, desc, seq in iter_fasta(path):
                if gid in seen:
                    dupes.append(gid)
                    continue
                seen[gid] = src
                kv = parse_kv_desc(desc)
                status = "ancestral" if src == "ntm2_anc" else kv.get("status", status_default)
                index_rows.append({
                    "genome_id": gid,
                    "source": src,
                    "cohort": ("ntm2_anc" if src == "ntm2_anc" else
                               "ntm2_ml_singleton" if status == "singleton"
                               else "ntm2_ml_reconstructed"),
                    "clade_id": (gid[len("ntm2_"):-len("_ML")] if src == "ntm2_ml"
                                 else gid[len("ntm2_"):-len("_ANCESTRAL")]),
                    "length": sum(len(s) for s in seq),
                    "status": status,
                    "host_clades": kv.get("host_clades", ""),
                    "origin_fa": os.path.basename(path),
                })
                # write record with the EXACT genome ID only (no description)
                out.write(f">{gid}\n")
                for s in seq:
                    out.write(s + "\n")

    if dupes:
        raise RuntimeError(f"duplicate genome IDs across release FASTAs: {dupes[:5]} "
                           f"({len(dupes)} total)")

    if dupes:
        raise RuntimeError(f"duplicate genome IDs across release FASTAs: {dupes[:5]} "
                           f"({len(dupes)} total)")

    # ---- round-trip validation: FASTA <-> index, exact and duplicate-free ---
    n_ml = sum(1 for r in index_rows if r["source"] == "ntm2_ml")
    n_anc = sum(1 for r in index_rows if r["source"] == "ntm2_anc")
    fa_ids = [gid for gid, *_ in iter_fasta(fa_out)]
    fa_set = set(fa_ids)
    idx_ids = [r["genome_id"] for r in index_rows]
    problems = []
    if len(fa_set) != len(fa_ids):
        problems.append(f"combined FASTA has duplicate IDs "
                        f"({len(fa_ids)} records, {len(fa_set)} unique)")
    if sorted(fa_ids) != sorted(idx_ids):
        problems.append("FASTA IDs do not round-trip to index IDs")
    if n_ml != EXPECT_ML:
        problems.append(f"expected {EXPECT_ML} ML records, got {n_ml}")
    if n_anc != EXPECT_ANC:
        problems.append(f"expected {EXPECT_ANC} ancestral records, got {n_anc}")
    if len(idx_ids) != len(set(idx_ids)):
        problems.append("index has duplicate genome_ids")
    if problems:
        raise RuntimeError("prepare validation FAILED: " + "; ".join(problems))

    # cross-check against source release headers (exact ID preservation)
    for src, path in (("ntm2_ml", ML_FA), ("ntm2_anc", ANC_FA)):
        src_ids = {gid for gid, *_ in iter_fasta(path)}
        want = {r["genome_id"] for r in index_rows if r["source"] == src}
        if src_ids != want:
            raise RuntimeError(f"ID set mismatch vs {path}")

    # header length self-check: index length == FASTA-computed length
    fa_len = {gid: sum(len(s) for s in seq) for gid, _d, seq in iter_fasta(fa_out)}
    for r in index_rows:
        if fa_len.get(r["genome_id"]) != r["length"]:
            raise RuntimeError(f"length mismatch for {r['genome_id']}")

    cols = ["genome_id", "source", "cohort", "clade_id", "length", "status",
            "host_clades", "origin_fa"]
    with open(idx_out, "w", newline="\n") as f:
        f.write("\t".join(cols) + "\n")
        for r in index_rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    log(f"prepare: wrote {fa_out} ({len(index_rows)} records, "
        f"{n_ml} ML + {n_anc} ancestral; round-trip exact, duplicate-free)")
    return {"n_total": len(index_rows), "n_ml": n_ml, "n_anc": n_anc,
            "input_fa": fa_out, "input_fa_sha256": sha256_file(fa_out),
            "genome_index": idx_out}


# --- stage commands ----------------------------------------------------------

def pharokka_cmd(root, threads):
    return [PHAROKKA_BIN, "run", "-m", "--mmseqs2_only",
            "--skip_extra_annotations", "--skip_mash", "-g", "prodigal-gv",
            "-i", os.path.join(root, "input", "all_v2_phage_genomes.fa"),
            "-o", os.path.join(root, "pharokka_out"),
            "-d", PHAROKKA_DB, "-t", str(threads), "--locustag", "NTMV2"]


def checkv_cmd(root, threads):
    return [CHECKV_BIN, "end_to_end",
            os.path.join(root, "input", "all_v2_phage_genomes.fa"),
            os.path.join(root, "checkv_out"),
            "-d", CHECKV_DB, "-t", str(threads)]


def report_cmd(root):
    return [sys.executable,
            os.path.join(HERE, "build_annotation_report_v2.py"),
            "--root", root]


STAGE_SPECS = {
    "pharokka": {
        "cmd": lambda root, t: pharokka_cmd(root, t),
        "outputs": [os.path.join("{root}", "pharokka_out",
                                 "pharokka_cds_final_merged_output.tsv"),
                    os.path.join("{root}", "pharokka_out",
                                 "pharokka_length_gc_cds_density.tsv")],
    },
    "checkv": {
        "cmd": lambda root, t: checkv_cmd(root, t),
        "outputs": [os.path.join("{root}", "checkv_out", "quality_summary.tsv")],
    },
    "report": {
        "cmd": lambda root, t: report_cmd(root),
        "outputs": [os.path.join("{root}", "report",
                                 "per_genome_functional_qc.tsv")],
    },
}


# --- restart safety ----------------------------------------------------------

def marker_path(root, stage):
    return os.path.join(root, "run_state", f"{stage}.done.json")


def stage_done(root, stage):
    m = marker_path(root, stage)
    if not os.path.exists(m):
        return False
    try:
        with open(m) as f:
            rec = json.load(f)
    except json.JSONDecodeError:
        return False
    if rec.get("exit_code") != 0:
        return False
    for tmpl in STAGE_SPECS[stage]["outputs"]:
        if not os.path.exists(tmpl.format(root=root)):
            return False
    return True


def write_marker(root, stage, cmd, threads, exit_code):
    os.makedirs(os.path.join(root, "run_state"), exist_ok=True)
    rec = {"stage": stage, "command": " ".join(cmd), "threads": threads,
           "exit_code": exit_code, "finished_utc": now_iso(),
           "pharokka_version": tool_version(PHAROKKA_BIN, ["version"]),
           "checkv_version": tool_version(CHECKV_BIN, ["--version"]),
           "pharokka_db": db_identity(PHAROKKA_DB),
           "checkv_db": db_identity(CHECKV_DB)}
    with open(marker_path(root, stage), "w", newline="\n") as f:
        json.dump(rec, f, indent=2, sort_keys=True)
        f.write("\n")
    return rec


def write_provenance(root, threads):
    os.makedirs(os.path.join(root, "run_state"), exist_ok=True)
    prov = {
        "generated_utc": now_iso(),
        "release_dir": RELEASE_DIR,
        "release_files": {
            os.path.basename(p): {"sha256": sha256_file(p),
                                  "size": os.path.getsize(p)}
            for p in (ML_FA, ANC_FA) if os.path.exists(p)},
        "expect": {"ml": EXPECT_ML, "ancestral": EXPECT_ANC,
                   "total": EXPECT_TOTAL},
        "pharokka_bin": PHAROKKA_BIN,
        "pharokka_version": tool_version(PHAROKKA_BIN, ["version"]),
        "pharokka_db": db_identity(PHAROKKA_DB),
        "checkv_bin": CHECKV_BIN,
        "checkv_version": tool_version(CHECKV_BIN, ["--version"]),
        "checkv_db": db_identity(CHECKV_DB),
        "max_threads_cap": MAX_THREADS,
        "threads": threads,
        "v1_annotation_dir_read_only": V1_ANNOTATION_DIR,
    }
    p = os.path.join(root, "run_state", "tool_provenance.json")
    with open(p, "w", newline="\n") as f:
        json.dump(prov, f, indent=2, sort_keys=True)
        f.write("\n")
    return prov


# --- safety ------------------------------------------------------------------

def guard_paths(root):
    root = os.path.realpath(root)
    v1 = os.path.realpath(V1_ANNOTATION_DIR)
    if root == v1 or root.startswith(v1 + os.sep):
        raise RuntimeError(f"v2 analysis root {root} must be OUTSIDE the v1 "
                           f"evidence tree {v1}")
    for b in (PHAROKKA_BIN, CHECKV_BIN):
        if not os.path.exists(b):
            raise FileNotFoundError(f"pinned tool missing: {b}")
    for d in (PHAROKKA_DB, CHECKV_DB):
        if not os.path.isdir(d):
            raise FileNotFoundError(f"pinned database missing: {d}")


# --- main --------------------------------------------------------------------

def run_stage(root, stage, threads, force):
    if stage == "prepare":
        if stage_done(root, "prepare") and not force:
            log("prepare: marker present, re-validating index (cheap)")
        info = stage_prepare(root, force=force)
        write_marker(root, "prepare", ["internal:stage_prepare"], threads, 0)
        return info

    if stage_done(root, stage) and not force:
        log(f"{stage}: already complete (marker + outputs present) — skipping; "
            f"use --force to redo")
        return {"skipped": True}

    cmd = STAGE_SPECS[stage]["cmd"](root, threads)
    log(f"{stage}: running: {' '.join(cmd)}")
    proc = subprocess.run(cmd)
    write_marker(root, stage, cmd, threads, proc.returncode)
    if proc.returncode != 0:
        raise RuntimeError(f"stage {stage} failed with exit code {proc.returncode}")
    return {"exit_code": proc.returncode}


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="NTM v2 functional-QC pipeline driver (Pharokka + CheckV)")
    ap.add_argument("--root", default=DEFAULT_ROOT,
                    help=f"v2 analysis root (default {DEFAULT_ROOT})")
    ap.add_argument("--threads", type=int, default=MAX_THREADS,
                    help=f"threads per tool invocation (cap {MAX_THREADS})")
    ap.add_argument("--stages", default=",".join(STAGES),
                    help=f"comma-separated subset of: {','.join(STAGES)}")
    ap.add_argument("--dry-run", action="store_true",
                    help="prepare + validate inputs, print pinned commands and "
                         "provenance, do NOT run pharokka/checkv/report")
    ap.add_argument("--force", action="store_true",
                    help="redo stages even if markers exist")
    args = ap.parse_args(argv)

    if args.threads > MAX_THREADS:
        raise SystemExit(f"--threads {args.threads} exceeds cap {MAX_THREADS}")
    if args.threads < 1:
        raise SystemExit("--threads must be >= 1")

    guard_paths(args.root)
    os.makedirs(args.root, exist_ok=True)
    prov = write_provenance(args.root, args.threads)
    log(f"provenance: pharokka={prov['pharokka_version']} "
        f"checkv={prov['checkv_version']} threads={args.threads} "
        f"(cap {MAX_THREADS})")

    stages = [s.strip() for s in args.stages.split(",") if s.strip()]
    bad = [s for s in stages if s not in STAGES]
    if bad:
        raise SystemExit(f"unknown stages: {bad}")

    if args.dry_run:
        if "prepare" not in stages:
            stages.insert(0, "prepare")
        info = stage_prepare(args.root, force=args.force)
        write_marker(args.root, "prepare", ["internal:stage_prepare (dry-run)"],
                     args.threads, 0)
        log("DRY RUN — would execute (restart-safe, markers in run_state/):")
        for st in ("pharokka", "checkv", "report"):
            if st in stages:
                log(f"  [ {' '.join(STAGE_SPECS[st]['cmd'](args.root, args.threads))} ]")
        log(f"inputs: {info['n_total']} unique records "
            f"({info['n_ml']} ML + {info['n_anc']} ancestral) — round-trip "
            f"exact, duplicate-free; input FASTA sha256={info['input_fa_sha256'][:16]}…")
        log(f"provenance written: {args.root}/run_state/tool_provenance.json")
        return 0

    for st in stages:
        run_stage(args.root, st, args.threads, args.force)
    log("all requested stages complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
