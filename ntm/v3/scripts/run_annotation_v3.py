#!/usr/bin/env python3
"""
run_annotation_v3.py — reproducible, restart-safe driver for the NTM v3
functional-QC pipeline (Pharokka + CheckV) on the v3 generated genomes:
1,304 ML + 893 ancestral = 2,197 records from
$NVME/ntm/v3/ml/ (all_ntm_v3_ml_phage_genomes.fa,
all_ntm_v3_ancestral_phage_genomes.fa).

Adapted from ntm/v2/scripts/run_annotation_v2.py with the same pinned tools,
databases and fast configuration:

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
ever READ (databases); every write goes to the separate v3 analysis root
(default /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/annotation). The driver
refuses to run if the v3 root is inside the v1 tree.

Layout (v3 analysis root):
  input/all_v3_phage_genomes.fa   2,197 records, header = exact genome ID only
  input/genome_index.tsv          genome_id, source (ntm3_ml|ntm3_anc), cohort,
                                  clade_id, length, status, host_clades, origin_fa
  pharokka_out/  checkv_out/  report/
  pharokka_out/per_genome_gff/    one GFF3 per genome + manifest.tsv
  run_state/<stage>.done.json     restart markers: command, versions, DB paths,
                                  thread count, timestamp, exit code
  run_state/tool_provenance.json  pinned versions + DB identity

Committed index (repo, one row per genome): genome_id -> gff path + gene count
  ntm/v3/annotation_gff_index.tsv

Stages (in order; each skipped if its marker exists and outputs are present):
  prepare  build combined input FASTA + index; FASTA/index round-trip check
           (exact, duplicate-free, 1,304 ML + 893 ancestral)
  pharokka annotate all records (single multi-FASTA run, meta mode)
  checkv   end_to_end on the same input FASTA
  split    one GFF3 per genome (from pharokka.gff, embedded ##FASTA omitted) +
           committed index TSV mapping genome_id -> gff path + gene count
  report   build_annotation_report_v3.py merge + graded classification

Usage:
  run_annotation_v3.py --dry-run          # prepare + validate + print commands
  run_annotation_v3.py                    # run everything (restart-safe)
  run_annotation_v3.py --stages pharokka  # run one stage
  run_annotation_v3.py --force --stages checkv   # redo a stage
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys
from collections import defaultdict

# --- pinned provenance -------------------------------------------------------
NVME = os.environ.get("NVME", "/mnt/nvme3n1/erikg/phind-genome-work")
RELEASE_DIR = os.path.join(NVME, "ntm", "v3", "ml")
ML_FA = os.path.join(RELEASE_DIR, "all_ntm_v3_ml_phage_genomes.fa")
ANC_FA = os.path.join(RELEASE_DIR, "all_ntm_v3_ancestral_phage_genomes.fa")
V1_ANNOTATION_DIR = os.path.join(NVME, "annotation")  # read-only
PHAROKKA_BIN = "/home/erikg/micromamba/envs/pharokka/bin/pharokka"
CHECKV_BIN = "/home/erikg/micromamba/envs/phage-annot/bin/checkv"
PHAROKKA_DB = os.path.join(V1_ANNOTATION_DIR, "pharokka_db")
CHECKV_DB = os.path.join(V1_ANNOTATION_DIR, "checkv_db", "checkv-db-v1.5")
DEFAULT_ROOT = os.path.join(NVME, "ntm", "v3", "annotation")
MAX_THREADS = 64
EXPECT_ML = 1304
EXPECT_ANC = 893
EXPECT_TOTAL = EXPECT_ML + EXPECT_ANC  # 2,197
ID_PREFIX = "ntm3_"

STAGES = ["prepare", "pharokka", "checkv", "split", "report"]

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT_BUILDER = os.path.join(HERE, "build_annotation_report_v3.py")
# committed per-genome GFF3 index (repo; GFF3 files live on NVMe under ROOT)
REPO_GFF_INDEX = os.path.normpath(os.path.join(HERE, "..", "annotation_gff_index.tsv"))


# --- small helpers -----------------------------------------------------------

def log(msg):
    print(f"[run_annotation_v3] {msg}", flush=True)


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


def clade_id_of(gid, src):
    if src == "ntm3_ml" and gid.endswith("_ML"):
        return gid[len(ID_PREFIX):-len("_ML")]
    if src == "ntm3_anc" and gid.endswith("_ANCESTRAL"):
        return gid[len(ID_PREFIX):-len("_ANCESTRAL")]
    raise RuntimeError(f"genome_id {gid!r} does not match source {src!r}")


def stage_prepare(root, force=False):
    """Build input/all_v3_phage_genomes.fa + input/genome_index.tsv.

    Genome IDs are preserved EXACTLY (first token of each release header).
    The combined FASTA is written with bare IDs so Pharokka/CheckV contig IDs
    round-trip 1:1 with the index. Duplicate IDs abort.
    """
    indir = os.path.join(root, "input")
    os.makedirs(indir, exist_ok=True)
    fa_out = os.path.join(indir, "all_v3_phage_genomes.fa")
    idx_out = os.path.join(indir, "genome_index.tsv")

    index_rows = []      # dicts in file order
    seen = {}
    dupes = []
    # single streaming pass: ML first, then ancestral; write bare-ID FASTA
    with open(fa_out, "w", newline="\n") as out:
        for src, path, status_default in (
            ("ntm3_ml", ML_FA, "ml"),
            ("ntm3_anc", ANC_FA, "ancestral"),
        ):
            if not os.path.exists(path):
                raise FileNotFoundError(path)
            for gid, desc, seq in iter_fasta(path):
                if gid in seen:
                    dupes.append(gid)
                    continue
                seen[gid] = src
                kv = parse_kv_desc(desc)
                status = "ancestral" if src == "ntm3_anc" else kv.get("status", status_default)
                index_rows.append({
                    "genome_id": gid,
                    "source": src,
                    "cohort": ("ntm3_anc" if src == "ntm3_anc" else
                               "ntm3_ml_singleton" if status == "singleton"
                               else "ntm3_ml_reconstructed"),
                    "clade_id": None,
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

    # validate every ID's suffix matches its source (strict, post-duplicate-check)
    for r in index_rows:
        r["clade_id"] = clade_id_of(r["genome_id"], r["source"])

    # ---- round-trip validation: FASTA <-> index, exact and duplicate-free ---
    n_ml = sum(1 for r in index_rows if r["source"] == "ntm3_ml")
    n_anc = sum(1 for r in index_rows if r["source"] == "ntm3_anc")
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
    for src, path in (("ntm3_ml", ML_FA), ("ntm3_anc", ANC_FA)):
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


# --- split stage -------------------------------------------------------------

GFF_INDEX_COLS = ["genome_id", "gff_path", "gene_count"]


def split_gff(root, index_out):
    """Split Pharokka's multi-contig pharokka.gff into one GFF3 per genome.

    Emits ROOT/pharokka_out/per_genome_gff/<genome_id>.gff (plus an on-NVMe
    manifest.tsv mirror) and the committed index (genome_id, gff_path,
    gene_count). Features are grouped by seqid (== bare genome ID). The
    embedded ##FASTA block is omitted (genomes live in the input FASTA). Any
    feature seqid not present in the prepared input FASTA aborts the stage.
    """
    src = os.path.join(root, "pharokka_out", "pharokka.gff")
    if not os.path.exists(src):
        raise FileNotFoundError(src)
    fa = os.path.join(root, "input", "all_v3_phage_genomes.fa")
    expected = [gid for gid, *_ in iter_fasta(fa)]
    expected_set = set(expected)

    regions = {}
    features = defaultdict(list)
    in_fasta = False
    with open(src) as f:
        for line in f:
            line = line.rstrip("\n")
            if in_fasta:
                continue
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    regions[parts[1]] = (parts[2], parts[3])
                continue
            if not line or line.startswith("#"):
                continue
            gid = line.split("\t", 1)[0]
            features[gid].append(line)

    unknown = sorted(set(features) - expected_set)
    if unknown:
        raise RuntimeError(f"pharokka.gff has {len(unknown)} seqid(s) not in the "
                           f"prepared input FASTA (first: {unknown[:3]})")

    outdir = os.path.join(root, "pharokka_out", "per_genome_gff")
    os.makedirs(outdir, exist_ok=True)
    index_rows = []
    n_features = 0
    for gid in expected:
        lines = features.get(gid, [])
        p = os.path.join(outdir, f"{gid}.gff")
        region = regions.get(gid)
        with open(p, "w", newline="\n") as out:
            out.write("##gff-version 3\n")
            if region:
                out.write(f"##sequence-region {gid} {region[0]} {region[1]}\n")
            for ln in lines:
                out.write(ln + "\n")
        index_rows.append((gid, p, len(lines)))
        n_features += len(lines)

    def _write_index(path):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w", newline="\n") as f:
            f.write("\t".join(GFF_INDEX_COLS) + "\n")
            for gid, p, n in index_rows:
                f.write(f"{gid}\t{p}\t{n}\n")

    manifest = os.path.join(outdir, "manifest.tsv")
    _write_index(manifest)
    _write_index(index_out)

    zero = [gid for gid, _p, n in index_rows if n == 0]
    log(f"split: wrote {len(index_rows)} per-genome GFF3 to {outdir} "
        f"({n_features} CDS features; {len(zero)} genomes with 0 features); "
        f"index -> {index_out}")
    return {"n_genomes": len(index_rows), "n_features": n_features,
            "n_zero_feature": len(zero), "gff_dir": outdir,
            "manifest": manifest, "index": index_out}


# --- stage commands ----------------------------------------------------------

def input_fa(root):
    return os.path.join(root, "input", "all_v3_phage_genomes.fa")


def pharokka_cmd(root, threads, tool_force=False):
    cmd = [PHAROKKA_BIN, "run", "-m", "--mmseqs2_only",
           "--skip_extra_annotations", "--skip_mash", "-g", "prodigal-gv",
           "-i", input_fa(root),
           "-o", os.path.join(root, "pharokka_out"),
           "-d", PHAROKKA_DB, "-t", str(threads), "--locustag", "NTMV3"]
    if tool_force:
        # restart safety: outdir exists without a done-marker = stale partial
        # run from a crashed attempt; pharokka refuses to overwrite without -f
        cmd.insert(cmd.index("run") + 1, "-f")
    return cmd


def checkv_cmd(root, threads):
    # checkv end_to_end resumes where it left off by default (--restart is
    # the documented default), so no force flag is needed for stale outdirs
    return [CHECKV_BIN, "end_to_end", input_fa(root),
            os.path.join(root, "checkv_out"),
            "-d", CHECKV_DB, "-t", str(threads)]


def report_cmd(root):
    return [sys.executable, REPORT_BUILDER, "--root", root]


STAGE_SPECS = {
    "prepare": {
        "cmd": lambda root, t, f=False: ["internal:stage_prepare"],
        "outputs": [os.path.join("{root}", "input", "all_v3_phage_genomes.fa"),
                    os.path.join("{root}", "input", "genome_index.tsv")],
    },
    "pharokka": {
        "cmd": lambda root, t, f=False: pharokka_cmd(root, t, f),
        "outputs": [os.path.join("{root}", "pharokka_out",
                                 "pharokka_cds_final_merged_output.tsv"),
                    os.path.join("{root}", "pharokka_out",
                                 "pharokka_length_gc_cds_density.tsv"),
                    os.path.join("{root}", "pharokka_out", "pharokka.gff")],
    },
    "checkv": {
        "cmd": lambda root, t, f=False: checkv_cmd(root, t),
        "outputs": [os.path.join("{root}", "checkv_out", "quality_summary.tsv")],
    },
    "split": {
        "cmd": lambda root, t, f=False: ["internal:split_gff"],
        "outputs": [os.path.join("{root}", "pharokka_out", "per_genome_gff",
                                 "manifest.tsv"),
                    REPO_GFF_INDEX],
    },
    "report": {
        "cmd": lambda root, t, f=False: report_cmd(root),
        "outputs": [os.path.join("{root}", "report",
                                 "per_genome_functional_qc.tsv")],
    },
}


def stage_env(bin_path):
    """subprocess env with the tool's own env bin prepended to PATH.

    Pharokka shells out to phanotate.py / other env binaries by name, so the
    pinned binary's env must be on PATH even when invoked by absolute path."""
    env = os.environ.copy()
    env["PATH"] = os.path.dirname(bin_path) + os.pathsep + env.get("PATH", "")
    return env


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
        if not os.path.exists(tmpl.format(root=root) if "{root}" in tmpl else tmpl):
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
        "committed_gff_index": REPO_GFF_INDEX,
        "report_builder": REPORT_BUILDER,
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
        raise RuntimeError(f"v3 analysis root {root} must be OUTSIDE the v1 "
                           f"evidence tree {v1}")
    for b in (PHAROKKA_BIN, CHECKV_BIN):
        if not os.path.exists(b):
            raise FileNotFoundError(f"pinned tool missing: {b}")
    for d in (PHAROKKA_DB, CHECKV_DB):
        if not os.path.isdir(d):
            raise FileNotFoundError(f"pinned database missing: {d}")
    if not os.path.exists(REPORT_BUILDER):
        raise FileNotFoundError(f"report builder missing: {REPORT_BUILDER}")


# --- main --------------------------------------------------------------------

def run_stage(root, stage, threads, force):
    if stage == "prepare":
        if stage_done(root, "prepare") and not force:
            log("prepare: marker present, re-validating index (cheap)")
        info = stage_prepare(root, force=force)
        write_marker(root, "prepare", ["internal:stage_prepare"], threads, 0)
        return info

    if stage == "split":
        if stage_done(root, "split") and not force:
            log("split: already complete (marker + outputs present) — skipping; "
                "use --force to redo")
            return {"skipped": True}
        info = split_gff(root, REPO_GFF_INDEX)
        write_marker(root, "split", ["internal:split_gff"], threads, 0)
        return info

    if stage_done(root, stage) and not force:
        log(f"{stage}: already complete (marker + outputs present) — skipping; "
            f"use --force to redo")
        return {"skipped": True}

    # restart safety: an existing outdir with no done-marker is a stale
    # partial run — pharokka needs -f to overwrite it
    tool_force = stage == "pharokka" and os.path.isdir(
        os.path.join(root, "pharokka_out"))
    cmd = STAGE_SPECS[stage]["cmd"](root, threads, tool_force)
    if tool_force:
        log(f"{stage}: stale output dir present without done-marker — "
            f"forcing overwrite")
    log(f"{stage}: running: {' '.join(cmd)}")
    env_bin = {"pharokka": PHAROKKA_BIN, "checkv": CHECKV_BIN,
               "report": sys.executable}.get(stage)
    env = stage_env(env_bin) if env_bin else None
    proc = subprocess.run(cmd, env=env)
    write_marker(root, stage, cmd, threads, proc.returncode)
    if proc.returncode != 0:
        raise RuntimeError(f"stage {stage} failed with exit code {proc.returncode}")
    return {"exit_code": proc.returncode}


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="NTM v3 functional-QC pipeline driver (Pharokka + CheckV)")
    ap.add_argument("--root", default=DEFAULT_ROOT,
                    help=f"v3 analysis root (default {DEFAULT_ROOT})")
    ap.add_argument("--threads", type=int, default=MAX_THREADS,
                    help=f"threads per tool invocation (cap {MAX_THREADS})")
    ap.add_argument("--stages", default=",".join(STAGES),
                    help=f"comma-separated subset of: {','.join(STAGES)}")
    ap.add_argument("--dry-run", action="store_true",
                    help="prepare + validate inputs, print pinned commands and "
                         "provenance, do NOT run pharokka/checkv/split/report")
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
        for st in ("pharokka", "checkv", "split", "report"):
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
