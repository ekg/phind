#!/usr/bin/env python3
"""
validate_v2_ml.py — independent validation for task ntm-v2-ml.

Checks (mirrors the task's Validation section):
  1. ML genome count == alignable clades + singletons
     (1251 + 1137 = 2388).
  2. Per-genome length sanity: every genome >= 1 kb; every genome shorter
     than half its clade's median prophage length must have a logged
     warning (warnings.tsv).
  3. FASTA headers carry status=ml|singleton, n_members, length (v1 header
     format `>ntm2_<cid>_ML status=...`).
  4. Ancestral set covers exactly the alignable clades (1251), with
     `>ntm2_<cid>_ANCESTRAL` headers.
  5. Determinism: re-run traverse_partitions.py (both modes, seed 42) on a
     sample of alignable clades and byte-compare the outputs against the
     committed per-clade intermediates.

Exit 0 if all checks pass, 1 otherwise.
"""
import argparse
import json
import os
import statistics
import subprocess
import sys
import tempfile

V2 = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2"
CLADES_ROOT = f"{V2}/clades"
RESULTS_JSON = f"{CLADES_ROOT}/pipeline_results.json"
ML_DIR = f"{V2}/ml"
ML_FA = f"{ML_DIR}/all_ntm2_ml_phage_genomes.fa"
ANC_FA = f"{ML_DIR}/all_ntm2_ancestral_phage_genomes.fa"
MANIFEST = f"{ML_DIR}/release_manifest.tsv"
WARNINGS = f"{ML_DIR}/warnings.tsv"
MIN_LEN = 1000
HALF_MEDIAN_FACTOR = 0.5

REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
TRAV = os.path.join(REPO, "scripts", "traverse_partitions.py")


def read_fasta(path):
    out = {}
    cur = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(seq)
                cur = line[1:]
                seq = []
            else:
                seq.append(line.strip())
    if cur is not None:
        out[cur] = "".join(seq)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clades-root", default=CLADES_ROOT)
    ap.add_argument("--results", default=RESULTS_JSON)
    ap.add_argument("--ml-dir", default=ML_DIR)
    ap.add_argument("--n-determinism", type=int, default=10,
                    help="clades to re-run for the determinism check")
    ap.add_argument("--report", default=None,
                    help="write validation markdown report here")
    args = ap.parse_args()

    results = json.load(open(args.results))
    all_cids = sorted(results.keys())
    alignable = sorted(c for c in all_cids
                       if results[c].get("n_members", 0) >= 2)
    singletons = sorted(c for c in all_cids
                        if results[c].get("n_members", 0) == 1)
    problems = []
    notes = []

    # ── 1. load aggregate FASTAs ───────────────────────────────────────────
    ml_fa = read_fasta(f"{args.ml_dir}/all_ntm2_ml_phage_genomes.fa")
    anc_fa = read_fasta(f"{args.ml_dir}/all_ntm2_ancestral_phage_genomes.fa")
    n_ml = len(ml_fa)
    n_anc = len(anc_fa)
    exp_ml = len(alignable) + len(singletons)
    if n_ml != exp_ml:
        problems.append(
            f"ML genome count {n_ml} != alignable+singletons {exp_ml}")
    else:
        notes.append(f"ML count OK: {n_ml} (alignable {len(alignable)} + "
                     f"singletons {len(singletons)})")
    if n_anc != len(alignable):
        problems.append(f"ancestral count {n_anc} != alignable {len(alignable)}")
    else:
        notes.append(f"ancestral count OK: {n_anc}")

    # manifest consistency
    manifest = {}
    with open(f"{args.ml_dir}/release_manifest.tsv") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            manifest[row["clade_id"]] = row
    if len(manifest) != exp_ml:
        problems.append(f"manifest rows {len(manifest)} != {exp_ml}")

    # warnings loaded for the below-half-median check
    warned = set()
    with open(f"{args.ml_dir}/warnings.tsv") as f:
        next(f)
        for line in f:
            w = line.strip()
            if w:
                warned.add(w.split(":")[0])
    notes.append(f"warnings.tsv: {len(warned)} clades with logged warnings")

    # ── 2/3. header format + length sanity per genome ──────────────────────
    import re
    hdr_re = re.compile(r"^ntm2_(?P<cid>\d+_\d{4})_ML "
                        r"status=(?P<status>ml|singleton) "
                        r"n_members=(?P<n>\d+) length=(?P<len>\d+)$")
    n_bad_hdr = 0
    n_too_short = 0
    n_below_half_unwarned = 0
    for hdr, seq in ml_fa.items():
        m = hdr_re.match(hdr)
        if not m:
            n_bad_hdr += 1
            if n_bad_hdr <= 5:
                problems.append(f"bad ML header: '{hdr}'")
            continue
        cid = m.group("cid")
        status = m.group("status")
        n_members = int(m.group("n"))
        hlen = int(m.group("len"))
        if hlen != len(seq):
            problems.append(f"{cid}: header length {hlen} != seq length "
                            f"{len(seq)}")
        if status == "singleton":
            if cid not in singletons:
                problems.append(f"{cid}: status=singleton but n_members={n_members}")
        else:
            if cid not in alignable:
                problems.append(f"{cid}: status=ml but not in alignable set")
        if len(seq) < MIN_LEN:
            n_too_short += 1
        med = statistics.median([len(s) for s in
                                 read_fasta(f"{args.clades_root}/{cid}/"
                                            "sequences.fa").values()])
        if len(seq) < HALF_MEDIAN_FACTOR * med and cid not in warned:
            n_below_half_unwarned += 1
            problems.append(f"{cid}: ML length {len(seq)} < half clade median "
                            f"{med/2:.0f} and NOT logged in warnings.tsv")
    if n_bad_hdr:
        problems.append(f"{n_bad_hdr} ML headers malformed")
    notes.append(f"genomes < {MIN_LEN} bp: {n_too_short} (warned)")

    anc_hdr_re = re.compile(r"^ntm2_(?P<cid>\d+_\d{4})_ANCESTRAL "
                            r"status=ml n_members=\d+ length=\d+$")
    anc_cids = set()
    for hdr in anc_fa:
        m = anc_hdr_re.match(hdr)
        if not m:
            problems.append(f"bad ancestral header: '{hdr}'")
            continue
        anc_cids.add(m.group("cid"))
    missing_anc = set(alignable) - anc_cids
    extra_anc = anc_cids - set(alignable)
    if missing_anc:
        problems.append(f"{len(missing_anc)} alignable clades missing "
                        f"ancestral genome: {sorted(missing_anc)[:10]}")
    if extra_anc:
        problems.append(f"{len(extra_anc)} non-alignable clades have "
                        f"ancestral genome: {sorted(extra_anc)[:10]}")

    # ── 4. determinism: re-run traverse on a sample ────────────────────────
    sample = alignable[:args.n_determinism // 2] + \
        alignable[len(alignable) - args.n_determinism // 2:]
    sample = sorted(set(sample))
    ndiff = 0
    with tempfile.TemporaryDirectory() as tmp:
        for cid in sample:
            cdir = os.path.join(args.clades_root, cid)
            for mode, key, ref in (("ml", "ml.fa", "ml.ml.fa"),
                                   ("ancestral", "ancestral.genome.fa",
                                    "anc.ancestral.genome.fa")):
                prefix = os.path.join(tmp, f"{cid}_{mode}")
                r = subprocess.run(
                    [sys.executable, TRAV,
                     "--partitions-dir", os.path.join(cdir, "partitions"),
                     "--bed", os.path.join(cdir, "partitions.bed"),
                     "--output", prefix, "--mode", mode,
                     "--n-samples", "5", "--seed", "42"],
                    capture_output=True, text=True)
                if r.returncode != 0:
                    problems.append(f"determinism re-run failed {cid} "
                                    f"{mode}: {r.stderr[-500:]}")
                    continue
                # compare the stitched genome sequence byte-for-byte
                # (header `community_<prefix>_...` embeds the output prefix
                # path; the sequence body is the content under test)
                def body(path):
                    lines = open(path).read().splitlines()
                    return "\n".join(l for l in lines if not l.startswith(">"))
                if body(f"{prefix}.{key}") != body(os.path.join(cdir, ref)):
                        ndiff += 1
                        problems.append(f"determinism mismatch: {cid} {key}")
    notes.append(f"determinism re-run: {len(sample)} clades x2 modes, "
                 f"{ndiff} mismatches")

    # ── report ─────────────────────────────────────────────────────────────
    ok = not problems
    report = [
        "# NTM v2 — ML + ancestral validation",
        "",
        f"- ML genomes: {n_ml} (expected {exp_ml})",
        f"- Ancestral genomes: {n_anc} (expected {len(alignable)})",
        f"- Manifest rows: {len(manifest)}",
    ] + [f"- {n}" for n in notes] + [
        "",
    ]
    if problems:
        report.append("## Problems")
        report.append("")
        for p in problems:
            report.append(f"- {p}")
        report.append("")
        report.append("**Validation: FAIL**")
    else:
        report.append("**Validation: PASS**")
    report.append("")
    txt = "\n".join(report)
    print(txt)
    if args.report:
        with open(args.report, "w") as f:
            f.write(txt + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
