#!/usr/bin/env python3
"""
build_v3_ml_genomes.py — NTM v3 ML + ancestral traversal over every prophage
clade (task ntm-v3-ml).

Runs ``scripts/traverse_partitions.py`` in both modes (``--mode ml`` and
``--mode ancestral``, ``--n-samples 5 --seed 42`` — v2 convention) over every
alignable (partitioned) clade under the NVMe v3 clades root, and passes
singletons through as-is (v1/v2 convention: first record of
``sequences.fa``, ``status=singleton``).

Per-clade intermediates (v2 convention) are written into the clade dir:
  <clade>/ml.ml.fa, ml.consensus.fa, ml.traversal.json, ml.coverage.tsv,
  ml.stats.json  (--mode ml)
  <clade>/anc.ancestral.genome.fa, anc.ancestral.fa, anc.trees.nwk, ...
  (--mode ancestral)

Resume: clades whose complete per-mode output sets already exist (non-empty)
are skipped and their existing genomes are reused — worker/provider deaths
never lose compute (same pattern as the ntm-v3-per partition driver).

Aggregates are written to <out> (default NVMe ntm/v3/ml):
  all_ntm_v3_ml_phage_genomes.fa         one ML genome per clade
                                         (alignable + singletons = 813)
  all_ntm_v3_ancestral_phage_genomes.fa  one ancestral genome per alignable
                                         clade (456)
  release_manifest.tsv                   per-clade metadata: member count,
                                         host species set, source prophages,
                                         genome files, sanity flags
  per_clade_stats.tsv                    per-clade traversal stats
  warnings.tsv                           length-sanity warnings
  ml_report.md                           run report

Header format (v2 convention, ntm3 prefix):
  >ntm3_<clade_id>_ML status=ml|singleton n_members=<N> length=<L>
  >ntm3_<clade_id>_ANCESTRAL status=ml n_members=<N> length=<L>

Usage:
  python3 ntm/v3/scripts/build_v3_ml_genomes.py \
      [--clades-root /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/clades] \
      [--clades ntm/v3/clades/tight_clades.json.gz] \
      [--prophage-manifest ntm/v3/inputs/v3_prophage_manifest.tsv.gz] \
      [--out /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/ml] \
      [--jobs 16] [--n-samples 5] [--seed 42] [--dry-run] [--no-resume]
"""
import argparse
import concurrent.futures
import gzip
import json
import os
import statistics
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
TRAV = os.path.join(REPO, "scripts", "traverse_partitions.py")

V3 = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3"
CLADES_ROOT = f"{V3}/clades"
TIGHT_CLADES = os.path.join(REPO, "ntm", "v3", "clades", "tight_clades.json.gz")
PROPHAGE_MANIFEST = os.path.join(REPO, "ntm", "v3", "inputs",
                                 "v3_prophage_manifest.tsv.gz")
OUT_DIR = f"{V3}/ml"
MIN_LEN = 1000            # validation: no genome shorter than 1 kb
HALF_MEDIAN_FACTOR = 0.5  # warn if genome < factor * clade median member len

ML_OUTPUTS = [".traversal.json", ".consensus.fa", ".ml.fa", ".coverage.tsv",
              ".stats.json"]
ANC_OUTPUTS = [".ancestral.fa", ".ancestral.genome.fa", ".trees.nwk",
               ".coverage.tsv", ".stats.json"]


def load_clades(path):
    """Return {clade_id: [member prophage ids]} from tight_clades json(.gz)."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        return json.load(f)


def load_species(path):
    """Return {accession: host species} from the v3 prophage manifest."""
    op = gzip.open if path.endswith(".gz") else open
    acc_species = {}
    with op(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        ai, si = header.index("accession"), header.index("species")
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > si:
                acc_species[cols[ai]] = cols[si]
    return acc_species


def member_acc(member_id):
    """GCA_002801175.1#1#NQTL01000003.1_prophage1 -> GCA_002801175.1"""
    return member_id.split("#")[0]


def first_record(fasta):
    """Return (header, sequence) of the first record of a FASTA file."""
    name = None
    seq = []
    with open(fasta) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0]
            else:
                seq.append(line.strip())
    if name is None:
        raise RuntimeError(f"no records in {fasta}")
    return name, "".join(seq)


def read_all_fasta(fasta):
    """Return {name: sequence}."""
    out = {}
    cur = None
    seq = []
    with open(fasta) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(seq)
                cur = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if cur is not None:
        out[cur] = "".join(seq)
    return out


def member_lengths(cdir):
    """Per-clade prophage lengths (bp) from sequences.fa."""
    seqs = read_all_fasta(os.path.join(cdir, "sequences.fa"))
    return [len(s) for s in seqs.values()]


def outputs_complete(prefix, suffixes):
    """True if every expected output file exists and is non-empty."""
    for suf in suffixes:
        p = prefix + suf
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            return False
    return True


def run_traverse(cid, cdir, mode, n_samples, seed, python):
    """Run traverse_partitions.py for one clade/mode; return output prefix."""
    prefix = os.path.join(cdir, "ml" if mode == "ml" else "anc")
    cmd = [python, TRAV,
           "--partitions-dir", os.path.join(cdir, "partitions"),
           "--bed", os.path.join(cdir, "partitions.bed"),
           "--output", prefix,
           "--mode", mode,
           "--n-samples", str(n_samples),
           "--seed", str(seed)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(
            f"traverse {mode} failed for {cid}: {r.stderr[-2000:]} "
            f"({r.stdout[-2000:]})")
    return prefix


def process_clade(args):
    """Worker: run both traverse modes for one alignable clade (with resume).

    Returns (cid, result_dict) with ML + ancestral genome sequences,
    per-mode stats and how many modes were recomputed (0, 1 or 2)."""
    cid, cdir, n_samples, seed, python, resume = args
    t0 = time.time()
    recomputed = 0
    ml_prefix = os.path.join(cdir, "ml")
    anc_prefix = os.path.join(cdir, "anc")
    if not (resume and outputs_complete(ml_prefix, ML_OUTPUTS)):
        run_traverse(cid, cdir, "ml", n_samples, seed, python)
        recomputed += 1
    if not (resume and outputs_complete(anc_prefix, ANC_OUTPUTS)):
        run_traverse(cid, cdir, "ancestral", n_samples, seed, python)
        recomputed += 1
    _, ml_seq = first_record(ml_prefix + ".ml.fa")
    _, anc_seq = first_record(anc_prefix + ".ancestral.genome.fa")
    ml_stats = json.load(open(ml_prefix + ".stats.json"))
    anc_stats = json.load(open(anc_prefix + ".stats.json"))
    return {
        "cid": cid,
        "ml_seq": ml_seq,
        "anc_seq": anc_seq,
        "ml_stats": ml_stats,
        "anc_stats": anc_stats,
        "recomputed": recomputed,
        "runtime_s": round(time.time() - t0, 1),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clades-root", default=CLADES_ROOT)
    ap.add_argument("--clades", default=TIGHT_CLADES,
                    help="tight_clades json(.gz): {clade_id: [members]}")
    ap.add_argument("--prophage-manifest", default=PROPHAGE_MANIFEST)
    ap.add_argument("--out", default=OUT_DIR)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--n-samples", type=int, default=5,
                    help="traversal draws per clade (v2 used 5)")
    ap.add_argument("--seed", type=int, default=42, help="RNG seed (v1/v2: 42)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-resume", action="store_true",
                    help="recompute even if per-clade outputs already exist")
    ap.add_argument("--clade", action="append",
                    help="restrict to these clade ids (debug)")
    args = ap.parse_args()

    clades = load_clades(args.clades)
    cids = args.clade or sorted(clades.keys())
    os.makedirs(args.out, exist_ok=True)

    alignable = [c for c in cids if len(clades[c]) >= 2]
    singletons = [c for c in cids if len(clades[c]) == 1]
    print(f"[build] clades={len(cids)} alignable={len(alignable)} "
          f"singletons={len(singletons)}", flush=True)

    # every clade must have a prepared clade dir (no clade silently missing)
    missing_dirs = [c for c in cids
                    if not os.path.isdir(os.path.join(args.clades_root, c))]
    if missing_dirs:
        print(f"[build] FATAL: {len(missing_dirs)} clade dirs missing: "
              f"{missing_dirs[:10]}", flush=True)
        return 1

    species_of = load_species(args.prophage_manifest)

    ml_rows = []      # manifest rows
    ml_genomes = []   # (header, seq)
    anc_genomes = []  # (header, seq)
    stats_rows = []   # per-clade traversal stats
    warnings = []
    failures = []

    if args.dry_run:
        n_skip = sum(1 for c in alignable
                     if outputs_complete(os.path.join(args.clades_root, c,
                                                       "ml"), ML_OUTPUTS)
                     and outputs_complete(os.path.join(args.clades_root, c,
                                                       "anc"), ANC_OUTPUTS))
        print(f"[dry-run] would traverse {len(alignable)} alignable clades "
              f"x2 modes (jobs={args.jobs}, n-samples={args.n_samples}, "
              f"seed={args.seed}, {n_skip} already complete/resume-skipped) "
              f"and pass through {len(singletons)} singletons")
        return 0

    python = sys.executable
    resume = not args.no_resume
    tasks = [(cid, os.path.join(args.clades_root, cid), args.n_samples,
              args.seed, python, resume) for cid in alignable]
    t_start = time.time()
    n_skipped_clades = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(process_clade, t): t[0] for t in tasks}
        done = 0
        for fut in concurrent.futures.as_completed(futs):
            cid = futs[fut]
            try:
                res = fut.result()
            except Exception as e:
                # retry once in-process (transient subprocess failures)
                try:
                    res = process_clade((cid, os.path.join(args.clades_root,
                                                           cid),
                                         args.n_samples, args.seed,
                                         python, False))
                except Exception as e2:
                    failures.append((cid, f"{e} | retry: {e2}"))
                    print(f"  [FAIL] {cid}: {e} | retry: {e2}", flush=True)
                    continue
            done += 1
            if res["recomputed"] == 0:
                n_skipped_clades += 1
            n = len(clades[cid])
            cdir = os.path.join(args.clades_root, cid)
            med = statistics.median(member_lengths(cdir))
            ml_len = len(res["ml_seq"])
            anc_len = len(res["anc_seq"])
            flags = []
            if ml_len < MIN_LEN:
                flags.append("too_short")
                warnings.append(
                    f"{cid}: ML genome {ml_len} bp < {MIN_LEN} bp minimum "
                    f"(clade median member {med} bp)")
            if ml_len < HALF_MEDIAN_FACTOR * med:
                flags.append("below_half_median")
                warnings.append(
                    f"{cid}: ML genome {ml_len} bp < half of clade median "
                    f"prophage length ({med / 2:.0f} bp)")
            if anc_len < MIN_LEN:
                flags.append("anc_too_short")
                warnings.append(
                    f"{cid}: ancestral genome {anc_len} bp < {MIN_LEN} bp "
                    f"minimum")
            hosts = sorted({species_of.get(member_acc(m), "?")
                            for m in clades[cid]})
            ml_rows.append((cid, n, "ml", "yes", ml_len, med,
                            "; ".join(hosts), ";".join(clades[cid]),
                            f"{cid}/ml.ml.fa", f"{cid}/anc.ancestral.genome.fa",
                            ",".join(flags)))
            stats_rows.append((cid, n,
                                res["ml_stats"]["n_partitions_total"],
                                res["ml_stats"]["n_partitions_with_seq"],
                                ml_len, res["ml_stats"]["best_sample_index"],
                                res["anc_stats"]["n_partitions_total"],
                                res["anc_stats"]["n_partitions_with_seq"],
                                anc_len,
                                res["anc_stats"]["best_sample_index"],
                                res["recomputed"], res["runtime_s"]))
            ml_genomes.append(
                (f"ntm3_{cid}_ML status=ml n_members={n} length={ml_len}",
                 res["ml_seq"]))
            anc_genomes.append(
                (f"ntm3_{cid}_ANCESTRAL status=ml n_members={n} "
                 f"length={anc_len}", res["anc_seq"]))
            if done % 50 == 0 or done == len(tasks):
                print(f"  [progress] {done}/{len(tasks)} alignable clades "
                      f"traversed ({time.time()-t_start:.0f}s, "
                      f"{n_skipped_clades} resume-skipped)", flush=True)

    # singletons: pass through as-is (v1/v2 convention)
    for cid in singletons:
        cdir = os.path.join(args.clades_root, cid)
        med = statistics.median(member_lengths(cdir))
        _, seq = first_record(os.path.join(cdir, "sequences.fa"))
        flags = []
        if len(seq) < MIN_LEN:
            flags.append("too_short")
            warnings.append(f"{cid}: singleton genome {len(seq)} bp "
                            f"< {MIN_LEN} bp minimum")
        hosts = sorted({species_of.get(member_acc(m), "?")
                        for m in clades[cid]})
        ml_rows.append((cid, 1, "singleton", "n/a", len(seq), med,
                        "; ".join(hosts), ";".join(clades[cid]),
                        f"{cid}/sequences.fa", "n/a", ",".join(flags)))
        stats_rows.append((cid, 1, 0, 0, len(seq), "", 0, 0, 0, "", 0, 0.0))
        ml_genomes.append(
            (f"ntm3_{cid}_ML status=singleton n_members=1 length={len(seq)}",
             seq))

    # sort rows/genomes by header (stable, deterministic: ntm3_<cid>_...)
    ml_genomes.sort(key=lambda x: x[0])
    anc_genomes.sort(key=lambda x: x[0])
    ml_rows.sort(key=lambda x: x[0])
    stats_rows.sort(key=lambda x: x[0])

    # ── write aggregates ────────────────────────────────────────────────────
    ml_fa = os.path.join(args.out, "all_ntm_v3_ml_phage_genomes.fa")
    with open(ml_fa, "w") as f:
        for hdr, seq in ml_genomes:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")

    anc_fa = os.path.join(args.out, "all_ntm_v3_ancestral_phage_genomes.fa")
    with open(anc_fa, "w") as f:
        for hdr, seq in anc_genomes:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")

    manifest = os.path.join(args.out, "release_manifest.tsv")
    with open(manifest, "w") as f:
        f.write("clade_id\tn_members\tstatus\tancestral\tlength_bp"
                "\tmedian_member_len_bp\thost_species\tsource_prophages"
                "\tml_genome_file\tancestral_genome_file\tflags\n")
        for row in ml_rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    stats_path = os.path.join(args.out, "per_clade_stats.tsv")
    with open(stats_path, "w") as f:
        f.write("clade_id\tn_members\tml_partitions_total"
                "\tml_partitions_with_seq\tml_length_bp\tml_best_sample"
                "\tanc_partitions_total\tanc_partitions_with_seq"
                "\tanc_length_bp\tanc_best_sample\tmodes_recomputed"
                "\truntime_s\n")
        for row in stats_rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    warn_path = os.path.join(args.out, "warnings.tsv")
    with open(warn_path, "w") as f:
        f.write("warning\n")
        for w in warnings:
            f.write(w + "\n")

    # ── report ──────────────────────────────────────────────────────────────
    ml_lens = [len(s) for _, s in ml_genomes]
    anc_lens = [len(s) for _, s in anc_genomes]
    n_short = sum(1 for r in ml_rows if "too_short" in r[10])
    n_halfmed = sum(1 for r in ml_rows if "below_half_median" in r[10])
    report = [
        "# NTM v3 — ML + ancestral traversal report",
        "",
        f"Generated: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} "
        f"(task ntm-v3-ml)",
        "",
        "## Run",
        "",
        f"`scripts/traverse_partitions.py` in both modes (`--mode ml`, "
        f"`--mode ancestral`) per alignable clade; `--n-samples "
        f"{args.n_samples} --seed {args.seed}` (v2 convention); "
        f"`--jobs {args.jobs}`; singletons passed through as-is (v1/v2 "
        f"convention). Driver: `ntm/v3/scripts/build_v3_ml_genomes.py` "
        f"(resume-capable; {n_skipped_clades} alignable clades reused "
        f"complete per-clade outputs).",
        "",
        "## Coverage",
        "",
        f"- Total prophage clades: {len(cids)}",
        f"- Alignable (n>=2, traversed x2 modes): {len(alignable)}",
        f"- Singletons (pass-through): {len(singletons)}",
        f"- ML genomes written: {len(ml_genomes)} "
        f"(== alignable + singletons: {len(alignable) + len(singletons)})",
        f"- Ancestral genomes written: {len(anc_genomes)} "
        f"(== alignable: {len(alignable)})",
        f"- Traversal failures: {len(failures)}",
        "",
        "## Length sanity",
        "",
        f"- ML lengths: n={len(ml_lens)} "
        + (f"min={min(ml_lens)} median={int(statistics.median(ml_lens))} "
           f"mean={int(statistics.mean(ml_lens))} max={max(ml_lens)}"
           if ml_lens else "(none)"),
        f"- Ancestral lengths: n={len(anc_lens)} "
        + (f"min={min(anc_lens)} median={int(statistics.median(anc_lens))} "
           f"mean={int(statistics.mean(anc_lens))} max={max(anc_lens)}"
           if anc_lens else "(none)"),
        f"- Genomes < 1 kb (MIN_LEN={MIN_LEN}): {n_short}",
        f"- ML genomes < half clade median member length: {n_halfmed}",
        f"- Warnings logged: {len(warnings)}",
        "",
    ]
    if failures:
        report.append("## Failures")
        report.append("")
        for cid, why in failures:
            report.append(f"- `{cid}`: {why}")
        report.append("")
    report.append("## Outputs")
    report.append("")
    report.append(f"- `{ml_fa}`")
    report.append(f"- `{anc_fa}`")
    report.append(f"- `{manifest}`")
    report.append(f"- `{stats_path}`")
    report.append(f"- `{warn_path}`")
    report.append(f"- per-clade intermediates: `<clade>/ml.*`, `<clade>/anc.*`")
    report.append("")
    report.append(f"Completed in {time.time() - t_start:.0f}s")
    report.append("")
    txt = "\n".join(report)
    print(txt, flush=True)

    report_path = os.path.join(args.out, "ml_report.md")
    with open(report_path, "w") as f:
        f.write(txt.rstrip("\n") + "\n")

    n_ml = sum(1 for r in ml_rows if r[2] == "ml")
    n_sing = sum(1 for r in ml_rows if r[2] == "singleton")
    print(f"[build] DONE ml={n_ml} singleton={n_sing} "
          f"ancestral={len(anc_genomes)} failures={len(failures)} "
          f"warnings={len(warnings)} skipped={n_skipped_clades} "
          f"elapsed={time.time()-t_start:.0f}s", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
