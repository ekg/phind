#!/usr/bin/env python3
"""
build_v2_ml_genomes.py — NTM v2 ML + ancestral traversal over every prophage
clade (task ntm-v2-ml).

Runs ``scripts/traverse_partitions.py`` in both modes (``--mode ml`` and
``--mode ancestral``, seed 42 as v1) over every alignable (partitioned) clade
under ``ntm/v2/clades/``, and passes singletons through as-is (v1
convention: first record of ``sequences.fa``, ``status=singleton``).

Per-clade intermediates (v1 convention) are written into the clade dir:
  <clade>/ml.ml.fa, ml.consensus.fa, ml.traversal.json, ml.coverage.tsv,
  ml.stats.json  (--mode ml)
  <clade>/anc.ancestral.genome.fa, anc.ancestral.fa, anc.trees.nwk, ...
  (--mode ancestral)

Aggregates are written to ntm/v2/ml/:
  all_ntm2_ml_phage_genomes.fa        one ML genome per clade
                                      (alignable + singletons = 2388)
  all_ntm2_ancestral_phage_genomes.fa one ancestral genome per alignable
                                      clade (1251)
  release_manifest.tsv                per-clade metadata + sanity flags
  ml_report.md                        run report
  warnings.tsv                        length-sanity warnings (short genomes)

Header format (v1 convention, ntm2 prefix):
  >ntm2_<clade_id>_ML status=ml|singleton n_members=<N> length=<L>
  >ntm2_<clade_id>_ANCESTRAL status=ml n_members=<N> length=<L>

Length sanity (validation):
  * every genome >= 1 kb
  * any genome shorter than half its clade's median prophage length
    (computed from the clade's sequences.fa) is logged as a warning

Determinism: the traversal is seeded (seed 42, as v1); same inputs + same
parameters -> byte-identical outputs.  Run with the same arguments twice and
the output files are identical.

Usage:
  python3 ntm/v2/scripts/build_v2_ml_genomes.py \
      [--clades /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/clades] \
      [--out /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/ml] \
      [--jobs 16] [--n-samples 5] [--seed 42] [--dry-run]
"""
import argparse
import concurrent.futures
import json
import os
import statistics
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
TRAV = os.path.join(REPO, "scripts", "traverse_partitions.py")

V2 = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2"
CLADES_ROOT = f"{V2}/clades"
RESULTS_JSON = f"{CLADES_ROOT}/pipeline_results.json"
OUT_DIR = f"{V2}/ml"
MIN_LEN = 1000            # validation: no genome shorter than 1 kb
HALF_MEDIAN_FACTOR = 0.5  # warn if genome < factor * clade median member len


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
    """Worker: run both traverse modes for one alignable clade.

    Returns (cid, result_dict) where result_dict carries the ML + ancestral
    genome sequences and diagnostics, or raises."""
    cid, cdir, n_samples, seed, python = args
    t0 = time.time()
    ml_prefix = run_traverse(cid, cdir, "ml", n_samples, seed, python)
    anc_prefix = run_traverse(cid, cdir, "ancestral", n_samples, seed, python)
    ml_hdr, ml_seq = first_record(ml_prefix + ".ml.fa")
    anc_hdr, anc_seq = first_record(anc_prefix + ".ancestral.genome.fa")
    return {
        "cid": cid,
        "ml_hdr": ml_hdr,
        "ml_seq": ml_seq,
        "anc_hdr": anc_hdr,
        "anc_seq": anc_seq,
        "runtime_s": round(time.time() - t0, 1),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clades-root", default=CLADES_ROOT)
    ap.add_argument("--results", default=RESULTS_JSON)
    ap.add_argument("--out", default=OUT_DIR)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--n-samples", type=int, default=5,
                    help="traversal draws per clade (v1 used 5)")
    ap.add_argument("--seed", type=int, default=42, help="RNG seed (v1: 42)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--clade", action="append",
                    help="restrict to these clade ids (debug)")
    args = ap.parse_args()

    results = json.load(open(args.results))
    cids = args.clade or sorted(results.keys())
    os.makedirs(args.out, exist_ok=True)

    alignable = [c for c in cids if results[c].get("n_members", 0) >= 2]
    singletons = [c for c in cids if results[c].get("n_members", 0) == 1]
    print(f"[build] clades={len(cids)} alignable={len(alignable)} "
          f"singletons={len(singletons)}", flush=True)

    ml_rows = []     # (cid, n_members, status, has_anc, length, med_len, flags)
    ml_genomes = []  # (header, seq)
    anc_genomes = []  # (header, seq)
    warnings = []
    failures = []

    t_start = time.time()
    if args.dry_run:
        print(f"[dry-run] would traverse {len(alignable)} alignable clades "
              f"x2 modes (jobs={args.jobs}, n-samples={args.n_samples}, "
              f"seed={args.seed}) and pass through {len(singletons)} "
              f"singletons")
        return 0

    python = sys.executable
    tasks = [(cid, os.path.join(args.clades_root, cid), args.n_samples,
              args.seed, python) for cid in alignable]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(process_clade, t): t[0] for t in tasks}
        done = 0
        for fut in concurrent.futures.as_completed(futs):
            cid = futs[fut]
            try:
                res = fut.result()
            except Exception as e:
                # retry once (transient subprocess failures)
                try:
                    res = process_clade((cid, os.path.join(args.clades_root,
                                                           cid),
                                         args.n_samples, args.seed, python))
                except Exception as e2:
                    failures.append((cid, f"{e} | retry: {e2}"))
                    print(f"  [FAIL] {cid}: {e} | retry: {e2}", flush=True)
                    continue
            done += 1
            n = results[cid]["n_members"]
            med = statistics.median(member_lengths(os.path.join(
                args.clades_root, cid)))
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
            ml_rows.append((cid, n, "ml", "yes", ml_len, med,
                            ",".join(flags)))
            ml_genomes.append(
                (f"ntm2_{cid}_ML status=ml n_members={n} length={ml_len}",
                 res["ml_seq"]))
            anc_genomes.append(
                (f"ntm2_{cid}_ANCESTRAL status=ml n_members={n} "
                 f"length={anc_len}", res["anc_seq"]))
            if done % 100 == 0 or done == len(tasks):
                print(f"  [progress] {done}/{len(tasks)} alignable clades "
                      f"traversed ({time.time()-t_start:.0f}s)", flush=True)

    # singletons: pass through as-is (v1 convention)
    for cid in singletons:
        cdir = os.path.join(args.clades_root, cid)
        med = statistics.median(member_lengths(cdir))
        _, seq = first_record(os.path.join(cdir, "sequences.fa"))
        flags = []
        if len(seq) < MIN_LEN:
            flags.append("too_short")
            warnings.append(f"{cid}: singleton genome {len(seq)} bp "
                            f"< {MIN_LEN} bp minimum")
        ml_rows.append((cid, 1, "singleton", "n/a", len(seq), med,
                        ",".join(flags)))
        ml_genomes.append(
            (f"ntm2_{cid}_ML status=singleton n_members=1 length={len(seq)}",
             seq))

    # sort rows/genomes by clade id (stable, deterministic)
    ml_genomes.sort(key=lambda x: x[0])
    anc_genomes.sort(key=lambda x: x[0])
    ml_rows.sort(key=lambda x: x[0])

    # ── write aggregates ────────────────────────────────────────────────────
    ml_fa = os.path.join(args.out, "all_ntm2_ml_phage_genomes.fa")
    with open(ml_fa, "w") as f:
        for hdr, seq in ml_genomes:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")

    anc_fa = os.path.join(args.out, "all_ntm2_ancestral_phage_genomes.fa")
    with open(anc_fa, "w") as f:
        for hdr, seq in anc_genomes:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")

    manifest = os.path.join(args.out, "release_manifest.tsv")
    with open(manifest, "w") as f:
        f.write("clade_id\tn_members\tstatus\tancestral\tlength_bp"
                "\tmedian_member_len_bp\tflags\n")
        for cid, n, status, anc, length, med, flags in ml_rows:
            f.write(f"{cid}\t{n}\t{status}\t{anc}\t{length}\t{med}"
                    f"\t{flags}\n")

    warn_path = os.path.join(args.out, "warnings.tsv")
    with open(warn_path, "w") as f:
        f.write("warning\n")
        for w in warnings:
            f.write(w + "\n")

    # ── report ──────────────────────────────────────────────────────────────
    ml_lens = [len(s) for _, s in ml_genomes]
    anc_lens = [len(s) for _, s in anc_genomes]
    n_short = sum(1 for r in ml_rows if "too_short" in r[6])
    n_halfmed = sum(1 for r in ml_rows if "below_half_median" in r[6])
    report = [
        "# NTM v2 — ML + ancestral traversal report",
        "",
        f"Generated: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} "
        f"(task ntm-v2-ml)",
        "",
        "## Run",
        "",
        f"`scripts/traverse_partitions.py` in both modes (`--mode ml`, "
        f"`--mode ancestral`) per alignable clade; `--n-samples "
        f"{args.n_samples} --seed {args.seed}` (seed 42 as v1); "
        f"`--jobs {args.jobs}`; singletons passed through as-is (v1 "
        f"convention).",
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
    n_anc = len(anc_genomes)
    print(f"[build] DONE ml={n_ml} singleton={n_sing} ancestral={n_anc} "
          f"failures={len(failures)} warnings={len(warnings)} "
          f"elapsed={time.time()-t_start:.0f}s", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
