#!/usr/bin/env python3
"""
in-silico-e: assign E. coli Clermont phylogroups (A, B1, B2, C, D, E, F,
clade I-V, cryptic) to every genome in the 26k cohort using ClermonTyping
(the unmodified upstream clermont.py, driven per-genome with NCBI BLAST+).

Method (faithful to the ClermonTyping pipeline):
  1. decompress the PanSN BGZF genome FASTA
  2. makeblastdb  (nucl)
  3. blastn of the ClermonTyping primers.fasta vs the genome (-perc_identity 90,
     task blastn, XML output)   [exactly the upstream recommended command]
  4. run upstream research/phylogroups/tool/clermont.py -x <xml>  -> phylogroup

The driver parallelises over genomes (chunked worker pool). For each genome it
records the full marker profile (pcr_products, quadruplex, specific) alongside
the final phylogroup.

Usage:
    python3 run_clermont_26k.py --fasta-glob <glob> --out <tsv> \
        [--threads N] [--scratch DIR] [--limit N] [--only <file-of-accessions>]
"""
import argparse, glob, multiprocessing as mp, os, subprocess, sys, tempfile, time

CLERMONT_PY_MISSING = None

def blast_bin():
    cand = [
        "/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/run-host-structure-1000/host-structure-1000-v1-run/tool-env/bin",
    ]
    for d in cand:
        if os.path.exists(os.path.join(d, "blastn")):
            return d
    raise RuntimeError("no blast tool-env found")

BLAST_ENV = blast_bin()

def argp():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta-glob", required=True)
    p.add_argument("--out", required=True, dest="tsv")
    p.add_argument("--threads", type=int, default=64)
    p.add_argument("--scratch", default="/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/phylogroup_scratch/pergenome")
    p.add_argument("--clermont-py", required=True)
    p.add_argument("--primers", required=True)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--only", default="", help="file listing accessions (one per line) to process")
    return p.parse_args()

def run_one(acc_fa):
    acc, fa = acc_fa
    clermont_py = os.environ['CLERMONT_PY']
    scratch = os.environ['SCRATCH']
    wd = os.path.join(scratch, acc)
    os.makedirs(wd, exist_ok=True)
    try:
        fasta = os.path.join(wd, "g.fa")
        with open(fasta, "w") as fh:
            r = subprocess.run(["gzip", "-dc", fa], stdout=fh, stderr=subprocess.DEVNULL, check=True)
        db = os.path.join(wd, "gdb")
        subprocess.run([os.path.join(BLAST_ENV,"makeblastdb"), "-in", fasta,
                        "-dbtype", "nucl", "-out", db],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        xml = os.path.join(wd, "g.xml")
        subprocess.run([os.path.join(BLAST_ENV,"blastn"), "-query",
                        os.environ['PRIMERS'], "-db", db, "-task", "blastn",
                        "-perc_identity", "90", "-outfmt", "5", "-out", xml],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        r = subprocess.run([sys.executable, clermont_py, "-x", xml],
                           capture_output=True, text=True, check=True)
        fields = r.stdout.rstrip("\n").split("\t")
        if len(fields) < 4:
            return (acc, "ERROR", r.stdout.strip(), "", "", "")
        pcr_products, quadruplex, specific, phylo = fields[0], fields[1], fields[2], fields[3]
        return (acc, phylo, pcr_products, quadruplex, specific, "")
    except subprocess.CalledProcessError as e:
        return (acc, "ERROR", "", "", "", f"cmd exit {e.returncode}: {e.stderr.strip()[:200]}")
    except Exception as e:
        return (acc, "ERROR", "", "", "", str(e)[:200])

def main():
    args = argp()
    os.environ['CLERMONT_PY'] = args.clermont_py
    os.environ['PRIMERS'] = args.primers
    os.environ['SCRATCH'] = args.scratch
    os.makedirs(args.scratch, exist_ok=True)

    fas = sorted(glob.glob(args.fasta_glob))
    pairs = []
    for fa in fas:
        acc = os.path.basename(fa).split(".pansn.fa.gz")[0]
        pairs.append((acc, fa))
    if args.only:
        only = set(l.strip() for l in open(args.only) if l.strip())
        pairs = [p for p in pairs if p[0] in only]
    if args.limit:
        pairs = pairs[:args.limit]
    print(f"[info] {len(pairs)} genomes to process, threads={args.threads}", file=sys.stderr)

    t0 = time.time()
    results = []
    if args.threads > 1:
        with mp.Pool(args.threads) as pool:
            results = pool.map(run_one, pairs, chunksize=8)
    else:
        results = [run_one(p) for p in pairs]
    wall = time.time() - t0

    with open(args.tsv, "w") as fh:
        fh.write("accession\tphylogroup\tpcr_products\tquadruplex\tspecific\tnote\n")
        for r in results:
            fh.write("\t".join(r) + "\n")

    n = len(results)
    ok = sum(1 for r in results if r[4] != "" and r[1] != "ERROR")
    assigned = sum(1 for r in results if r[1] not in ("ERROR", "") and r[1] not in ("Non Escherichia","Unknown","Non Escherichia"))
    from collections import Counter
    phylo_counter = Counter(r[1] for r in results)
    print(f"[result] processed {n} in {wall:.1f}s", file=sys.stderr)
    print(f"[result] success (clermont returned) {ok}/{n}", file=sys.stderr)
    print(f"[result] phylogroup distribution: {dict(phylo_counter)}", file=sys.stderr)
    print(f"[done] wrote {args.tsv}", file=sys.stderr)

if __name__ == "__main__":
    main()