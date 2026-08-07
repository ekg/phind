#!/usr/bin/env python3
"""
call_ntm_prophages.py — geNomad prophage calling across all NTM genomes.

Resume-safe, parallel: each genome -> genomad end-to-end (find-proviruses on),
8 threads/genome, --jobs concurrent genomes (32 x 8 = 256 cores on this box).
Completion marker per genome = <acc>_find_proviruses/<acc>_provirus.tsv
(created even with 0 prophages, so it is a valid done-signal).

After this finishes, consolidate with workflow/ntm/normalize_genomad.py ->
ntm/v1/prophages/ntm_prophages.csv (E. coli-compatible schema).

Usage:
  export PATH=/home/erikg/micromamba/envs/geomad/bin:$PATH
  export PYTHONNOUSERSITE=1
  python3 ntm/scripts/call_ntm_prophages.py --jobs 32 --threads 8
  python3 ntm/scripts/call_ntm_prophages.py --limit 5   # smoke test
"""
import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

GENOMAD = "/home/erikg/micromamba/envs/geomad/bin/genomad"
DB = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1/geomad_db/genomad_db"
GENOMES = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1/genomes/canonical_objects"
OUT = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1/prophages/per_genome"
ENV_PATH = "/home/erikg/micromamba/envs/geomad/bin"


def call_one(args_tuple):
    acc, threads, timeout = args_tuple
    gz = f"{GENOMES}/{acc}/{acc}.pansn.fa.gz"
    odir = f"{OUT}/{acc}_genomad"
    # geNomad derives output prefix from the INPUT file stem (-> "{acc}.pansn")
    stem = os.path.basename(gz)
    if stem.endswith(".fa.gz"):
        stem = stem[: -len(".fa.gz")]
    marker = f"{odir}/{stem}_find_proviruses/{stem}_provirus.tsv"
    if os.path.exists(marker):
        return acc, "SKIP"
    if not os.path.exists(gz):
        return acc, "NO_FASTA"
    env = dict(os.environ,
               PATH=ENV_PATH + ":" + os.environ.get("PATH", ""),
               PYTHONNOUSERSITE="1")
    try:
        r = subprocess.run(
            [GENOMAD, "end-to-end", gz, odir, DB, "-t", str(threads)],
            capture_output=True, text=True, env=env, timeout=timeout)
        if r.returncode != 0:
            tail = (r.stderr or r.stdout or "")[-300:].replace("\n", " ")
            return acc, f"FAIL({r.returncode}): {tail}"
        if not os.path.exists(marker):
            return acc, "NO_MARKER"
        return acc, "OK"
    except subprocess.TimeoutExpired:
        return acc, "TIMEOUT"
    except Exception as e:
        return acc, f"ERR: {e!r}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jobs", type=int, default=32, help="concurrent genomes")
    ap.add_argument("--threads", type=int, default=8, help="threads per genome")
    ap.add_argument("--timeout", type=int, default=3600, help="per-genome timeout (s)")
    ap.add_argument("--limit", type=int, default=0, help="only first N (smoke test)")
    a = ap.parse_args()

    accs = sorted(d for d in os.listdir(GENOMES)
                  if os.path.exists(os.path.join(GENOMES, d, d + ".pansn.fa.gz")))
    if a.limit:
        accs = accs[:a.limit]
    os.makedirs(OUT, exist_ok=True)
    work = [(acc, a.threads, a.timeout) for acc in accs]
    print(f"{len(accs)} genomes | jobs={a.jobs} threads/genome={a.threads} "
          f"timeout={a.timeout}s", flush=True)

    done = skip = fail = 0
    t0 = time.time()
    faillog = open(os.path.join(OUT, "_failures.log"), "a")
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(call_one, w): w[0] for w in work}
        for i, fut in enumerate(as_completed(futs), 1):
            acc, st = fut.result()
            if st == "OK":
                done += 1
            elif st == "SKIP":
                skip += 1
            else:
                fail += 1
                faillog.write(f"{acc}\t{st}\n"); faillog.flush()
            if i % 50 == 0 or i == len(work):
                el = time.time() - t0
                rate = i / el if el else 0
                eta = (len(work) - i) / rate if rate else 0
                print(f"[{i}/{len(work)}] ok={done} skip={skip} fail={fail} "
                      f"elapsed={el:.0f}s rate={rate:.2f}/s eta={eta/60:.0f}min",
                      flush=True)
    faillog.close()
    print(f"FINISHED ok={done} skip={skip} fail={fail} total={len(work)}", flush=True)
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
