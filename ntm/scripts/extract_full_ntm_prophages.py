#!/usr/bin/env python3
"""
extract_full_ntm_prophages.py — extract every NTM prophage (geNomad provirus.tsv
coords -> samtools faidx) into a single FASTA, single-line sequences (the format
per_clade_alignment_pipeline.py's offset-index reader expects).

Output: ntm/v1/full_prophages.fa  (headers = {acc}_prophage_{N}, matches ntm_prophages.csv)
"""
import csv
import glob
import os
import subprocess
import time

ROOT = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1"
PG = f"{ROOT}/prophages/per_genome"
GENOMES = f"{ROOT}/genomes/canonical_objects"
OUT = f"{ROOT}/full_prophages.fa"


def main():
    t0 = time.time()
    n = ok = 0
    lens = []
    with open(OUT, "w") as fout:
        for pv in sorted(glob.glob(f"{PG}/*_genomad/*_find_proviruses/*_provirus.tsv")):
            acc = os.path.basename(pv).replace(".pansn_provirus.tsv", "")
            gz = f"{GENOMES}/{acc}/{acc}.pansn.fa.gz"
            if not os.path.exists(gz):
                continue
            for i, r in enumerate(csv.DictReader(open(pv), delimiter="\t"), 1):
                n += 1
                seq_name = r["source_seq"]
                beg, end = int(r["start"]), int(r["end"])
                if beg > end:
                    beg, end = end, beg
                region = f"{seq_name}:{beg}-{end}"
                res = subprocess.run(
                    ["samtools", "faidx", gz, region],
                    capture_output=True, text=True)
                if res.returncode != 0 or not res.stdout.strip():
                    continue
                lines = res.stdout.splitlines()
                seq = "".join(l for l in lines[1:] if l and not l.startswith(">"))
                if not seq:
                    continue
                pid = f"{acc}_prophage_{i}"
                fout.write(f">{pid}\n{seq}\n")
                ok += 1
                lens.append(len(seq))
                if ok % 1000 == 0:
                    print(f"  {ok} extracted ({time.time()-t0:.0f}s)", flush=True)
    import statistics
    print(f"extracted {ok}/{n} prophages -> {OUT}", flush=True)
    print(f"length: min {min(lens)} median {int(statistics.median(lens))} "
          f"mean {int(statistics.mean(lens))} max {max(lens)} total {sum(lens)}",
          flush=True)


if __name__ == "__main__":
    main()
