#!/usr/bin/env python3
"""
abscessus_ml_genomes.py — run ML + ancestral traversal over all abscessus
tight clades -> first NTM ML phage genomes + release FASTA/manifest.

Per clade:
  - partitions present: scripts/traverse_partitions.py --mode ml (and ancestral)
  - singleton / no partitions: ML genome = the member sequence itself
Assemble all -> subset_abscessus/release/all_abscessus_ml_phage_genomes.fa
                  + release_manifest.tsv
"""
import glob
import json
import os
import subprocess

ROOT = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1/subset_abscessus"
TRAV = "/home/erikg/phind/scripts/traverse_partitions.py"
tc = json.load(open(f"{ROOT}/clades/0/tight_clades.json"))
REL = f"{ROOT}/release"
os.makedirs(REL, exist_ok=True)


def first_seq(fasta):
    seq = []
    name = None
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].strip()
            else:
                seq.append(line.strip())
    return name, "".join(seq)


def run(args):
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  WARN {args[-3:]}: {r.stderr[-300:]}", flush=True)
    return r.returncode


rows = []
combined = []
n_ml = n_anc = n_sing = 0
for cid, members in sorted(tc.items()):
    cdir = f"{ROOT}/clades/{cid}"
    bed = f"{cdir}/partitions.bed"
    pdir = f"{cdir}/partitions"
    has_part = os.path.exists(bed) and len(glob.glob(f"{pdir}/*.maf")) > 0
    src_genomes = ",".join(sorted({m.rsplit("_prophage_", 1)[0] for m in members}))
    if has_part:
        run(["python3", TRAV, "--partitions-dir", pdir, "--bed", bed,
             "--output", f"{cdir}/ml", "--mode", "ml", "--n-samples", "5",
             "--seed", "42"])
        run(["python3", TRAV, "--partitions-dir", pdir, "--bed", bed,
             "--output", f"{cdir}/anc", "--mode", "ancestral", "--n-samples", "5",
             "--seed", "42"])
        mlfa = f"{cdir}/ml.ml.fa"
        ancfa = f"{cdir}/anc.ancestral.genome.fa"
        status = "ml"
        if os.path.exists(mlfa):
            n_ml += 1
            name, seq = first_seq(mlfa)
            combined.append((f"abscessus_{cid}_ML status=ml n_members={len(members)} length={len(seq)} genomes={src_genomes}", seq))
        else:
            status = "ml_failed"
        anc_ok = os.path.exists(ancfa)
        if anc_ok:
            n_anc += 1
        rows.append((cid, len(members), status, "yes" if anc_ok else "no",
                     len(seq) if 'seq' in dir() else 0, src_genomes))
    else:
        # singleton or no-alignment clade: genome = member sequence
        n_sing += 1
        name, seq = first_seq(f"{cdir}/sequences.fa")
        combined.append((f"abscessus_{cid}_ML status=singleton n_members={len(members)} length={len(seq)} genomes={src_genomes}", seq))
        rows.append((cid, len(members), "singleton", "n/a", len(seq), src_genomes))

# write release FASTA
with open(f"{REL}/all_abscessus_ml_phage_genomes.fa", "w") as f:
    for hdr, seq in combined:
        f.write(f">{hdr}\n{seq}\n")
with open(f"{REL}/release_manifest.tsv", "w") as f:
    f.write("clade_id\tn_members\tstatus\tancestral\tlength\tsource_genomes\n")
    for r in rows:
        f.write("\t".join(str(x) for x in r) + "\n")

lens = [len(s) for _, s in combined]
import statistics
print(f"clades={len(tc)} ml={n_ml} ancestral={n_anc} singleton={n_sing}", flush=True)
print(f"genome lengths: n={len(lens)} min={min(lens)} median={int(statistics.median(lens))} "
      f"mean={int(statistics.mean(lens))} max={max(lens)}", flush=True)
print(f"release -> {REL}/all_abscessus_ml_phage_genomes.fa", flush=True)
