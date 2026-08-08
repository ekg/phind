#!/usr/bin/env python3
"""
build_ntm_release.py — ML + ancestral traversal over all NTM tight clades,
assemble the full release FASTA, and JOIN by NTM host clade -> the deliverable:
ML phage genomes labeled per NTM clade, with host-range.

Inputs:
  ntm/v1/mash_clades/clades/0/tight_clades.json   (913 clades)
  ntm/v1/mash_clades/clades/<cid>/{partitions.bed, partitions/, sequences.fa}
  ntm/v1/host_clades/host_clades.tsv              (accession -> host_clade_id, species)

Outputs (ntm/v1/release/):
  all_ntm_ml_phage_genomes.fa        one ML genome per clade (header carries host clades)
  all_ntm_ancestral_phage_genomes.fa one ancestral genome per partitioned clade
  release_manifest.tsv               per-clade metadata
  host_range.tsv                     ml_genome -> host_clade_ids + species it spans
  per_ntm_clade_summary.tsv          host_clade -> # ml genomes, # prophages
"""
import csv
import glob
import json
import os
import statistics
import subprocess

ROOT = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1"
CL = f"{ROOT}/mash_clades/clades"
HOST = f"{ROOT}/host_clades/host_clades.tsv"
REL = f"{ROOT}/release"
TRAV = "/home/erikg/phind/scripts/traverse_partitions.py"


def first_seq(fasta):
    name, seq = None, []
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
    return r.returncode


def main():
    os.makedirs(REL, exist_ok=True)
    tc = json.load(open(f"{CL}/0/tight_clades.json"))
    # host lookup
    host = {}
    with open(HOST) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            host[row["accession"]] = (row["host_clade_id"], row["species"])

    ml_rows, anc_rows, combined, ancestral = [], [], [], []
    per_host = {}  # host_clade_id -> [ml ids]
    host_range = []
    n_ml = n_anc = n_sing = 0

    for cid, members in sorted(tc.items()):
        cdir = f"{CL}/{cid}"
        bed = f"{cdir}/partitions.bed"
        pdir = f"{cdir}/partitions"
        has_part = os.path.exists(bed) and len(glob.glob(f"{pdir}/*.maf")) > 0
        src_genomes = sorted({m.rsplit("_prophage_", 1)[0] for m in members})
        hc_ids = sorted({host.get(g, ("?", "?"))[0] for g in src_genomes})
        hc_spp = sorted({host.get(g, ("?", "?"))[1] for g in src_genomes})
        host_range.append((cid, len(members), ",".join(hc_ids), ",".join(hc_spp)))

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
            seq = ""
            if os.path.exists(mlfa):
                _, seq = first_seq(mlfa)
                n_ml += 1
            else:
                status = "ml_failed"
            if os.path.exists(ancfa):
                _, aseq = first_seq(ancfa)
                ancestral.append((f"ntm_{cid}_ANCESTRAL host_clades={','.join(hc_ids)} species={','.join(hc_spp[:3])}", aseq))
                n_anc += 1
            ml_rows.append((cid, len(members), status, "yes" if os.path.exists(ancfa) else "no",
                            len(seq), len(src_genomes), ",".join(hc_ids), ",".join(hc_spp[:3])))
        else:
            n_sing += 1
            _, seq = first_seq(f"{cdir}/sequences.fa")
            status = "singleton"
            ml_rows.append((cid, len(members), status, "n/a", len(seq), len(src_genomes),
                            ",".join(hc_ids), ",".join(hc_spp[:3])))

        if seq:
            hdr = (f"ntm_{cid}_ML status={status} n_members={len(members)} length={len(seq)} "
                   f"host_clades={','.join(hc_ids)} species={','.join(hc_spp[:3])}")
            combined.append((hdr, seq))
            for hcid in hc_ids:
                per_host.setdefault(hcid, []).append(cid)

    # write FASTAs
    with open(f"{REL}/all_ntm_ml_phage_genomes.fa", "w") as f:
        for hdr, seq in combined:
            f.write(f">{hdr}\n{seq}\n")
    with open(f"{REL}/all_ntm_ancestral_phage_genomes.fa", "w") as f:
        for hdr, seq in ancestral:
            f.write(f">{hdr}\n{seq}\n")
    with open(f"{REL}/release_manifest.tsv", "w") as f:
        f.write("clade_id\tn_members\tstatus\tancestral\tlength\tn_source_genomes\thost_clades\tspecies\n")
        for r in ml_rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    with open(f"{REL}/host_range.tsv", "w") as f:
        f.write("clade_id\tn_members\thost_clade_ids\tspecies\n")
        for r in host_range:
            f.write("\t".join(str(x) for x in r) + "\n")
    with open(f"{REL}/per_ntm_clade_summary.tsv", "w") as f:
        f.write("host_clade_id\tn_ml_genomes\tn_clades\n")
        for hcid, cids in sorted(per_host.items(), key=lambda kv: -len(kv[1])):
            f.write(f"{hcid}\t{len(set(cids))}\t{len(cids)}\n")

    lens = [len(s) for _, s in combined]
    print(f"clades={len(tc)} ml={n_ml} ancestral={n_anc} singleton={n_sing}", flush=True)
    print(f"ML lengths: n={len(lens)} min={min(lens)} median={int(statistics.median(lens))} "
          f"mean={int(statistics.mean(lens))} max={max(lens)}", flush=True)
    print(f"host clades represented: {len(per_host)}", flush=True)
    print(f"release -> {REL}/", flush=True)


if __name__ == "__main__":
    main()
