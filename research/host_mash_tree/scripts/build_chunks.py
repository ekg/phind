#!/usr/bin/env python3
"""build_chunks.py — build chunk fastas + ids.txt directly from the 26k PanSN
bgzip genomes, one sketch/leaf per assembly accession.

Source: /mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/26k/canonical_objects/<ACC>/*.pansn.fa.gz
Each PanSN file is one genome, possibly multi-contig
(header `>ACC#hap#chr`). We concatenate all contigs of a genome into a single
record named by the assembly accession (GCF_/GCA_, stripped of `#...`).

Output:
  workdir/chunks/chunk_NNN.fa  (CHUNK_SIZE genomes each, last partial)
  workdir/../ids.txt           (accession per line; line index = row/col index)
  workdir/all.txt              (list of accessions as sketched, same order)
"""
import argparse, gzip, glob, os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True,
                    help="canonical_objects dir containing ACC/*.pansn.fa.gz")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--ids", required=True)
    ap.add_argument("--chunk-size", type=int, default=5000)
    args = ap.parse_args()

    pansn = sorted(glob.glob(os.path.join(args.src, "*", "*.pansn.fa.gz")))
    print(f"[build] {len(pansn)} PanSN genomes", flush=True)

    # dedupe by directory (each genome dir has exactly one pansn file)
    bydir = {}
    for p in pansn:
        acc = os.path.basename(os.path.dirname(p))
        bydir[acc] = p
    accs = sorted(bydir.keys())
    print(f"[build] {len(accs)} unique accessions", flush=True)

    os.makedirs(os.path.join(args.workdir, "chunks"), exist_ok=True)
    with open(args.ids, "w") as f:
        for a in accs:
            f.write(a + "\n")

    # build chunk fastas
    chunks = {}
    buff = []
    for acc in accs:
        if len(buff) >= args.chunk_size:
            chunks[len(chunks)] = buff
            buff = []
        buff.append(acc)
    if buff:
        chunks[len(chunks)] = buff
    nchunks = len(chunks)
    print(f"[build] {nchunks} chunks of up to {args.chunk_size}", flush=True)

    total_len = 0
    for ci, caccs in sorted(chunks.items()):
        cf = os.path.join(args.workdir, "chunks", f"chunk_{ci:03d}.fa")
        with open(cf, "w") as out:
            for acc in caccs:
                p = bydir[acc]
                # concatenate all contig sequences
                seq_chars = []
                with gzip.open(p, "rt") as fh:
                    for line in fh:
                        if line.startswith(">"):
                            continue
                        seq_chars.append(line.rstrip())
                seq = "".join(seq_chars).upper()
                total_len += len(seq)
                out.write(f">{acc}\n")
                # wrap at 100 for readability; mash doesn't need specific wrap
                for i in range(0, len(seq), 100):
                    out.write(seq[i:i+100] + "\n")
        print(f"[build] chunk {ci:03d}: {len(caccs)} genomes -> {cf}", flush=True)

    print(f"[build] total sequence length (all genomes) ~{total_len/1e6:.1f} Mbp", flush=True)

if __name__ == "__main__":
    main()