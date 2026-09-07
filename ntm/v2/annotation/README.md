# NTM v2 per-genome GFF3 files (Pharokka)

Per-phage GFF3 annotation for all **3,639** NTM v2 genomes (2,388 ML + 1,251
ancestral), split from the single merged `pharokka.gff` of the 2026-08-17
Pharokka 1.10.1 run (pyrodigal-gv 0.3.2 CDS calls; PHROG/terL annotation in
attributes). Each file: `##gff-version 3` + one `##sequence-region` + all
CDS features for that genome. Embedded `##FASTA` omitted — genome sequences
live in `ml/all_ntm2_ml_phage_genomes.fa` (NVMe root, see RELEASE).

## Locations

- Per-genome GFFs: `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/annotation/pharokka_out/per_genome_gff/<genome_id>.gff`
- Single-file tarball (4.0 MB, all 3,639):
  `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/annotation/pharokka_out/ntm_v2_per_genome_gff.tar.gz`
  sha256 `fc95e9c52021eb1df879570d1aad564e84697b4399793e22e2d908eb7da777c1`
- Index committed here: `manifest.tsv` (genome, features, bytes, sha256),
  `manifest_meta.json` (totals: 3,639 files, 206,249 CDS, 0 empty).

Full merged outputs (GFF/GBK/TBL, pharokka CDS tables, prodigal FFN/FAA) stay
on NVMe under `.../annotation/pharokka_out/`.
