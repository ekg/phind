# NTM v2 — ML + ancestral traversal report

Generated: 2026-08-16T15:17:24Z (task ntm-v2-ml)

## Run

`scripts/traverse_partitions.py` in both modes (`--mode ml`, `--mode ancestral`) per alignable clade; `--n-samples 5 --seed 42` (seed 42 as v1); `--jobs 16`; singletons passed through as-is (v1 convention).

## Coverage

- Total prophage clades: 2388
- Alignable (n>=2, traversed x2 modes): 1251
- Singletons (pass-through): 1137
- ML genomes written: 2388 (== alignable + singletons: 2388)
- Ancestral genomes written: 1251 (== alignable: 1251)
- Traversal failures: 0

## Length sanity

- ML lengths: n=2388 min=237 median=20268 mean=28284 max=150000
- Ancestral lengths: n=1251 min=237 median=30277 mean=37031 max=150000
- Genomes < 1 kb (MIN_LEN=1000): 23
- ML genomes < half clade median member length: 7
- Warnings logged: 39

## Outputs

- `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/ml/all_ntm2_ml_phage_genomes.fa`
- `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/ml/all_ntm2_ancestral_phage_genomes.fa`
- `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/ml/release_manifest.tsv`
- `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/ml/warnings.tsv`
- per-clade intermediates: `<clade>/ml.*`, `<clade>/anc.*`

Completed in 93s

