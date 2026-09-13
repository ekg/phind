# NTM v3 — prophage MASH + tight clades report

Generated: 2026-09-13T17:55:20Z (task ntm-v3-prophage)

## Prophage set

- `full_prophages.fa`: **9446 prophages**, 206,284,318 bp (NVMe `ntm/v3/full_prophages.fa`, sha256 `5adfd2aee6e1e03d109a49a3b1fbba837a21649851fc39ee9c9ff07217f0480c` — matches the extract_report receipt)
- provenance: task ntm-v3-unified unified manifest (BV-BRC phigaro primary, v2 calls only for genomes absent from the export); 9,446 of 36,940 manifest rows are extractable — the 27,494 run-assembly rows have no FASTA on this host yet (collaborator delivery pending, `ntm/v3/download_report.md`), none dropped
- Mash sketch: `-i -k 21 -s 10000` (identical to v1/v2)

## Commands (reproducible)

```bash
# driver (every stage resumable from disk state; heavy phase run detached under nohup setsid)
python3 ntm/v3/scripts/full_ntm_mash_clades_v3.py --phase heavy   # ids, sketch, triangle, f32, labels, tree, index
python3 ntm/v3/scripts/full_ntm_mash_clades_v3.py --phase light   # clades, summary, spotcheck, report
# the two mash invocations the driver runs (cwd .../ntm/v3/mash_clades):
mash sketch -i -k 21 -s 10000 -p 64 -o v3_prophages ../full_prophages.fa   # -> v3_prophages.msh
mash triangle -k 21 -s 10000 -p 64 v3_prophages.msh > v3_prophages.triangle.txt
# driver-internal conversion, then clustering with the reused-as-is E. coli machinery (unmodified):
python3 scripts/build_tight_clades.py --threshold 0.25 --max-size 100 --communities 0 \
    --outdir /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/clades --ids-file /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/mash_clades/ids.txt \
    --triangle /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/mash_clades/v3_prophages.dist --labels-csv /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/mash_clades/labels.csv
# independent re-validation
python3 ntm/v3/scripts/validate_mash_clades_v3.py
```

## Clades (threshold 0.25, max-size 100, community 0)

- clade count: **813**
- alignable (>=2 members): **456**
- singletons: **357** (43.9% of clades; each singleton is its own clade — v1/v2 convention, pass-through, no ML reconstruction for n=1)
- assignment check: sum of clade sizes = 9446 == 9446 prophages, all member ids unique (every prophage assigned exactly once)
- max clade size: 100 (cap 100)

### Clade size distribution

| size band | clades |
|---|---:|
| 1 | 357 |
| 2-10 | 310 |
| 11-50 | 81 |
| 51-100 | 65 |
| >100 | 0 |

Top 10 largest clades: `0_0000` (n=100, median=0.1491), `0_0001` (n=100, median=0.0198), `0_0002` (n=100, median=0.1027), `0_0003` (n=100, median=0.1751), `0_0004` (n=100, median=0.0497), `0_0005` (n=100, median=0.0629), `0_0006` (n=100, median=0.1716), `0_0007` (n=100, median=0.0417), `0_0008` (n=100, median=0.0440), `0_0009` (n=100, median=0.1352)

## Internal similarity (per-clade pairwise MASH)

- median of per-clade median distances: **0.0669**
- per-clade median distribution: min 0.000000, Q1 0.000633, median 0.066909, Q3 0.173604, max 0.248966
- clades with internal median exactly 0.0: 38 (identical/near-identical prophages recurring across isolates of the same species)
- clades with median <= 0.25: 456 of 456 alignable (build_tight_clades.py `tighten_clades` enforces the median criterion; any violation is flagged here)
- clades with median > 0.25: **0**

## Prophage length band check

- lengths (bp, recomputed from `full_prophages.fa` by this task): min 237, median 18,276, mean 21,838, max 93,713, total 206,284,318
- outside [1,000, 100,000] bp: **91** (91 below 1,000, 0 above 100,000) — matches the extract_report flag count (91) and its min/max (237 / 93,713); all retained (kept, not dropped — v2/v3 extract convention)

## Comparison with v1 / v2

| version | prophages | clades | alignable | singletons | median internal |
|---|---:|---:|---:|---:|---:|
| v1 (geNomad, 7,352 NTM assemblies) | 10,438 | 913 | 472 | 441 | 0.074 |
| v2 (collab manifest, NCBI-only) | 8,502 | 2,388 | 1,251 | 1,137 | 0.0060 |
| **v3 (unified manifest, extractable set)** | **9446** | **813** | **456** | **357** | **0.0669** |

Cohort note: the task brief anticipated "~4-6x more prophages than v2", but the extractable v3 set is 9446 = 1.11x v2's 8,502 — the other 27,494 manifest rows are run-assembly prophages whose FASTAs are blocked pending collaborator delivery (ntm/v2/run_assemblies/REQUEST.md PENDING). Clade count 813 is in the plausible range for this input size (v2: 2,388 clades / 8,502 prophages; v3 adds BV-BRC phigaro calls for genomes v2 never held). When the run assemblies land, re-running this driver over the extended FASTA reproduces the identical recipe.

## Mechanics note (triangle row order)

`ids.txt` is written in `sorted()` order, not FASTA order. `build_tight_clades.py` (reused as-is) sorts community members with `sorted(set(members))` and its `read_community_matrix` documents *"members must be sorted by triangle row index"*; it fills D[i][j] only for member pairs whose member-list order agrees with triangle row order. In the v1/v2 invocations ids.txt was in FASTA order, so 35.2% of v2's within-community pairs were never filled (NaN, measured on the v2 triangle: 23,407,628 of 36,137,751) and silently behaved as maximally distant — missed joins, hence fragmented (inflated) clade counts. Writing ids.txt pre-sorted satisfies the contract with the script unmodified and fills 100% of pairs (verified: no NaN in the community matrix). All downstream consumers join by prophage id, so the row order is internal to this step. v1/v2 clade counts are therefore upper bounds.

## Validation

- every prophage assigned to exactly one tight clade: sum of clade sizes 9446 == FASTA records 9446 (9446 unique ids) — PASS
- clade internal similarity: every non-singleton clade median pairwise MASH <= 0.25: 0 violations — PASS
- v2-comparable stats reported above (clade count, alignable/singleton split, median internal distance, size distribution)
- commands reproducible: exact commands above and in `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/mash_clades/commands.log`; independent re-validation: `python3 ntm/v3/scripts/validate_mash_clades_v3.py`
- distance integrity: 200 sampled pairs triangle-text vs float32 (0 mismatches); 50 pairs re-measured with a fresh `mash sketch -k 21 -s 10000` + `mash dist` (worst |delta| 1e-08, 0 failures)

## Outputs

- NVMe `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/mash_clades/`: `v3_prophages.msh`, `v3_prophages.triangle.txt` (mash triangle text), `v3_prophages.dist` (float32 upper triangle, 44,608,735 values == n*(n-1)/2, 178,434,940 bytes), `ids.txt`, `prophage_lengths.tsv`, `labels.csv` (community 0), `full_prophages.idx.json`, `prophages_tree.nwk`, `tree_stats.json`, `spotcheck.tsv`, `commands.log`, `driver.log`
- NVMe `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/clades/`: `0/` (tight_clades.json, clade_similarity.json, members.json, distances.npz, commands.log), `tight_clades_summary.json`, `clade_summary.tsv`, `alignable_clades.tsv`, `singletons.tsv`
- repo `ntm/v3/clades/`: `tight_clades.json.gz` (clade definitions), `clade_similarity.json.gz` (per-clade internal similarity stats), `clade_summary.tsv`, `alignable_clades.tsv`, `singletons.tsv`, `SHA256SUMS` (gzip artifacts deterministic, gzip -n)
- FASTA / sketch / triangle stay on NVMe (repo holds only code, manifests and small reports — v1/v2 rule)
