# NTM v3.1 — index

`v3.1` is a **method revision of the v3 ML step**, plus the analysis run on top of
it. v3's clades, host clades, partitions and annotation are reused unchanged; v3's
artifacts on disk are untouched.

Start here:

| doc | what |
|---|---|
| [`REPORT.md`](REPORT.md) | what v3.1 changes, the A/B, full-set results, chimerism, honest verdict |
| [`NOVELTY.md`](NOVELTY.md) | novelty vs 6,060 actinobacteriophages, with controls, k-sensitivity and caveats |
| [`SPEC.md`](SPEC.md) | the original change spec (defects, flags, validation protocol) |
| [`../v3/RELEASE.md`](../v3/RELEASE.md) | **the v3 correction** — v3's ML genomes are mosaics; read this first |

## The four things worth knowing

1. **v3's ML genomes were mosaics, and CheckV hid it.** The walk only *biased*
   toward adjacency, so it accumulated blocks to the 150 kb budget. 57% of clades
   exceeded 1.5× the median member length; **387 of 395 CheckV "high-quality"
   genomes were mosaics**, because concatenating every common block includes every
   marker gene. CheckV completeness rises monotonically with the mosaic ratio
   (15.6% at ≤1.5× → 100% at >3×) and is therefore unusable alone as an
   acceptance criterion.
2. **The fix works.** `--path-mode observed` walks only real adjacencies to a
   natural terminus; with a CheckV-derived per-clade target it lands at **1.00×
   the expected genome length** (vs ~2.7× for v3). Median ML length falls from
   ~150 kb to 14.3 kb.
3. **Honest completeness is lower, and that is the point.** 124 high-quality
   genomes versus v3's 395 — but 98% of those 395 were chimeras. High completeness
   and true genome length are mutually exclusive here.
4. **Novelty is substantial but must be claimed carefully.** 86.7% of ML genomes
   have no relative within MASH 0.30 of 6,060 actinobacteriophages (median 2/10,000
   shared hashes at k=21, 55 at k=13; controls validate the machinery). Caveat:
   PhagesDB is biased toward fast-growing-host isolates, and this is a
   *nucleotide* measure — protein-level comparison is still required.

Plus the synthesis-relevant result: **45.7% of reconstructions follow a single
member genome's partition order exactly** (median longest single-member run 85.7%);
the rest join 2+ genome contexts and contain at least one junction no single genome
has. `chimerism_verdict` is carried in the manifest.

## Deliverables

`$NVME/ntm/v3.1/` (`$NVME = /mnt/nvme3n1/erikg/phind-genome-work`):

| artifact | |
|---|---|
| `catalog/all_ntm_v31_ml_phage_genomes.fa.gz` (+ `.fai`/`.gzi`) | 1,304 ML genomes, bgzip + faidx |
| `catalog/all_ntm_v31_ancestral_phage_genomes.fa.gz` | 893 ancestral genomes |
| `catalog/release_manifest.tsv` | per-clade: target, ratio, CheckV expected/completeness/quality, terminated_by, chimerism verdict |
| `chimerism.tsv` | per-clade single-member vs chimeric + longest run |
| `targets.tsv` | per-clade genome-length target + source (CheckV AAI vs fallback) |
| `ab/` | the A/B run (25 clades, free vs observed) |
| `member_checkv/`, `checkv_catalog/` | CheckV over prophages and over the v3.1 catalog |

## Drivers

`ab_subset.py`, `build_v31_targets.py`, `run_v31_ml.py`, `chimerism.py`.
Changes to the shared pipeline are backward-compatible: `--split-mode cap`
(default) and `--path-mode free` (default) reproduce v3 byte-for-byte.

## Open items

1. **697/1,304 clades used the `member_seq_median` fallback target** (0.78× the
   true expected length) — those are ~22% short, and it is most of the mid-clade
   completeness drop. A two-pass retarget (reconstruct → CheckV AAI → retarget)
   is the obvious improvement.
2. **Protein-level novelty** vs the panel not yet done (Pharokka finds phage
   hallmark genes, so they are phage-like at the protein level).
3. The `--split-mode medoids` clustering change is implemented and tested but
   **not applied** — the mosaic defect was in the ML step, not the clustering.
