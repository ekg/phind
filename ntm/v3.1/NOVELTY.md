# NTM v3.1 — novelty against the actinobacteriophage reference collection

**Date:** 2026-10-07 · **Panel:** PhagesDB `Actinobacteriophages-All.fasta`
(6,060 actinobacteriophage genomes, 378 MB, downloaded 2026-10-07),
sketched per-sequence (`mash sketch -i -s 10000`).

## Method and controls

`mash dist` between the panel sketch and (a) the 1,304 v3.1 ML genomes,
(b) a 3,000-record sample of the raw extracted prophages.

Controls, all passing:

| control | result |
|---|---|
| panel vs itself (per-sequence) | self-match dist 0, shared **10000/10000** |
| known phages vs panel (21 control phages) | best non-self match shared **345–9,990**, dist 0.00–0.13 |
| **our prophages vs each other** (3,000 × 3,000) | median shared **9,902**, median dist **0.000** |
| our prophages vs panel | median shared **4** |

The third control is the important one: our sequences share hashes freely with
related sequences, so the near-zero overlap with the panel is a property of the
data, not a broken comparison. (`mash dist` with a FASTA query merges the file
into a single sketch — the analysis uses `.msh` on both sides.)

## Result

| novelty band (min MASH to any of 6,060 phages) | ML genomes |
|---|---:|
| near (0.01–0.05) | 18 (1.4%) |
| related (0.05–0.15) | 19 (1.5%) |
| distant (0.15–0.30) | 137 (10.5%) |
| far (0.30–0.50) | 868 (66.6%) |
| novel (>0.50) | 262 (20.1%) |

- median min-distance **0.373**; **86.7% (1,130/1,304)** have no close relative (>0.30)
- best-hit shared hashes: median **2**, p25 1, p75 6 (of 10,000)
- closest relatives found: `ntm31_0_0534_ML` ≈ **Dori** (0.027),
  `ntm31_0_0064_ML` ≈ **P3MA** (0.036)
- the same result holds on the **raw prophages** (median 4 shared hashes), so it is
  not a reconstruction artifact

Novelty by clade size is flat (90% / 89% / 84% / 80% for singleton → large), so it
is not driven by short fragments alone.

## Sensitivity to k (is the divergence a k-mer artifact?)

| | k=21 | **k=13** (more sensitive) |
|---|---:|---:|
| median best-hit distance | 0.373 | **0.347** |
| median best-hit shared hashes | 4 | 55 |
| control: our prophages vs each other | 9,902 | 9,903 |

Lowering k recovers more shared hashes (as expected) but the **best-hit distance
barely moves** (0.373 → 0.347), and the within-dataset control stays near saturation.
The divergence is therefore not a k-mer sensitivity artifact.

## Interpretation and caveats

**Claim supported:** the NTM prophage population in this cohort is
**nucleotide-divergent from the curated actinobacteriophage collection** — 87% have
no relative within MASH 0.30.

**Caveats that must travel with it:**

1. `mash` k=21 detects exact 21-mers, i.e. roughly ≥95% local identity. Low
   shared-hash counts mean *low nucleotide identity*, not necessarily an unrelated
   phage. Protein-level homology is a separate question — and the Pharokka
   annotation does find phage hallmark genes in these genomes (terminase 35%,
   capsid 29%, tail 47%), so they are phage-like at the protein level.
2. **PhagesDB is biased** toward phages isolated on fast-growing hosts
   (M. smegmatis) and toward curated, cluster-assigned isolates. Our cohort is
   dominated by *Mycobacteroides abscessus* (76% of prophages), whose phages are
   poorly represented. Some of the apparent novelty is reference-coverage bias.
3. These are **integrated prophages**, whereas most PhagesDB entries are isolated
   phage isolates — a different sampling frame.
4. A definitive novelty claim should be made at the **protein/gene-content** level
   (e.g. PHROG/pham profile against the panel), not nucleotides alone.

**Recommended phrasing for the release:** "87% of the reconstructed prophage
lineages have no close nucleotide relative (MASH > 0.30) among 6,060 curated
actinobacteriophages; this is consistent with substantial novelty, but the
reference collection is biased toward fast-growing-host isolates, so
protein-level comparison is required before claiming biological novelty."
