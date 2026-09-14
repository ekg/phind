# NTM v2 — ML + ancestral validation

- ML genomes: 2388 (expected 2388)
- Ancestral genomes: 1251 (expected 1251)
- Manifest rows: 2388
- ML count OK: 2388 (alignable 1251 + singletons 1137)
- ancestral count OK: 1251

> **Erratum (2026-09-13, task verify-v2-clade-order):** the clade counts
> above — 2,388 total (1,251 alignable + 1,137 singletons) — come from the
> frozen v2 clade run, which was inflated by NaN fragmentation (35.23% of
> within-community distance-triangle pairs were never filled and behaved as
> maximally distant). Re-derivation with a pre-sorted `ids.txt` gives **767
> clades = 413 alignable + 354 singletons, median internal mash 0.0639** —
> see [`clades_sorted_report.md`](clades_sorted_report.md). Historical
> numbers in this report are unchanged.

- warnings.tsv: 30 clades with logged warnings
- genomes < 1000 bp: 23 (warned)
- determinism re-run: 10 clades x2 modes, 0 mismatches

**Validation: PASS**
