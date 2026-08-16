# ENA probe findings — run-assembly acquisition (2026-08-16)

## Result

The chat agent's earlier 5/5 → 0-hit probe **does not represent the true
rate**. A deterministic scale probe (n=100, seed=42) of the collaborator
run-assembly cohort returned:

| metric | value |
|---|---|
| runs probed | 100 |
| run → sample pivot OK | 100 (100%) |
| with any ENA analysis | 79 (79%) |
| with retrievable `submitted_ftp` | 79 (79%) |
| pivot/query errors | 0 |

**79% >> the 5% threshold — the hits are itemized below and in
`ena_probe_results.tsv`.**

Method (per task spec): `result=read_run&fields=accession,sample_accession`
(run→sample pivot), then
`result=analysis&query=sample_accession=<sample>&fields=analysis_accession,submitted_ftp`.
Probe script: `ntm/v2/scripts/probe_ena_runs.py`. Raw results:
`ntm/v2/run_assemblies/ena_probe_results.tsv`; sampled accessions:
`ntm/v2/run_assemblies/ena_probe_sample.txt`.

## Key findings

1. **79/100 runs have retrievable ENA submitted assemblies**
   (`analysis_type=SEQUENCE_ASSEMBLY`, `submitted_ftp` confirmed HTTP 200).
   Of the 79 hits, **49 are prophage-bearing** (62% of the 63 sampled
   prophage-bearing runs).
2. **The 21 misses cluster in one study's accessions** — 18/21 are `ERR541*`
   (and 3 others: `ERR5412579`, `SRR11839050`, `ERR5411974`), consistent
   with a single run-group whose samples were not deposited as analyses.
3. **ENA assemblies are NOT the collaborator assemblies.** Spot-checked 4
   ENA FASTA files against the manifest `assembly_contigs`/
   `assembly_length_bp` columns — **none match**:

   | accession | ENA contigs | ENA bp | manifest contigs | manifest bp |
   |---|---|---|---|---|
   | ERR3566243 | 48 | 6,085,029 | 98 | 6,109,636 |
   | ERR2759527 | 48 | 5,198,350 | 23 | 5,193,530 |
   | ERR3177943 | 11 | 5,099,854 | 9 | 5,098,843 |
   | ERR330884 | 61 | 5,207,668 | 12 | 5,209,116 |

   The manifest prophage coordinates (`prophage_contig/start/end`) and the
   `assembly_contigs`/`assembly_length_bp` numbers were computed on the
   collaborator's `contigs.fasta`. Prophage-bearing genome coverage (5,907/9,434)
   depends on the collaborator files; **the collaborator request stands as
   the authoritative source**. ENA is a fallback only for individual
   unavailable genomes, and only with full re-verification (contig set and
   prophage re-calling) before use.

## Itemization (79 retrievable runs)

See machine-readable `ena_probe_results.tsv` (columns:
`run, sample_accession, analysis_accessions, submitted_ftp, has_prophage,
n_hits`). Summary of hit accessions by prefix: 45 ERR + 34 SRR.

<details>
<summary>79 hit run accessions (click to expand)</summary>

```
ERR2759527 ERR3177943 ERR3178009 ERR3198400 ERR3198441 ERR3198446
ERR330884 ERR3337404 ERR3337600 ERR337844 ERR340492 ERR340525 ERR3468808
ERR3468825 ERR3468954 ERR350139 ERR350185 ERR350236 ERR3524842 ERR3566236
ERR3566243 ERR3566406 ERR369247 ERR374163 ERR374164 ERR376915 ERR4020091
ERR4022321 ERR4022377 ERR4022390 ERR4326533 ERR459912 ERR460023 ERR494832
ERR494856 ERR494885 ERR494931 ERR5262857 ERR5384766 ERR5384787 ERR567856
ERR567870 ERR6415236 ERR8266678 ERR852853
SRR10112072 SRR11838929 SRR11839001 SRR11839109 SRR12051556 SRR12051714
SRR12051731 SRR12051739 SRR12051742 SRR12599291 SRR12781297 SRR12781330
SRR13395927 SRR13395966 SRR13395978 SRR13399106 SRR14719138 SRR14719213
SRR14719304 SRR14719315 SRR14863305 SRR315392 SRR350465 SRR6045287
SRR6045468 SRR6045702 SRR6046238 SRR6046274 SRR6046896 SRR6388760
SRR7800440 SRR7800442 SRR7800580 SRR8236425
```
</details>

## Recommended follow-up

Given 79% availability, an **ENA backfill** track is feasible for any
run whose collaborator file never arrives — but it requires re-verifying
content against the manifest and re-calling prophages on the ENA contigs
before it feeds `ntm-v2-extract`. This is intentionally scoped as a
*proposal* here, not executed, to avoid diverging from the collaborator
path without a decision. See task log `ntm-v2-run`.