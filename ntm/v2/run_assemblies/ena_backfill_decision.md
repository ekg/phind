# ENA backfill feasibility — decision note (2026-08-16)

Status: **feasibility study complete** (input: `ena_probe_results.tsv`, 79 hit
accessions). This note answers the scope question — *can ENA backfill
substitute for collaborator files, or only serve as fallback?*

## Bottom line

**ENA backfill CANNOT substitute for the collaborator files. It can only
serve as a fallback for individual genomes the collaborator cannot deliver —
and only with full re-verification + prophage re-calling before use.**

Evidence (batch = all 79 retrievable probe hits):

| metric | value |
|---|---|
| accessions attempted | 79 (49 prophage-bearing) |
| downloads OK (https ftp.sra.ebi.ac.uk) | 79/79 (100%) |
| content MATCH vs manifest (`assembly_contigs`/`assembly_length_bp`) | **0/79 (0%)** |
| content MISMATCH | 79/79 |
| PanSN bgzip ingested + `samtools faidx` valid | 79/79 |
| quarantined (QUARANTINE.txt marker) | 79/79 |

The probe's 4/4 spot-check finding is reproduced exactly at scale (e.g.
ERR2759527 ENA 48 contigs/5,198,350 bp vs manifest 23/5,193,530 bp;
ERR3177943 ENA 11/5,099,854 vs manifest 9/5,098,843). **No ENA assembly is
byte-identical to the collaborator assembly it would replace.**

## Closeness analysis (how far off are the mismatches?)

The mismatches are *related-but-different*, not unrelated:

- 3/79 have **equal contig count** but different total bp
  (ERR374164 29c/0.14%, SRR350465 60c/0.16%, SRR8236425 267c/0.34%)
- 70/79 (89%) within 1% total length; 35/79 (44%) within 0.1%
- median |total-bp delta| = 0.153%; best near-match ERR374163 (0.0049%,
  4.9 kb on 5.4 Mb) but contig sets differ (78 vs 58 contigs)

Interpretation: ENA submitted assemblies are almost certainly the same
underlying genome but a **different assembly version/refinement** (or a
re-assembly from the same reads). Sequence identity at the contig level is
*plausible* but unproven — it would require per-genome alignment/sketch
comparison (e.g. MASH or minimap2) against the collaborator FASTA, which is
impossible while the collaborator files are still pending. Without the
collaborator file as a reference, "close" is not "equal".

## Why substitution is rejected (not just a policy choice)

1. **Contig tokens differ.** Manifest `prophage_contig/start/end`
   coordinates are exact tokens of the collaborator `contigs.fasta` (e.g.
   `NODE_6_length_418114_cov_19.254028`). ENA headers are ENA-style
   (`SAMEA104138914.contig00001`, `ERZ...`). None of the 79 ENA contig sets
   can address the manifest's prophage coordinates.
2. **Prophage-bearing coverage is the priority cohort.** All 49
   prophage-bearing accessions in the batch mismatch. Using ENA contigs
   without re-calling prophages would silently drop/misplace the 5,907
   prophage-bearing genomes the pipeline depends on.
3. **Content check is the ingestion gate.** The existing
   `ingest_run_assemblies.py` hard-fails on contig/length mismatch — ENA
   assemblies would not pass the collaborator ingestion path by design.

## Fallback path (approved mechanism, per-genome only)

If the collaborator cannot deliver a specific accession, ENA backfill for
that one genome requires, in order:

1. **Download** via `ena_backfill_feasibility.py` (PanSN bgzip + faidx
   already validated — output lives in
   `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/genomes/ena_backfill/`,
   quarantined with `QUARANTINE.txt`).
2. **Re-verify contig set** against any available reference for that genome
   (manifest totals as sanity bound; if the collaborator file later
   arrives, re-ingest and drop the ENA copy).
3. **Re-call prophages** on the ENA contigs (genomad `find_proviruses`)
   for prophage-bearing accessions — manifest coordinates cannot be
   transferred. This is a follow-up task
   (`ntm-v2-ena` → re-call on ENA contigs), **not** part of this
   feasibility pass.
4. **Promote** to the pipeline-visible PanSN dir only after (2)+(3) pass.

## Recommendation to coordinator

- **Keep the collaborator request as the authoritative source.** The
  REQUEST.md upload stands; ENA does not change the 9,543-genome gap.
- **Do not bulk-switch any cohort segment to ENA.** 0/79 content match is
  a hard blocker for substitution.
- **Adopt ENA as per-genome fallback** for accessions the collaborator
  explicitly cannot deliver (e.g. the 21 probe misses have no ENA record at
  all, so those stay collaborator-only either way).
- **Follow-ups:** (a) prophage re-calling on ENA contigs for fallback
  genomes before ntm-v2-extract consumes them; (b) if a nearly-identical
  ENA assembly (delta <0.01%, e.g. ERR374163) is ever needed, sequence
  identity check vs the collaborator file is mandatory before acceptance.

## Artifacts

- Script: `ntm/v2/scripts/ena_backfill_feasibility.py`
- Batch report (79-row per-genome table): `ntm/v2/run_assemblies/ena_backfill_report.md`
- Summary: `ntm/v2/run_assemblies/ena_backfill_summary.json`
- Downloaded + ingested data (quarantined): `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/genomes/ena_backfill/` (113 MB)
