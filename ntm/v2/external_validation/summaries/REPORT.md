# Public mycobacteriophage reference panel — build report

Panel size: **6266** retained records (3065 mycobacteriophages / 3200 other actinobacteriophages from the PhagesDB leg).

## Source coverage and reconciliation

| source | metric | count |
|---|---|---:|
| PhagesDB | metadata rows (sequenced phages) | 6018 |
| PhagesDB | bulk FASTA records | 6026 |
| PhagesDB | retained entries | 6009 |
| INPHARED | table rows (release 7Apr2026) | 36025 |
| INPHARED | in-scope mycobacteriophage rows | 3309 |
| INPHARED | out-of-scope rows (non-myco hosts; by design) | 32716 |
| INPHARED | in-scope with sequence | 3308 |
| NCBI | esearch uids | 3857 |
| NCBI | backfill accessions resolved by id (esearch index lag) | 2824 |
| NCBI | records with FASTA | 6674 |
| NCBI | NCBI-only records retained after screening | 42 |
| merge | accession-union entries | 6988 |
| merge | exact-sequence duplicates merged | 722 |
| merge | retained panel records | 6266 |

Exclusions with explicit reasons: 554 (see `exclusions.tsv`).

NCBI versioned accessions used as primary identifiers for **5542** records; **5544** of all primary accessions carry an explicit `.N` version. Records whose primary accession is a PhagesDB-only phage (no GenBank deposit) keep the PhagesDB name as identifier.

## Source-combination distribution

| sources | records |
|---|---:|
| `ncbi;phagesdb` | 2893 |
| `ncbi;inphared;phagesdb` | 2392 |
| `phagesdb` | 724 |
| `ncbi;inphared` | 249 |
| `ncbi` | 8 |

## Host-label distribution (top genera)

| host genus (reported) | records |
|---|---:|
| Mycobacterium | 3037 |
| Gordonia | 926 |
| Microbacterium | 840 |
| Arthrobacter | 767 |
| Streptomyces | 435 |
| Rhodococcus | 78 |
| Propionibacterium | 57 |
| Curtobacterium | 54 |
| Corynebacterium | 37 |
| Mycolicibacterium | 25 |
| Brevibacterium | 3 |
| (unresolved) | 2 |
| Tsukamurella | 2 |
| Mycobacteroides | 2 |
| Tetrasphaera | 1 |

Host labels resolved: 6264; unresolved: 2 (see completeness caveats).

## Evidence class (cultured vs MAG/provirus)

| evidence_class | records |
|---|---:|
| isolated_sequenced | 6009 |
| deposited_isolate_sequence | 257 |

`isolated_sequenced` = PhagesDB plaque-purified, sequenced isolate. `deposited_isolate_sequence` = GenBank/INPHARED deposit without culture status in machine-readable form. `predicted_prophage` / `metagenome_assembled` are keyword-derived from source descriptions (heuristic — see caveats).

## Completeness / topology

| completeness label | records |
|---|---:|
| finished_sequence | 5948 |
| complete_genome_record | 289 |
| (none) | 20 |
| length_reported_only | 9 |

## Caveats

1. **Reported host ≠ verified host range.** Every host field is a *reported label* (`host_evidence` column records its provenance: PhagesDB metadata, INPHARED table, or a GenBank `/host=` qualifier). None of these represent experimentally verified host range.
2. **RefSeq/GenBank dual membership.** Many phages carry both a GenBank and a RefSeq accession for the same sequence; these are merged as exact duplicates and both accessions remain in the crosswalk (`relation=same_sequence`).
3. **PhagesDB terminal overhangs.** PhagesDB bulk FASTA includes terminal overhang / extended ends; the GenBank (NCBI) version of the same phage can be shorter. Where both exist and bytes differ, the NCBI record is the canonical carrier and the PhagesDB variant is recorded in the crosswalk (`relation=sequence_variant`) with its own sha256.
4. **Exact-duplicate dedup only.** Dedup is by accession alias and by exact (uppercase) sequence hash; it is not reverse-complement or near-duplicate aware.
5. **Evidence-class heuristics.** MAG/provirus classification is keyword-derived from source descriptions where PhagesDB culture metadata is absent; treat as provisional.
6. **INPHARED scope.** Only Mycobacteriaceae-host (or mycobacteriophage-named) INPHARED rows are included; other actino hosts come from the PhagesDB leg only (task scope).
7. **Completeness labels are heterogeneous** (PhagesDB 'finished', NCBI title 'complete genome', INPHARED length-only) and are not comparable across sources without the `completeness_source` column.

## Licensing / redistribution

- **PhagesDB**: public research database; bulk downloads provided for research use (no explicit dataset license statement; site terms at phagesdb.org). Redistribute the cached raw FASTA with attribution.
- **INPHARED (Millard Lab, 7Apr2026 release)**: aggregates GenBank/RefSeq/ENA records plus PhagesDB; redistribution constraints follow the underlying repositories (mostly public-domain/CC0 for RefSeq; GenBank submitters retain rights but records are published openly).
- **NCBI GenBank/RefSeq**: U.S. Government public domain for RefSeq; GenBank records are openly available; submitters retain copyright on submitted sequences.
- The **panel FASTA is stored externally** at `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/external_validation/panel/`; the repository commits code, manifests, summaries and provenance only.

## Determinism

All outputs are pure functions of `cache/` (see PROVENANCE.md). Rebuilds from the same cache are byte-identical; `verify_panel.py --determinism` enforces this.
