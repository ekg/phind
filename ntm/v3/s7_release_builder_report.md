# S7 — NTM v3 unified release builder (revised: clade-level scope)

**Status:** implemented + validated (dry-run on full set + 5-clade subset write). No full release run.

**Revision:** applied the orchestrator's decision to derive `host_scope` at the **host-clade** level
(majority MTC vote) instead of per genome by species prefix alone. All expected counts now match.

## Deliverable

- Script: `ntm/v3/scripts/build_v3_release.py` (new; read `ntm/v2/scripts/build_v2_release.py` and
  `ntm/scripts/build_ntm_release.py` as prior art).
- Change on top of the previous revision: `load_host_clades()` now also derives a per-host-clade scope
  by majority species-prefix vote; every prophage member inherits its host clade's scope. Added
  `host_clade_scope.tsv` and expanded `release_scope_statement.md`.

## Full invocation (for the parent)

```bash
python3 /home/erikg/phind/ntm/v3/scripts/build_v3_release.py \
  --ml-dir       /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/ml \
  --host-clades  /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/host_clades/host_clades.tsv \
  --clade-summary /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/clades/clade_summary.tsv \
  --out          /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/release
```
All flags default to those paths, so a bare `python3 ntm/v3/scripts/build_v3_release.py` is equivalent.
`--dry-run` computes/reports without writing. Debug subset: `--clade-ids a,b,c` or `--limit N`.
Reads only `ntm/v3/ml` + `ntm/v3/host_clades`; modifies no upstream artifact.

### Outputs written to `$NVME/ntm/v3/release/`

| file | content |
|---|---|
| `all_ntm_v3_ml_phage_genomes.fa` | 1304 ML genomes; header carries `host_clade_ids`, `host_scope`, `species` |
| `all_ntm_v3_ancestral_phage_genomes.fa` | 893 ancestral genomes; same header convention |
| `release_manifest.tsv` | input manifest + `host_scope`, `host_clade_ids`, `n_host_clades`, `source_breakdown` |
| `per_scope_summary.tsv` | NTM / MTC / MIXED roll-up |
| `cross_boundary_report.tsv` | MIXED clades; minority scope, `median_mash`, `confidence` |
| `host_clade_scope.tsv` | per host clade: derived scope, genome counts, MTC fraction (audit of the derivation) |
| `release_scope_statement.md` | release-note scoping statement (rules + 9 clade-derived assignments + both accountings) |

## Scope rule (as implemented)

1. **Per-genome MTC classifier (kept):** `host_clades.tsv.species` startswith
   `Mycobacterium tuberculosis|canetti|orygis|africanum|bovis|microti` → MTC, else NTM
   (`startswith` excludes `Mycobacterium avium subsp. paratuberculosis`).
2. **Per-host-clade scope (new):** a host clade is MTC when **more than half** of its genomes carry an
   MTC prefix, else NTM. Only `host_clade_0002` qualifies (7273/7467 = 97.4% MTC; all other 442 host
   clades NTM).
3. Every prophage member inherits its host clade's scope. **Per-genome `host_scope` ∈ {NTM, MTC}** (binary).
4. **Per-prophage-clade label ∈ {NTM, MTC, MIXED}**: MTC/NTM when all members inherit that scope,
   MIXED when both occur. Label only — no filtering.
5. Join unchanged: member `canonical_acc#1#prophage_id` → `host_clades.tsv` via `canonical_acc`; all
   36,857/36,857 members join, 0 unjoined.
6. Cross-boundary `confidence`: `low` if only 1 minority member; else `high` if `median_mash ≤ 0.10 and
   n_members ≥ 10`; else `medium` if `median_mash ≤ 0.25 and n_members ≥ 5`; else `low`.

## Dry-run output (full set, verbatim)

```
clades=1304 ml_genomes=1304 ancestral_genomes=893
members by member-scope (clade-derived): NTM=36466 MTC=391
  [species-prefix-only comparison] clades: NTM=1285 MTC=10 MIXED=9
  host_scope=NTM: clades=1285 members=36228 host_clades=289
  host_scope=MTC: clades=12 members=202 host_clades=1
  host_scope=MIXED: clades=7 members=427 host_clades=26
  host_scope=UNJOINED: clades=0 members=0 host_clades=0
cross-boundary (MIXED) clades: 7
   0_0306 n=100 NTM=7 MTC=93 minority=NTMx7 median_mash=0.035272 confidence=high
   0_0338 n=100 NTM=10 MTC=90 minority=NTMx10 median_mash=0.045108 confidence=high
   0_0365 n=53 NTM=51 MTC=2 minority=MTCx2 median_mash=0.026023 confidence=high
   0_0252 n=74 NTM=73 MTC=1 minority=MTCx1 median_mash=0.200342 confidence=low
   0_0413 n=79 NTM=78 MTC=1 minority=MTCx1 median_mash=0.021322 confidence=low
   0_0435 n=18 NTM=17 MTC=1 minority=MTCx1 median_mash=0.113401 confidence=low
   0_0657 n=3 NTM=2 MTC=1 minority=MTCx1 median_mash=0.239237 confidence=low
dry-run: no files written
```

## New counts vs expectation, and diff vs the previous revision

| quantity | prev (species-prefix) | **new (clade-level)** | expected | match |
|---|---:|---:|---:|---|
| MTC prophage members | 382 | **391** | 391 | ✅ |
| wholly-MTC clades | 10 | **12** | 12 | ✅ |
| MIXED clades | 9 | **7** | 7 | ✅ |
| fully-NTM clades | 1,285 | **1,285** | 1,285 | ✅ |
| clades touched by ≥1 MTC member | 19 | **19** | 19 | ✅ |
| MTC host genomes | 7,273 | 7,273 | 7,273 | ✅ |
| clades / ML / ancestral genomes | 1304/1304/893 | 1304/1304/893 | — | ✅ |

Diff = 9 members reclassified NTM→MTC; wholly-MTC +2 (`0_0332`, `0_0336`); MIXED −2; NTM clades unchanged.
Independent recompute (standalone script over `host_clades.tsv` + `release_manifest.tsv`) reproduced
391 / 12 / 7 / 1285 / 19 exactly.

Clade-derived MTC assignments (9 members; the `Mycobacterium sp.` hosts inside `host_clade_0002`):
`0_0306`: 1, `0_0332`: 5, `0_0336`: 2, `0_0338`: 1.

`host_clade_scope.tsv` MTC row: `host_clade_0002  MTC  7467  7273  0.9740`.

## Subset test evidence (`--clade-ids 0_0000,0_0306,0_0332,0_0252,0_0657 --out /tmp/ntm_v3_rel_test`)

- Wrote all 7 artifacts. Subset result: `0_0332` now `host_scope=MTC` (was MIXED); `0_0306` MIXED
  with MTC=93 (was 92); subset clade-derived MTC = 6 (1+5).
- `release_manifest.tsv`: input columns 1–11 + new 12–15; e.g.
  `0_0332 → MTC, host_clade_0002, n_host_clades=1, NTM=0;MTC=100`.
- Header example: `>ntm3_0_0332_ML status=ml n_members=100 length=149730 host_clade_ids=host_clade_0002 host_scope=MTC species=Mycobacterium sp.,Mycobacterium tuberculosis`.
- `release_scope_statement.md` renders both accountings + the clade-derived table + namespaces.
- Sequences remain byte-identical to source (verified earlier by assert; header rewrite only).

## Residual risks / open items

1. **MIXED third value for the per-clade label.** PLAYBOOK says `host_scope ∈ {NTM, MTC}`; the decision
   confirms per-genome is binary and the **per-prophage-clade** label may be `MIXED`. Implemented accordingly.
2. **`Mycobacterium sp.` label resolution** now handled by clade majority; if the parent later resolves those
   9 genomes' species the classification is unchanged (they remain MTC via `host_clade_0002`).
3. Mixed clades `0_0306`/`0_0336`/`0_0338` carry saturated pairwise MASH (max=1.0); calls rest on the clade
   median (median-tight rule), same caveat as `s3_clade_validation.md`.
4. No full run performed (parent will); loading the 75 MB + 70 MB FASTAs is the only full-run cost (~seconds).

## Acceptance

- one changed file, `ntm/v3/scripts/build_v3_release.py`; no other source touched; no staged files.
- evidence: compile OK, full-set dry-run with matching counts, 5-clade subset write, independent recompute.
