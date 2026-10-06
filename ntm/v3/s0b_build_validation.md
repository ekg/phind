# S0b Build Verification — NTM v3 canonical PanSN objects

**Scope:** read-only verification of `BUILD_DONE built=17920 reused=0 failed=0`.
**NVME:** `/mnt/nvme3n1/erikg/phind-genome-work`
**Objects:** `$NVME/ntm/v3/genomes/canonical_objects/<run>/<run>.pansn.fa.gz{,.fai,.gzi}`
**Verdict: PASS** — all six checks passed; independent recomputation matches the driver claim. Two minor factual deviations from the task prompt are noted (not defects, see Notes).

## Checks summary

| # | Check | Expected | Observed | Status |
|---|-------|----------|----------|--------|
| 1 | Run objects present, 3 files non-empty | 17,920 dirs, 0 missing/empty | 17,920 dirs (ERR 8370 / SRR 8679 / DRR 871); missing=0, empty=0 | **PASS** |
| 2 | `.fai` first col PanSN `<run>#1#<contig>` | all sampled match | 24 runs (8/namespace), 0 deviations | **PASS** |
| 3 | `.fai` line count == master `n_contigs` | 0 mismatches | 17,920/17,920 compared; **0 mismatches** | **PASS** |
| 4 | `bgzip -t` integrity | 0 failures | 24 runs, 0 failures | **PASS** |
| 5 | Coordinate scaffolds resolve in `.fai` bare names | 0 unresolved | 24 runs (8/namespace), 0 unresolved | **PASS** |
| 6 | Non-run objects unchanged | 34,163 total / 16,243 non-run | exact match; pre-existing GCA_/GCF_ `.fai` valid; mtime 2026-09-12 | **PASS** |

## Check 1 — object inventory
```
total dirs:        34163
run dirs:          17920   (regex ^(ERR|SRR|DRR)[0-9]+$)
non-run dirs:      16243
```
Per-file scan over all 17,920 run dirs for `<run>.pansn.fa.gz`, `.pansn.fa.gz.fai`, `.pansn.fa.gz.gzi`:
`missing=0 empty=0`. Namespace split: ERR 8370, SRR 8679, DRR 871 (sum 17920).

## Check 2 — PanSN naming (sample, 24 runs)
Sampled 8 ERR + 8 SRR + 8 DRR. Every line in each sampled `.fai` began with `<run>#1#`; 0 deviations.
Example (`ERR5412094`):
```
ERR5412094#1#NODE_1_length_309207_cov_21.060488	309207	49	60	61
ERR5412094#1#NODE_2_length_286865_cov_20.251410	286865	314459	60	61
```

## Check 3 — exhaustive contig-count comparison (strongest check)
For all 17,920 runs, `.fai` line count (`wc -l`) vs master `NTM_master_table_complete_20260922.tsv` `n_contigs` (field 28, `source==ASSEMBLY`):
```
joined rows:            17920
mismatch count:         0
run missing in master:  0
master row unused:      0
```
Set equality also holds: master ASSEMBLY `genome_id` set == run-dir set (0-only-in-dirs, 0-only-in-master).
Aggregate contig totals agree: master/of-summary 10,311,626 == observed 10,311,626.

## Check 4 — bgzip integrity (sample, 24 runs)
`bgzip -t` on 24 objects (8 per namespace): failures=0.

## Check 5 — scaffold resolution (sample, 24 runs)
For 24 runs present in `ntm/v3/inputs/ntm_qc_passed_phigaro_coordinates_20260909.csv`, every coordinate `scaffold` (col 5) was matched against the object's `.fai` bare names (`cut -f1 | awk -F'#' '{print $NF}'`): unresolved=0.
(There are 11,323 distinct run genome_ids with coordinates; 24 sampled.)

## Check 6 — pre-existing non-run objects untouched
- Total 34,163, non-run 16,243, run 17,920 — exactly as expected.
- 8 sampled `GCA_`/`GCF_` objects all have valid non-empty `.fai` (e.g. `GCA_001207945.1` 271 contigs, `GCA_000172115.1` 353).
- Directory mtimes: sampled non-run dir `2026-09-12 22:17Z` vs run dir `2026-10-06 16:44Z` — non-run predates the build.

## Exact commands (reproducible)
```bash
NVME=/mnt/nvme3n1/erikg/phind-genome-work
OBJ=$NVME/ntm/v3/genomes/canonical_objects
ls $OBJ | wc -l
ls $OBJ | grep -Ec '^(ERR|SRR|DRR)[0-9]+$'
ls $OBJ | grep -Evc '^(ERR|SRR|DRR)[0-9]+$'

# check 1: all run dirs, per-file presence/size
while read r; do for f in "$r.pansn.fa.gz" "$r.pansn.fa.gz.fai" "$r.pansn.fa.gz.gzi"; do
  p="$OBJ/$r/$f"; [ -e "$p" ] || echo "MISSING $r $f"; [ -s "$p" ] || echo "EMPTY $r $f"; done; done

# check 3: exhaustive
awk -F'\t' 'NR>1 && $2=="ASSEMBLY"{print $1"\t"$28}' NTM_Data/NTM_master_table_complete_20260922.tsv | sort > m.s
while read r; do echo -e "$r\t$(wc -l < "$OBJ/$r/$r.pansn.fa.gz.fai")"; done < rundirs.txt | sort > o.s
join -t$'\t' m.s o.s | awk -F'\t' '$2!=$3' | wc -l   # -> 0

# check 4
bgzip -t "$OBJ/<run>/<run>.pansn.fa.gz"

# check 5
awk -F',' -v R=<run> 'NR>1 && $1==R{print $5}' coordinates.csv | sort -u > sc.txt
cut -f1 "$OBJ/<run>/<run>.pansn.fa.gz.fai" | awk -F'#' '{print $NF}' | sort -u > bare.txt
comm -23 sc.txt bare.txt
```

## Notes / prompt-vs-reality
- `.fai`/`.gzi` filenames are `<run>.pansn.fa.gz.fai` / `.gzi` (suffix appended to the full `.fa.gz`), not `<run>.pansn.fai`. Paths in the prompt abbreviated them; actual naming above.
- The prompt said "17,920 delivered SRA-run genomes" but the objects cover ERR+SRR+DRR namespaces; all 17,920 are accounted for (ERR 8370, SRR 8679, DRR 871).

## Residual risks / uncertainties
- bgzip integrity and PanSN naming were sampled (24 runs), not exhaustive; the exhaustive contig-count check bounds but does not fully prove per-record content integrity. A full `bgzip -t` sweep of all 17,920 would close this.
- `.fai` contig counts were compared to the master table; the master `n_contigs` itself was not independently re-derived from raw assemblies.
- Coordinate scaffold resolution was verified for 24 of 11,323 run genome_ids present in the coordinates file. No evidence of a systemic gap, but a full sweep is not done here.
- Non-run "unchanged" is supported by count + mtime + spot-check, not by content hashing.

## Overall verdict
**PASS.** The driver's `BUILD_DONE built=17920 reused=0 failed=0` is independently confirmed: complete object inventory, PanSN naming, exhaustive per-run contig-count agreement (0 mismatches), bgzip integrity, scaffold resolution, and untouched pre-existing objects. Stop/ask condition (a discrepancy) was not triggered.
