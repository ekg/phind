# S8 — NTM v3 Unified Release Final Verification

Read-only independent recomputation from primary artifacts.
cwd `/home/erikg/phind`; `NVME=/mnt/nvme3n1/erikg/phind-genome-work`.

**Overall verdict: PASS for the release artifacts** (counts, scope reconciliation,
cross-boundary, annotation, sequence integrity all reproduce). **One documentation
drift finding** (§7: `partition_report.md` / `partition_validation.txt` are stale —
data unaffected).

---

## Check table (observed vs expected)

| # | Check | Expected | Observed | Verdict |
|---|---|---|---|---|
| 1 | release ML FASTA records | 1304 | 1304 | PASS |
| 1 | release ancestral FASTA records | 893 | 893 | PASS |
| 1 | ML records missing host_clade_ids / host_scope / species | 0 / 0 / 0 | 0 / 0 / 0 | PASS |
| 1 | ANC records missing host_clade_ids / host_scope / species | 0 / 0 / 0 | 0 / 0 / 0 | PASS |
| 2 | ML sequence byte-identical to source (all 1304; ≥25 sample) | all | 1304/1304, sample 25/25 | PASS |
| 2 | ANC sequence byte-identical to source (all 893; ≥25 sample) | all | 893/893, sample 25/25 | PASS |
| 3 | release_manifest.tsv data rows | 1304 | 1304 | PASS |
| 3 | required cols host_scope / host_clade_ids / n_host_clades / source_breakdown | present | all present | PASS |
| 3 | referenced ml_genome_file exists (base `$NVME/ntm/v3/clades/`) | all | 1304/1304 | PASS |
| 3 | referenced ancestral_genome_file exists (ancestral=yes rows) | all 893 | 893/893 | PASS |
| 3 | host_scope values ⊆ {NTM,MTC} for genomes | yes | NTM 1285, MTC 12, **MIXED 7** | PASS w/ nuance (see §3) |
| 4 | MTC members | 391 | 391 | PASS |
| 4 | wholly-MTC prophage clades | 12 | 12 | PASS |
| 4 | MIXED clades | 7 | 7 | PASS |
| 4 | fully-NTM clades | 1285 | 1285 | PASS |
| 4 | clades touched (any MTC member) | 19 | 19 | PASS |
| 4 | host_clade_scope.tsv vs recomputed host-clade majority | 0 mismatch | 0 mismatch (443 clades, 1 MTC) | PASS |
| 4 | per_scope_summary.tsv (all 15 numbers) | as published | all 15 match | PASS |
| 5 | cross_boundary_report.tsv rows | 7 | 7 | PASS |
| 5 | minority scope + member counts (per clade) | match recompute | 7/7 match | PASS |
| 5 | confidence column consistent with documented rule | consistent | 7/7 consistent | PASS |
| 6 | per-genome GFF count | 2197 | 2197 | PASS |
| 6 | annotation_gff_index.tsv rows | 2197 | 2197 | PASS |
| 6 | sum of gene_count over index | report | **272,974** | INFO |
| 6 | every index genome_id ∈ ML∪ancestral set | all | 2197/2197 (0 missing) | PASS |
| 6 | index gene_count == actual CDS count in GFF | all | 2197/2197 | PASS |
| 6 | 3 spot-checked GFFs: `##gff-version 3` + coords ≤ genome length | yes | 3/3 valid | PASS |
| 7 | release counts vs ml_report.md | match | match (1304/893/411) | PASS |
| 7 | release counts vs mash_clades_report.md | match | match (36857/1304/893/411) | PASS |
| 7 | release counts vs partition_report.md | match | **813/456/357, 9,446 prophages** | **FAIL (stale doc)** |
| 7 | release counts vs s3_clade_validation.md | match | match (391/12/7/1285) | PASS |

---

## Scope reconciliation table (independent)

Rule applied: per-genome MTC iff `host_clades.tsv.species` `startswith` one of
`Mycobacterium tuberculosis | canetti | orygis | africanum | bovis | microti`
(the binding prefix set from the documented rule; genus included, per §"MTC_PREFIXES").
Host clade = MTC iff strictly >50% of its genomes are MTC. Members inherit host-clade scope.

| metric | independently recomputed | artifact claim | match |
|---|---:|---:|---|
| host clades | 443 | 443 (`host_clade_scope.tsv`) | ✅ |
| MTC host clades | 1 (`host_clade_0002`, 7273/7467 = 0.974) | 1 | ✅ |
| MTC members | 391 | 391 | ✅ |
| wholly-MTC prophage clades | 12 | 12 | ✅ |
| MIXED clades | 7 | 7 | ✅ |
| fully-NTM clades | 1285 | 1285 | ✅ |
| clades touched by MTC | 19 (=12+7) | 19 | ✅ |
| host-clade scope mismatches vs `host_clade_scope.tsv` | 0 | — | ✅ |

`per_scope_summary.tsv` independent recomputation (all match):

| host_scope | n_clades | n_ml | n_anc | n_members | n_host_clades |
|---|---:|---:|---:|---:|---:|
| NTM | 1285 | 1285 | 878 | 36228 | 289 |
| MTC | 12 | 12 | 8 | 202 | 1 |
| MIXED | 7 | 7 | 7 | 427 | 26 |

Statement cross-checks also reproduce: clade-derived MTC members = 9
(`0_0306`:1, `0_0332`:5, `0_0336`:2, `0_0338`:1); species-prefix-only alternate
accounting = 382 MTC members / 10 wholly-MTC / 9 MIXED / 1285 fully-NTM.

## Sample-integrity results (byte-identical sequences)

Compared release FASTA sequences to source `$NVME/ntm/v3/ml/*` by record id
(headers differ only in the appended `host_clade_ids/host_scope/species` fields).

| set | source records | release records | ids only-in-release | ids only-in-source | all-sequence mismatches | random 25-sample |
|---|---:|---:|---:|---:|---:|---:|
| ML | 1304 | 1304 | 0 | 0 | 0 | 25/25 identical |
| ANC | 893 | 893 | 0 | 0 | 0 | 25/25 identical |

Cross-boundary (7/7 match):
`0_0306` NTM×7/93, `0_0338` NTM×10/90, `0_0365` MTC×2/51, `0_0252` MTC×1/73,
`0_0413` MTC×1/78, `0_0435` MTC×1/17, `0_0657` MTC×1/3.
Confidence rule `minority≤1→low; median≤0.10 & n≥10→high; median≤0.25 & n≥5→medium`
reproduces all 7 labels (2 high, 5 low).

---

## §3 nuance — `host_scope` MIXED in manifest/headers

`release_manifest.tsv` `host_scope` is the **per-prophage-clade** label and contains
`MIXED` for the 7 cross-boundary clades (`release_scope_statement.md`: per-genome ∈
{NTM,MTC}; per-clade ∈ {NTM,MTC,MIXED}). The 7 release FASTA headers likewise carry
`host_scope=MIXED`. The **per-genome binary** scope (NTM/MTC only) is in
`release/host_clade_scope.tsv` — 442 NTM + 1 MTC host clades, no MIXED. A literal
reading of check 3 ("only in {NTM,MTC} for genomes") would flag 7 rows, but the
manifest column is documented per-clade; treat as PASS-under-documented-semantics.

Minor naming note: `release_scope_statement.md` calls `host_clade_scope.tsv`
"Per-genome (host) scope", but the file is keyed per **host clade** (443 rows),
not per release genome.

---

## §7 drift finding — stale partition report

`ntm/v3/partition_report.md` (dated 2026-09-14) and `ntm/v3/partition_validation.txt`
(executed 2026-09-14) state a prophage/clade universe that predates the final
2026-10-06 expansion:

| quantity | partition_report.md | release / final | note |
|---|---:|---:|---|
| tight clades | 813 | 1304 | disagrees |
| alignable | 456 | 893 | disagrees |
| singletons | 357 | 411 | disagrees |
| prophages | 9,446 | 36,857 | disagrees |
| partition_summary.tsv rows | 456 (doc) | 893 (actual NVMe file) | actual file regenerated |

The actual NVMe artifacts agree with the release: `clades/` has 1304 clade dirs,
893 `partitions.bed`, 893 `ml.ml.fa`, 893 `anc.ancestral.genome.fa`,
`partition_summary.tsv` = 893 rows; repo `clades/alignable_clades.tsv` = 893,
`singletons.tsv` = 411, `clade_summary.tsv` = 1304. So this is **documentation
drift only** — the partition step was re-run for the expanded set but the two
reports were not regenerated. `ml_report.md` and `mash_clades_report.md` both
match the release. Upstream note: `mash_clades_report.md` self-reports a
length-band mismatch vs `extract_report` (263 vs 91 outside [1k,100k]); not a
release-level claim.

---

## Exact recomputation commands

```bash
NVME=/mnt/nvme3n1/erikg/phind-genome-work

# (1) counts + header fields
grep -c '^>' $NVME/ntm/v3/release/all_ntm_v3_ml_phage_genomes.fa          # 1304
grep -c '^>' $NVME/ntm/v3/release/all_ntm_v3_ancestral_phage_genomes.fa   # 893
python3 - <<'PY'   # verify each header has host_clade_ids=/host_scope=/species=
import re
for p in ["release/all_ntm_v3_ml_phage_genomes.fa","release/all_ntm_v3_ancestral_phage_genomes.fa"]:
    miss=0
    for l in open(f"/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/{p}"):
        if l.startswith(">") and not all(re.search(rf"(^|\s){k}=",l) for k in ("host_clade_ids","host_scope","species")):
            miss+=1
    print(p,"missing header fields:",miss)
PY

# (2) byte-identical sequence check (parse both FASTAs, compare by record id) -> /tmp/s8/integrity.py
#   result: ML 1304/1304, ANC 893/893, 0 mismatches; sample 25/25 each

# (3) manifest
tail -n +2 $NVME/ntm/v3/release/release_manifest.tsv | wc -l        # 1304
head -1  $NVME/ntm/v3/release/release_manifest.tsv | tr '\t' '\n'
# file refs (base = $NVME/ntm/v3/clades/): all 1304 ml + 893 ancestral exist

# (4) scope reconciliation (independent) -> /tmp/s8/scope.py / scope2.py
#   391 / 12 / 7 / 1285 / 19; 0 host-clade-scope mismatches; per_scope_summary 15/15

# (5) cross-boundary -> /tmp/s8/cross.py   (7/7 match, confidence 7/7)

# (6) annotation
grep -c '^genome_id' /home/erikg/phind/ntm/v3/annotation_gff_index.tsv   # 2197
ls $NVME/ntm/v3/annotation/pharokka_out/per_genome_gff/*.gff | wc -l     # 2197
awk -F'\t' 'NR>1{s+=$3}END{print s}' /home/erikg/phind/ntm/v3/annotation_gff_index.tsv  # 272974
# every index genome_id ∈ ML∪anc; gene_count == CDS count in GFF (2197/2197)
# 3 GFF spot-checks valid ##gff-version 3, max coord <= genome length
```

---

## Residual uncertainty / risks

- **Stale partition docs** (finding §7): `partition_report.md` + `partition_validation.txt`
  misstate the final clade universe; regenerate them or add a superseded notice.
- **host_scope semantics**: manifest/header `host_scope` is per-clade (MIXED allowed).
  Any downstream consumer expecting strict {NTM,MTC} per row must use per-genome scope
  (`host_clade_scope.tsv`, keyed per host clade).
- **`host_clade_scope.tsv` naming**: statement calls it per-genome but it is per-host-clade.
- **Ancestral headers** carry `status=ml` (inherited from the source FASTA); cosmetic only.
- **Scope classifier sensitivity**: reconciliation reproduces exactly only with the full
  genus-included prefix set (`Mycobacterium ...`); a bare-epithet interpretation yields
  0 MTC and must not be used.
- Upstream `mash_clades_report.md` length-band mismatch (263 vs 91) is outside the release
  artifact set and was not adjudicated here.

```acceptance-report
{
  "criteriaSatisfied": [
    {
      "id": "criterion-1",
      "status": "satisfied",
      "evidence": "All 7 checks recomputed from primary artifacts; release artifacts PASS. Check table + scope reconciliation table + sample integrity included. One documentation drift finding reported (stale partition_report.md/partition_validation.txt); residual risks listed."
    }
  ],
  "changedFiles": [],
  "testsAddedOrUpdated": [],
  "commandsRun": [
    {
      "command": "grep -c '^>' release ML/ancestral FASTA + header-field parse",
      "result": "passed",
      "summary": "ML=1304, ANC=893; 0 records missing host_clade_ids/host_scope/species"
    },
    {
      "command": "python3 integrity compare release vs source FASTAs by id",
      "result": "passed",
      "summary": "ML 1304/1304 and ANC 893/893 byte-identical; 25/25 sample each identical"
    },
    {
      "command": "python3 manifest checks (rows, cols, file existence, scope values)",
      "result": "passed",
      "summary": "1304 rows; required cols present; all ml/ancestral file refs exist; NTM 1285/MTC 12/MIXED 7"
    },
    {
      "command": "python3 independent scope reconciliation from host_clades.tsv",
      "result": "passed",
      "summary": "391/12/7/1285/19 match expected; 0 host-clade-scope mismatches; per_scope_summary 15/15"
    },
    {
      "command": "python3 cross_boundary recompute + confidence rule",
      "result": "passed",
      "summary": "7/7 rows match minority scope/counts; confidence 7/7 consistent"
    },
    {
      "command": "python3 annotation checks (count, index, gene_count, gff validity)",
      "result": "passed",
      "summary": "2197 GFFs, 2197 index rows, gene_count sum 272974, 2197/2197 ids in ML+anc, 3/3 spot-checks valid"
    }
  ],
  "validationOutput": [
    "Release ML=1304, ancestral=893; headers complete.",
    "Sequence integrity: 100% byte-identical to source (all records, not just sample).",
    "Scope reconciliation reproduces 391/12/7/1285/19 and per_scope_summary exactly.",
    "Cross-boundary 7/7 match incl. confidence rule.",
    "Annotation: 2197 GFFs, index 2197 rows, gene_count sum 272974, all ids valid, GFF format valid.",
    "Drift: ml_report and mash_clades_report match; partition_report.md stale (813/456/357, 9446 prophages)."
  ],
  "residualRisks": [
    "partition_report.md and partition_validation.txt are stale (pre-2026-10-06 clade set of 813/456/357, 9446 prophages) and could mislead downstream; actual NVMe partition artifacts agree with release (893 alignable, 1304 clades).",
    "release_manifest.tsv / release FASTA header host_scope uses MIXED for 7 cross-boundary clades (per-clade label); per-genome binary scope is in host_clade_scope.tsv (per host clade).",
    "release_scope_statement.md calls host_clade_scope.tsv 'Per-genome (host) scope' though it is keyed per host clade.",
    "Ancestral FASTA headers carry status=ml (inherited from source), cosmetic.",
    "Scope classifier must use the full genus-included prefix set; a bare-epithet reading yields 0 MTC.",
    "Upstream mash_clades_report.md self-reported length-band mismatch (263 vs 91) not adjudicated (outside release artifacts)."
  ],
  "noStagedFiles": true,
  "diffSummary": "No repository files modified; read-only verification. Only the routed output report was written.",
  "reviewFindings": [
    "blocker: none for release artifacts",
    "finding: ntm/v3/partition_report.md: 'Total tight clades: 813' / alignable 456 / singletons 357 / 9,446 prophages disagrees with final release 1304/893/411/36,857 and with regenerated NVMe artifacts; documentation drift only",
    "finding: ntm/v3/partition_validation.txt: 'clades total=813 alignable=456 singletons=357' same stale set",
    "note: release_manifest.tsv host_scope has 7 MIXED rows (documented per-clade semantics)"
  ],
  "manualNotes": "All release artifacts verified PASS by independent recomputation. The only issue is stale partition documentation; the underlying partition outputs on NVMe were re-run for the 1304-clade set (893 partitions.bed / ml.ml.fa / anc.ancestral.genome.fa, partition_summary.tsv 893 rows) and are consistent with the release. Reports used only as claims; all counts recomputed from primary artifacts."
}
```
