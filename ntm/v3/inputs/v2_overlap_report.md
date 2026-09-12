# NTM v3 inputs ↔ v2 holdings overlap report

Measured: 2026-09-12 (task `import-ntm-v3`), by direct comparison of the v3
source CSVs against the frozen v2 master manifest
(`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/inputs/NTM_QC_passed_prophage_master_manifest.tsv`,
sha256 `c88686661928c66ab9638da77dd0d615e87bc70975ad859a81e7644142c83386`).
All numbers below were recomputed from the files post-move; the headline counts
match the chat-agent pre-check of 2026-09-12.

## 1. Bridge methodology (contig accessions, unversioned)

Species labels **cannot** be used to join v3 to v2 (v3 uses BV-BRC type-strain
names, v2 uses plain taxonomy — same genome, different strings). The bridge is
therefore **contig-level**:

1. Take every `prophage_contig` value from the v2 master manifest (23,575
   unique, e.g. `NZ_JBFUXV010000022.1` for assemblies, `NODE_…` SPAdes names for
   run assemblies).
2. Take every `scaffold` value from the v3 coordinates CSV per genome.
3. Normalize both sides **unversioned**: strip any leading `accn|` prefix (BV-BRC
   rows) and strip a trailing `.<digits>` version suffix. `NZ_`/`NW_` prefixes
   are kept as-is. SPAdes `NODE_…` contig names carry no version and match
   literally — this is what bridges the run assemblies (the same ENA assembly
   yields identical NODE names in v2 and v3).
4. A v3 prophage-bearing genome is **bridged** if any of its scaffolds matches
   any v2 prophage contig.

Caveat: this bridge only observes prophage-bearing v3 genomes — a genome whose
v3 phigaro calls produced no prophages (or whose v2 record had no prophage)
contributes no contigs to either side, so the 9,382 zero-prophage v3 genomes are
**unbridgeable by this method** (lower bound on true overlap). ID-level matching
(§3) covers the rest.

## 2. Contig-bridge results (prophage-bearing genomes only)

| v3 prophage-bearing genomes | count |
|---|---:|
| bridged to v2 via shared prophage contig | **10,747** |
| **new** (no scaffold seen in v2 prophage contigs) | **6,370** |
| total | 17,117 |

Bridged by v3 `source`: 5,942 run assemblies (`ASSEMBLY`, via matching `NODE_…`
names), 4,675 `NCBI` assemblies, 130 `BV-BRC`. The remaining 9,382 v3 genomes
(zero-prophage) cannot bridge by contigs.

## 3. ID-level overlap (whole cohort, independent of prophage status)

Matching identifiers directly (run accessions literally; `GCA/GCF` by
`<numeric>.<version>` part, i.e. GCA/GCF twins collapse):

| set | count |
|---|---:|
| v2 cohort: 13,122 unique NCBI assemblies (numeric) + 9,543 run assemblies | 22,665 |
| v3 export: 5,451 unique NCBI assemblies + 17,920 run assemblies + 248 BV-BRC `taxid.version` | 26,499 rows / 23,619 unique assemblies |
| v2 runs present in the v3 export | 9,408 of 9,543 |
| v2 numeric assemblies present in the v3 export | 4,901 of 13,122 |
| **v2-only genomes (must be retained from v2): 8,221 numeric + 135 runs** | **8,356** |
| v3-only genomes (new to v3): 550 numeric + 8,512 runs | 9,062 |

Note the correction to the task framing: the v3 export is **not** a pure BV-BRC
snapshot — 17,920 of its 26,499 rows are `source=ASSEMBLY` run assemblies, and
9,408 of v2's 9,543 runs are already inside it. But 135 v2 runs and 8,221 v2
numeric assemblies are **not** in the export, so v2 local holdings must be
retained and unioned in.

## 4. Union rationale (v3 = all-inclusive NTM cohort)

v3 is defined as the **union** of the 2026-09-09 BV-BRC phigaro QC-passed export
and the existing v2 local holdings:

| union component | unique assemblies |
|---|---:|
| NCBI assemblies (numeric, GCA/GCF-deduped): 13,122 (v2) + 550 (v3-only) | 13,672 |
| run assemblies (ERR/SRR/DRR): 9,543 (v2) + 8,512 (v3-only) | 18,055 |
| BV-BRC `taxid.version` genomes (v3 only, separate ID namespace) | 248 |
| **v3 union total** | **31,975** |

Why the union, not the export alone:

1. **8,356 v2 genomes are absent from the export** (8,221 NCBI assemblies +
   135 runs) — dropping v2 would silently shrink the cohort.
2. v2's run-assembly holdings are **already-downloaded local assets**
   (`ntm/v2/run_assemblies/`, 9,543 ENA assemblies) — 9,408 of them are shared
   with the export and must be **linked, not re-downloaded**; the 8,512
   export-only runs are the true download delta on the run side.
3. The 248 BV-BRC `taxid.version` genomes have no NCBI identifier in the export
   at all — acquisition must resolve them through BV-BRC or by contig-accession
   matching (§1), not by assuming an NCBI twin.

## 5. Consequences for `acquire-ntm-v3`

- Acquisition set = export-only identifiers (9,062) + BV-BRC resolvable
  identifiers, minus anything already on NVMe from v2/v1; shared identifiers
  (14,309 by ID: 9,408 runs + 4,901 numerics) should be hard-linked from
  `ntm/v2` holdings where the local file exists.
- The contig bridge (§2) is the prophage-level crosswalk: 10,747 bridged
  genomes' v2 records can be reused for host-clade/annotation context; 6,370
  new prophage-bearing genomes need full v3-side processing.
