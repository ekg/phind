# AMENDMENTS — local-ntm-logan run

Amendments are logged per RUN_PLAN §0. None of these change thresholds,
budgets, bait panel, or the hit definition.

## A1 — 2026-08-19T14:2xZ — EUtils compound-query defect for the `prophage` term

**Observed (receipts in `logs/eutils.log` + `logs/eutils_calibration.log`):**

| query | 14:09Z calibration | 14:2xZ fetch |
|---|---|---|
| `"Mycobacterium smegmatis"[Organism] AND prophage[All Fields]` (preregistered A2) | 1,590 | **0** (XML and JSON, 5 repeats) |
| `"Mycobacterium smegmatis"[Organism] AND prophage[Title]` | — | 0 |
| `"Mycobacterium smegmatis"[Organism] AND prophages[All Fields]` | — | 0 |
| `(... AND prophage[All Fields]) OR prophages[All Fields]` | — | 931 |
| `prophage[All Fields]` (bare) | — | 891 (stable) |
| `"Mycobacterium smegmatis"[Organism]` | 2,541 | 2,541 (stable) |
| `... AND mitomycin[All Fields]` (A1) | 46 | 46 (stable) |

NCBI EUtils currently returns 0 for any compound query intersecting
`prophage[All Fields]` with an organism clause, while the bare term and the
organism term each resolve normally — a service-side posting-list defect,
not an empty result set (the A2 calibration run 20 min earlier returned
1,590 and the bare term currently yields 891 SRA entries mentioning
prophage).

**Remediation (semantics-preserving, no scope change):** execute A2 as the
bare term `prophage[All Fields]` (uid list + sha256 recorded exactly like
the preregistered queries, id `A2bare`), fetch run metadata via
`efetch runinfo`, and apply the preregistered A2 organism restriction
(*Mycobacterium smegmatis*, scientific-name substring match) **locally**.
This reproduces precisely the intended intersection
`"Mycobacterium smegmatis"[Organism] AND prophage[All Fields]`; the only
difference is where the AND is evaluated. `tierA_manifest.json` records
`A2` with `"status": "amended_see_AMENDMENTS_A1"`.

**A4 NTM filter (preregistered in RUN_PLAN §5 as "post-filter to NTM"):**
scientific_name genus in {Mycobacterium, Mycolicibacterium, Mycobacteroides,
Mycolicibacillus} (modern NTM genera after the sensu-lato split — EUtils
returns e.g. `Mycolicibacterium smegmatis` for renamed runs) and does not
contain `tuberculos|bovis|africanum|canettii|leprae|caprae|pinnipedii|mungi|orycis`.
Applied uniformly to every Tier A query at manifest-freeze time (no-op for
A1/A3 which are smegmatis-only by construction).

## A2 — 2026-08-19T14:2xZ — run metadata via esummary instead of efetch runinfo

`efetch.fcgi?db=sra&rettype=runinfo&retmode=csv` returns a **headerless** CSV
with no uid column, making exact uid-to-run attribution impossible in batches.
Switched the metadata fetch to `esummary.fcgi?db=sra&retmode=json` (uid,
`runs` blob with Run accs, expxml with ScientificName / LIBRARY_STRATEGY /
LIBRARY_SOURCE / Platform / Study). Provenance-affecting tooling change only;
no change to any query string, filter, threshold, or budget.

## A3 — 2026-08-19T14:3xZ — Tier A completion with stable induction terms

Rationale: the preregistered A1–A4 set under-captures the task's Tier A
definition ("mitomycin-C / prophage-induction studies") because (a) the
`prophage`/`bacteriophage`/`lysogeny` compound queries are currently broken
on EUtils (Amendment A1) and (b) A1–A4 missed classic induction vocabulary.
Counts observed at 14:3xZ (stable across repeats, unlike the broken terms):

| id | term | count |
|---|---|---|
| A5 | `"Mycobacterium smegmatis"[Organism] AND ciprofloxacin[All Fields]` | 16 |
| A6 | `"Mycobacterium smegmatis"[Organism] AND induction[All Fields]` | 25 |
| A7 | `"Mycobacterium smegmatis"[Organism] AND induced[All Fields]` | 48 |
| A8 | `"Mycobacterium smegmatis"[Organism] AND phage[All Fields]` | 66 |

Amendment executed **before any Tier A download**; result set is the union
A1–A8 (NTM-filtered per Amendment A1's A4 rule). No change to thresholds,
budgets, panel, or hit definition.
