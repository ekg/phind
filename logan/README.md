# Logan Search / IndexThePlanet — safe pilot client

**Task:** `implement-safe-logan` (WG). **Documentation retrieval date for every
live fact below: 2026-08-17 (UTC)** unless stated otherwise.

This directory contains a bounded, cached, rate-limited client for
[Logan Search](https://logan-search.org) and the retrieval/confirmation
workflow for sequence-bait discovery across SRA assemblies (NTM prophage
pilot). It implements the "safe client" contract: documented interfaces only,
hard bounds, full request/response provenance, resume, dry-run.

## 1. Official interfaces — what is (and is not) programmable

| Interface | Status (2026-08-17) | Evidence |
|---|---|---|
| Logan Search web dashboard `https://logan-search.org/dashboard` | **Only documented query-submission surface.** Human/browser-driven. Results retained **one month**; query/session IDs look like `kmviz-<uuid>` (e.g. `kmviz-08a28c6c-9691-4eca-a8aa-ef62f098f62c`). | docs.logan-search.org ("Submitting a query", "History") |
| kmviz REST **query** API (`POST /api/query`, `POST /api/query/<db>`) | **NOT registered on the public instance.** Live probes: `POST https://logan-search.org/api/query` → `405` (Dash GET-only catch-all); `GET /api` → Dash SPA HTML, not the kmviz infos JSON. The kmviz framework *can* host this API (`api.enabled`, `with_query`; kmviz `kmviz/api.py::_register`), but the public instance does not enable query routes. | kmviz source (github.com/tlemane/kmviz, `kmviz/api.py`, `kmviz/core/config.py` `capi`) + live probes 2026-08-17 |
| kmviz REST **download** API `GET /api/download/<session>` | **REGISTERED and working on the public instance** (i.e. `api.enabled=true, with_query=false`). Returns a ZIP with one TSV per query sequence + `session.json` + `README.md`. Unknown/not-ready session → HTTP 400, body `An error occured while processing your request` (kmviz generic API error string). | Live probe with bogus session → kmviz error signature; semantics from kmviz `api.py::_make_download_result` |
| `api.indextheplanet.com` | **Does not resolve** (NXDOMAIN). `indextheplanet.com`/`www.` CNAME to `logan-search.org` (A 20.62.90.26). The service has moved to `logan-search.org`. | DNS lookups 2026-08-17 |
| Logan data on AWS S3 (`logan-pub`) | Unrestricted, free; `https://s3.amazonaws.com/logan-pub/u/<acc>/<acc>.unitigs.fa.zst` and `.../c/<acc>/<acc>.contigs.fa.zst`. HEAD probes OK (e.g. `SRR17555654` unitigs = 2,290,735,510 bytes; `DRR000016` contigs = 189,571 bytes). | docs.logan-search.org (Introduction, download commands); IndexThePlanet/Logan `Unitigs.md`, `Contigs.md`; live HEAD 2026-08-17 |
| Local per-accession search tools | `back_to_sequences` (splits query into k-mers, reports matching unitig/contig sequences, `--output-mapping-positions` for coordinates; bioconda + b2s-doc.readthedocs.io), `rcgrep` (single k-mer, strand-aware), `LoganBlaster` (multi-sample alignment). | IndexThePlanet/Logan `Sequence_Search.md` |

### 1.1 Live-vs-docs discrepancies (recorded, not resolved by us)

* **Groups.** docs.logan-search.org lists 5 groups (`All`, `All_No_viral_human`,
  `Fast`, `Fast_No_human`, `GenBank_RefSeq`). The **live dashboard combobox
  (observed 2026-08-17)** lists **9**: the 5 documented plus
  `Fast_No_RefSeq`, `Transcriptomic`, `Metatranscriptomic`, `Metagenomic`.
  The client accepts all 9 and flags which are doc-documented.
* **Search-index snapshot vs dataset release.** docs.logan-search.org describes
  the search index as "all Logan unitigs, SRA up until 2023 - 23.4 million
  total samples" (+ GenBank/RefSeq references, ~45k). The GitHub repo announces
  dataset **v1.2** (SRA freeze **2025-12-31**; 38.1M accessions with unitigs,
  37.3M with contigs; 4.19 PB compressed unitigs / 623 TB contigs —
  `Stats-v1.2.md`, retrieved 2026-08-17). We could not verify from the docs
  whether the live search index covers v1.2 or still v1; treat the *index*
  snapshot as "SRA ≤ 2023 unless the service states otherwise" and record the
  group + date with every query.
* **Threshold default.** Live dashboard slider default = **0.5** (range 0.25–1.0),
  matching the documented bounds; docs do not state a default.

## 2. Documented constraints encoded in the client

From docs.logan-search.org ("Submitting a query"), retrieved 2026-08-17:

* **One FASTA sequence per submission**, maximum length **2.5 kb** (2500 nt).
  Minimum one k-mer ⇒ ≥ 31 nt (`k=31`).
* **k = 31**; the index is over **unitigs**; engine is **kmindex**
  (Bloom-filter k-mer search, built with kmtricks-logan).
* **Threshold** = minimum proportion of query k-mers shared with a sample,
  **0.25–1.0**; "no impact on the query time".
* **Result fields**: per-accession `kmer_coverage` (shared-k-mer ratio),
  `ANI_estimation` (Mash-screen estimator `(|Qk∩Sk|/|Qk|)^(1/k)`), Poisson
  `p-value`/`e-value`, SRA metadata (organism, bioproject, assay, …), and
  Logan assembly stats (contigs/unitigs N50, nbseq, maxlen, sumlen).
* **Sub-index ("group") choice**: `GenBank_RefSeq` (~45k samples) is the
  smallest/bounded group — used for smoke; the NTM pilot ladder uses
  `GenBank_RefSeq` first, then `Fast` (77.59% of samples, ~99.5% of non-viral)
  only if gates pass; `All` (2,869 sub-indexes, 23.4M samples) is NOT used in
  the bounded pilot.
* **k-mer semantics**: kmindex/kmtricks canonicalizes k-mers (strand-agnostic);
  unitigs contain every SRA k-mer with abundance ≥ 2 (docs Introduction);
  contigs only guarantee k-mers present in reads.

## 3. Service terms

* Logan sequences "are provided free of charge and are available for
  unrestricted download" (docs.logan-search.org, Introduction, retrieved
  2026-08-17).
* If using the Logan Search service, cite: Chikhi et al., *Logan:
  Planetary-Scale Genome Assembly Surveys Life's Diversity*, bioRxiv 2024,
  doi:10.1101/2024.07.30.605881 (stated on the submission page).
* S3 bucket `logan-pub` is on the AWS Registry of Open Data
  (`registry.opendata.aws/pasteur-logan/`): license follows NCBI/NIH data
  policies; requested citation: "Logan Unitigs and Contigs of the Sequence
  Read Archive (SRA) on AWS was accessed on DATE from
  https://registry.opendata.aws/pasteur-logan" (retrieved 2026-08-17).
* No published rate limits or quotas exist for Logan Search; we therefore
  self-impose: ≥ 30 s between any HTTP requests (60 s while polling), ≤ 40
  HTTP requests per run, one dashboard submission at a time, and no automated
  retry storms (bounded exponential backoff 5→10→20→40 s, cap 300 s).

## 4. Client design (`logan_search_client.py`, `run_pilot.py`)

* **Validation before any network**: single bare sequence, ACGTN only,
  31 ≤ len ≤ 2500, threshold ∈ [0.25, 1.0] rounded to 2 dp, group ∈ live set,
  unique (bait, group, threshold); plan size ≤ hard cap (default 40) —
  over-limit sequences/batches are **refused** with `LoganValidationError`.
* **Deterministic normalization**: uppercase, whitespace/`-`/`*` stripped,
  header rejection; result TSVs lower-cased columns, numeric coercion,
  sorted by `kmer_coverage` DESC then `acc` ASC; normalized TSV output has a
  stable column order (sorted) — byte-identical for identical inputs.
* **Cache**: raw ZIP per session id + sidecar meta (sha256, size, URL,
  params, fetch time); cache reads verify sha256 (tamper → `CacheIntegrityError`);
  cached sessions are never re-downloaded.
* **Rate limiting**: persistent `rate_limit_state.json`; minimum interval
  enforced across processes/restarts (protects resumed runs).
* **Ledger**: append-only `ledger.jsonl`; every request records ts, url,
  session, attempt, status, size, sha256, duration, rate-limit wait, params;
  failures and backoff sleeps are logged too.
* **Retries/backoff**: deterministic exponential backoff; HTTP 400 during a
  polling window (`--poll N`) is a soft "not ready", spaced ≥ min-interval.
* **Resume**: `manifest.json` rows move `pending → submitted → complete`
  (or `failed` with error); `fetch` touches the network only for `submitted`
  rows; completed rows are skipped; failed rows keep their error for re-run.
* **Dry-run**: `python3 run_pilot.py dryrun --fasta ... --group ...` prints
  the full plan (params per bait, index snapshot, bounds) with **zero**
  network side effects.
* **Submission path**: dashboard-assisted (human or browser driver) — the
  client writes `submissions/<bait>.fa` + `INSTRUCTIONS.md`; the returned
  `kmviz-<uuid>` is recorded offline via `record`; only then does `fetch`
  use the documented download endpoint.

## 5. Smoke query (executed 2026-08-17, single bounded submission)

* Query: 150 nt of *E. coli* K-12 MG1655 **NC_000913.3:1,000,000–1,000,149**
  (public BSL-1 reference; non-sensitive), retrieved from NCBI E-utilities
  2026-08-17; seq sha256 `ac58a6648da0ca4ac7f14728e50f3e110bc3e1c25cdc52b544650ec2b4c801ed`.
* Parameters: group `GenBank_RefSeq`, threshold 0.5 (dashboard default),
  no email. Submitted once via the dashboard (browser-driven); session
  `kmviz-9b877a85-62a8-46da-94d7-b056c3fda36d`.
* Results: fetched via `GET /api/download/<session>` (rate-limited), cached,
  parsed, normalized under `smoke/run-2026-08-17/`; raw ZIP + sha256 retained
  as test fixture `tests/fixtures/smoke_ecoli_k12.zip`.
* Expected sanity: hits in *E. coli* reference genomes at high k-mer coverage
  (positive control behavior). No biological claims beyond parsing mechanics.

## 6. Pilot execution plan (for `execute-bounded-ntm`)

See `PILOT_PLAN.md` for the full bounded plan: query-count and download-volume
estimates, threshold ladder (0.5/0.7/0.9 — all within documented bounds),
stop/go gates, S3 retrieval + `back_to_sequences`/minimap2 confirmation
workflow.

## 7. Files

| Path | Purpose |
|---|---|
| `logan_search_client.py` | library: validation, cache, rate limit, ledger, retries, parsing/normalization |
| `run_pilot.py` | CLI: `dryrun`, `prep`, `record`, `fetch`, `status` |
| `tests/test_logan_client.py` | unit tests (fixtures only; no network) |
| `tests/fixtures/` | cached real smoke response + synthetic fixtures |
| `smoke/run-2026-08-17/` | smoke run dir: manifest, ledger, cache, normalized results |
| `PILOT_PLAN.md` | bounded pilot execution plan with gates |

## 8. Reproduce

```bash
cd logan
python3 -m pytest tests/ -q                  # unit tests (offline)
python3 run_pilot.py dryrun --fasta baits.fa --group GenBank_RefSeq
python3 run_pilot.py prep    --run-dir runs/pilot1 --fasta baits.fa \
    --group GenBank_RefSeq --thresholds 0.5 0.7 0.9
# submit each runs/pilot1/submissions/*.fa via the dashboard (INSTRUCTIONS.md)
python3 run_pilot.py record  --run-dir runs/pilot1 --bait <id> --session kmviz-<uuid>
python3 run_pilot.py fetch   --run-dir runs/pilot1 --poll 6 --min-interval 60
python3 run_pilot.py status  --run-dir runs/pilot1
```
