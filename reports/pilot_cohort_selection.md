# Frozen nested 1,000-assembly pilot cohort

**Verdict: PASS (metadata-only selection).** Immutable release `pilot-cohorts-v1-8afc0ea03d9e50dc` freezes exactly nested rungs of 10, 100, 250, 500, and 1,000 distinct exact *E. coli* assembly revisions. It reused the validated Stage-B ten without changing bytes, identity, or order. This task made zero network requests, downloaded zero sequence payloads, and performed zero biological analyses.

External release:

```text
/home/erikg/phind-data/ecoli26k/v1/releases/select-freeze-1k/pilot-cohorts-v1-8afc0ea03d9e50dc/
```

Its `SHA256SUMS` digest is `0f34ba98b248cac9b8f9c954aeb8a5971ed6b3340b708d60267937de2bd58682`; that digest is the content of the last-written, fsynced `COMPLETE` marker. The release was promoted by same-filesystem atomic directory rename. Consumers must reject an absent `COMPLETE`, changed inventory, non-PASS applicable gate, or changed release/manifest byte, row, digest, or order contract.

## Immutable inputs and automatic gates

The workflow automatically verified, without a human wait or override:

- collection release `collection-v1-f7494b4b89d1382b`, tracked `release.json` SHA-256 `59c6907e2c053e9d8ac3df8d5eb820bab0097030a9259ca2c9354c47cb6642bf`, external inventory SHA-256 `bd27dcca6c1e33eda8ccaafd30224aae5a2b5770ee6a87246474d6d55844f0e9`, and all 26,077 unique exact-version rows;
- the 2,246-byte, 10-row Stage-B manifest SHA-256 `0d179cbafce2ba1fa14d1929a4acd6621810a335f25bcd7ec67dd2083eb101f6` in both collection and acquisition handoffs;
- validated acquisition release `canonical-cohort-010-v1-e71484de9994fc28`, tracked release SHA-256 `4cf1e5f7abb11d13dbae886543a343b0a57a389b46aa3df4ebc4fb14d280ff23`, and external inventory SHA-256 `96a40035c15684d4c3c12c88f8134c32c4df421eb9d138119581ab7473badc44`;
- engineering-control source reference `prophage-semantics-v1-f5619e221ff272ae`, release SHA-256 `6a8de2063e4e12c0f0f363ebe41aba03a2f463e8a45711ccbe7ebcdae581b728`, and inventory SHA-256 `ab86ba983f4e5e44823e321dbafb84b19b1c31667abd3fc6cc0102ab1bdcef31`.

The source-semantics engineering gates are PASS, but `extraction_eligibility=EXTRACTION_BLOCKED` and `consumer_action=REJECT`. That release was used only to label post-selection controls. No coordinate extraction was started or authorized.

Both immutable root files matched before and after publication and were never edited:

| Input | SHA-256 |
|---|---|
| `26k_ecoli_accession.txt` | `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5` |
| `26k_prophage1.csv` | `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996` |

The external release contains exact copies of the compressed 26,077-row collection assembly input, Stage-B manifest, and three pinned predecessor release JSON files.

## Selection design frozen before downstream results

The predecessor documents 25,291 current rows, 714 authoritative terminal suppressed rows, and 72 exact-check-valid rows with detailed metadata unavailable. Suppressed rows are retained losslessly in the full mapping but are ineligible, leaving 25,363 eligible exact revisions.

The frozen design has three strata:

| Stratum | Rows | Rule |
|---|---:|---|
| Stage-B certainty | 10 | Exact validated Stage-B bytes and order; inclusion probability 1 |
| Main phage-blind random | 25,353 | All other non-suppressed exact revisions |
| Terminal suppressed, ineligible | 714 | Retained in frame; inclusion probability 0 |

The main ordering is the lexical order of

```text
SHA256("pilot-cohorts-v1-main-srs-sha256-seed-20260724" NUL assembly_id NUL exact_accession.version)
```

with `assembly_id` as a collision tie-breaker. Seed SHA-256 is `0adacc0cabf18bc167495380b788e23c3bbcc538d5652702cc750771a03d9fe4`. The first ten positions are the exact Stage-B order; the next 990 are the first hash-ranked main units. Every rung is an exact cohort-order prefix.

Only frozen accession/assembly identity, terminal assembly status, input mapping, and Stage-B membership enter selection. Prophage presence, counts, transposable/taxonomy flags, clusters, and all other phage traits are explicitly forbidden by `selection_policy.json` and absent from selection code paths. A test mutates forbidden traits and reverses input order without changing the cohort.

Exact main-stratum inclusion probabilities are:

| Rung | Stage-B certainty | Main random | Suppressed |
|---:|---:|---:|---:|
| 10 | `1/1` | `0/25353` | `0/1` |
| 100 | `1/1` | `90/25353` | `0/1` |
| 250 | `1/1` | `240/25353` | `0/1` |
| 500 | `1/1` | `490/25353` | `0/1` |
| 1,000 | `1/1` | `990/25353` | `0/1` |

Each selected row records its exact inverse-probability design weight. Host-clade or prevalence inference must use one recorded design weight per assembly; synthetic controls are excluded and observed controls receive no extra multiplicity.

## Frozen manifests

| Manifest | Rows | Bytes (plain external) | SHA-256 (plain external) |
|---|---:|---:|---|
| `cohort-0010.tsv` | 10 | 3,983 | `83f8925eb261ea915192e3923fe8d858f09abfa6791a27dbda20530da44ec146` |
| `cohort-0100.tsv` | 100 | 36,573 | `13e203961a9fcec18a8a09e690582652d8085b2a386811e6c6a03184b9489182` |
| `cohort-0250.tsv` | 250 | 91,475 | `ba2cf2909ccf62a0c1944a76b522edc5600953511ec355479117b4a419acbc9f` |
| `cohort-0500.tsv` | 500 | 182,605 | `bb6497bff230ecb6987dc5cda865307524a7fd63653e12e0d54d0808afe15ecb` |
| `cohort-1000.tsv` | 1,000 | 365,970 | `265a1e7784a4d5db3ea3577892feba8173290518b6c621f7e5091dbad66bfe77` |
| `engineering-controls.tsv` | 26 | 10,126 | `b65494d6a9bb31e5112034a0a52b5f43fc7900cc1a785824008b678d5ee68eb1` |
| `frame.tsv` | 26,077 | 10,624,176 | `dcb64433e032fd85a07cddf21b8cc7b87bdea7495c57e75cf1f83cf489db7379` |

The ordered accession list digest is `d4b03eb45c00f0607bd56085f599b1320baec732890a733694f020b2280f7392`. Boundary identities at cohort orders 10/100/250/500/1,000 are respectively `GCF_000498835.2`, `GCF_003730815.1`, `GCF_020089695.1`, `GCF_006238575.1`, and `GCF_017570625.1`.

Tracked handoff `manifests/pilot-cohorts-v1/` contains byte-identical rung/control manifests, deterministic `frame.tsv.gz` SHA-256 `8f54da187b2cca1db7046a4eb7c6e42a5790d56c641ec1d6b5dc5db3fa5257c2`, release/policy/input/evidence JSON, the exact external inventory, and tracked `SHA256SUMS` SHA-256 `8b27bd12b54424950ddc8a3d7249f94ce021e8897665043ecf0b498754a05139`. It is 4.8 MiB total; no tracked file exceeds 10 MiB.

## Separate bounded engineering controls

Control labels were computed only after cohort order was frozen and have `selection_effect=NONE_POST_SELECTION_LABEL`. The 26 rows cover:

- observed exact revision ambiguity (version >1);
- observed Stage-B minimum/maximum assembly size and contig count;
- observed assembled-molecule and unplaced-scaffold source-contig roles;
- a bounded synthetic unsafe-ID fixture (`ctg#A`);
- all ten already-selected Stage-B assemblies as a bounded future interval-sentinel scope;
- synthetic start-edge, end-edge, short, long, and circular-wrap interval fixtures.

Unsafe-ID and interval fixtures are excluded from biological inference. Observed controls remain ordinary cohort members and may be used only once with their frozen design weight. Every interval control is `BLOCKED_EXTRACTION_SEMANTICS` until a later immutable release provides unqualified extraction GO; none may be treated as coordinate evidence now.

## Atomic restart and deterministic evidence

The first execution self-sent `SIGKILL` after writing/fsyncing staged `SHA256SUMS` and before `COMPLETE`. It exited 137; neither final release nor final `COMPLETE` existed. Restart under the same run ID validated and reused only the checksum-complete selection unit, independently validated then discarded the incomplete publication stage, rebuilt it, wrote `COMPLETE` last, and atomically renamed it. `restart_evidence.json` and the append-only state ledger record all four PASS conditions.

A separate pre-publication rerun regenerated every selection manifest and matched every byte digest. The independent validator then re-derived the selection, all probabilities, exact nesting, row hashes, inventories, resource records, source block, compact handoff, global cap, and root hashes: 17/17 check groups PASS.

## Resource and global-cap evidence

Four live preflight records (selection/publication across kill and restart) contain `findmnt`, ownership/mode, write probes, free bytes/inodes, explicit allocations, and unfinished-write reservations. Assigned RAM was 4,294,967,296 bytes; peak RSS was 146,800,640 bytes (3.42%), with zero swap growth. The pre-seal stage was 14,610,135 bytes / 16 files against a 100 MB upper-95 predicted peak, 1 GB durable allocation, 4 TB scratch allocation, and 100,000-inode allocation. Durable 2 TB/1M-inode and scratch 4 TB/5M-inode preflight floors, 2× unfinished-write reserves, 70% disk/RAM, and 50% inode bounds all passed.

The graph-wide sequence-bearing union at start and finish remains exactly the validated Stage-B ten, all contained in this frozen cohort. This task added zero sequence-bearing identities; the exact global union is 10 and the hard cap is 1,000. BGZF/index/name, coordinate extraction, and scale-slope gates are explicitly not applicable to metadata-only selection, not mislabeled PASS.

## Reproduction

```bash
python -m unittest -v workflow.pilot_selection.test_selection
python workflow/pilot_selection/validate_release.py --repo-root . \
  --external-release /home/erikg/phind-data/ecoli26k/v1/releases/select-freeze-1k/pilot-cohorts-v1-8afc0ea03d9e50dc
```

`workflow/pilot_selection/validation.json` is the compact independent PASS record. Any checksum, identity, row accounting, source-status, resource, deterministic, restart, nesting, probability, cap, or root-input failure is a hard `NO_GO`; there is no human waiver path.
