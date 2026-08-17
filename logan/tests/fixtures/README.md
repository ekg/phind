# Fixtures

- `smoke_ecoli_k12.zip` — RAW response bytes of the bounded live smoke query
  (E. coli K-12 MG1655 NC_000913.3:1000000-1000149, 150 nt; group
  GenBank_RefSeq; threshold 0.5; session recorded in
  ../../smoke/run-2026-08-17/manifest.json), fetched once via the documented
  GET /api/download/<session> endpoint on 2026-08-17 and cached verbatim.
  sha256 recorded in the run ledger and the .meta.json sidecar.
- All other test responses are synthesized in-test (see make_zip in
  tests/test_logan_client.py); no test ever contacts the network.
