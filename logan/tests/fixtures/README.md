# Fixtures

- `smoke_ecoli_k12.zip` — **RESERVED, currently absent.** Intended: raw ZIP
  bytes of the first successful live smoke retrieval via the documented
  `GET /api/download/<session>` endpoint. The 2026-08-17 smoke submission's
  session never became retrievable (14 bounded polls, all HTTP 400 with the
  kmviz unknown-session error signature; see
  `../../smoke/run-2026-08-17/OUTCOME.md`), so no live bytes exist and none
  were fabricated. `TestRealSmokeFixture` in `tests/test_logan_client.py`
  skips until this file appears (Stage 0 of `PILOT_PLAN.md`).
- All other test responses are synthesized in-test (see make_zip in
  `tests/test_logan_client.py`); no test ever contacts the network.
