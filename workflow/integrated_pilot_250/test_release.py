#!/usr/bin/env python3
"""Tests for the integrated 250-genome scale-bearing workflow.

Covers the mandated task-scoped failure modes:
  * manifest/checksum mismatch
  * resource refusal (blank/non-positive allocation, >70%/>50%)
  * interrupted promotion/resume (mixed-resume refusal)
plus scale-trend math and a real integration build that must reach PASS/GO_500.
"""

import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from workflow.integrated_pilot_250 import release as R  # noqa: E402

# The live preflight requires >=4 TB scratch (nvme3n1) and >=2 TB durable (root FS).
# Tests must place scratch on the real nvme3n1 volume; durable may stay on tmp_path.
_TEST_SCRATCH_ROOT = Path("/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/_wg_test_250")


@pytest.fixture()
def nvme_scratch(tmp_path):
    _TEST_SCRATCH_ROOT.mkdir(parents=True, exist_ok=True)
    d = _TEST_SCRATCH_ROOT / f"run-{os.getpid()}-{tmp_path.name}"
    d.mkdir(parents=True, exist_ok=True)
    yield d
    shutil.rmtree(d, ignore_errors=True)


# --------------------------------------------------------------------------- #
# Pure helper / failure-mode tests
# --------------------------------------------------------------------------- #
class TestRootInputImmutability:
    def test_root_input_constants_match_disk(self):
        assert R.sha_file(REPO / "26k_ecoli_accession.txt") == R.ACCESSIONS_SHA
        assert R.sha_file(REPO / "26k_prophage1.csv") == R.SOURCE_SHA

    def test_global_cap_is_1000(self):
        assert R.GLOBAL_CAP == 1000


class TestChecksumMismatch:
    def test_verify_exact_missing_raises(self, tmp_path):
        with pytest.raises(R.GateError, match="required input missing"):
            R.verify_exact(tmp_path / "nope.json", "0" * 64, "missing")

    def test_verify_exact_mismatch_raises(self, tmp_path):
        p = tmp_path / "f.json"
        p.write_text("not the bytes")
        with pytest.raises(R.GateError, match="checksum mismatch"):
            R.verify_exact(p, "0" * 64, "wrong")

    def test_inventory_mismatch_raises(self, tmp_path):
        root = tmp_path / "rel"
        (root).mkdir()
        (root / "a.txt").write_text("hello")
        sums = tmp_path / "SHA256SUMS"
        sums.write_text(f"{'0' * 64}  rel/a.txt\n")
        with pytest.raises(R.GateError, match="inventory mismatch"):
            R.verify_inventory(root.parent, sums)


class TestResourceRefusal:
    def test_blank_allocation_refused(self):
        a = R.Allocations(0, 1, 1, 1, 1, 1, 1, 1)
        with pytest.raises(R.GateError, match="blank or non-positive"):
            a.validate()

    def test_durable_peak_over_70pct_refused(self):
        a = R.Allocations(100, 100, 10_000, 1_000_000, 80, 1, 1, 1)  # 80 > 70% of 100
        with pytest.raises(R.GateError, match="durable upper-95% peak exceeds 70"):
            a.validate()

    def test_files_over_50pct_inodes_refused(self):
        a = R.Allocations(100, 10_000, 10_000, 100, 1, 1, 60, 1)  # 60 > 50% of 100
        with pytest.raises(R.GateError, match="projected files exceed 50% inode"):
            a.validate()


class TestInterruptedResume:
    def test_mixed_resume_refused_on_digest_mismatch(self, tmp_path):
        stage = tmp_path / "stage"
        stage.mkdir()
        (stage / "state.jsonl").write_text("")
        unit = stage / "unit.json"
        unit.write_text('{"verdict": "OLD"}\n')
        new_bytes = R.canonical_bytes({"verdict": "NEW"})
        with pytest.raises(R.GateError, match="refusing mixed resume"):
            R.write_static_unit(stage, "unit.json", new_bytes)

    def test_clean_resume_validates_identical_unit(self, tmp_path):
        stage = tmp_path / "stage"
        stage.mkdir()
        (stage / "state.jsonl").write_text("")
        payload = R.canonical_bytes({"verdict": "SAME"})
        R.write_static_unit(stage, "unit.json", payload)
        # second call with identical bytes must not raise and must log RESUME
        R.write_static_unit(stage, "unit.json", payload)
        events = [json.loads(l)["event"] for l in (stage / "state.jsonl").read_text().splitlines() if l.strip()]
        assert "UNIT_COMMITTED" in events and "RESUME_UNIT_VALIDATED" in events


# --------------------------------------------------------------------------- #
# Scale-trend math
# --------------------------------------------------------------------------- #
class TestScaleTrendMath:
    def _fake_allocations(self):
        return R.Allocations(
            assigned_ram_bytes=68_719_476_736,
            durable_allocation_bytes=214_748_364_800,
            scratch_allocation_bytes=2_199_023_255_552,
            inode_allocation=1_000_000,
            predicted_durable_peak_bytes=107_374_182_400,
            predicted_scratch_peak_bytes=1_099_511_627_776,
            predicted_files=10_000,
            unfinished_write_bytes=10_737_418_240,
        )

    def test_linear_scaling_passes(self, monkeypatch, tmp_path):
        # Deterministic per-assembly QC walls: prior=0.050s, current=0.125s (~2.5x objects, ~linear)
        monkeypatch.setattr(R, "_timed_qc_wall", lambda *a, **k: 0.050)
        scale, go = R.compute_scale_trend(
            REPO, 0.125, 115_343_360, 64_000, 28, self._fake_allocations(), 7200.0,
        )
        assert go == "GO_500"
        assert scale["verdict"] == "PASS"
        assert scale["time_exponent"]["current_n100_to_n250"] <= 1.3
        gating = scale["last_two_rung_per_base_slopes"]["gating_metric"]
        assert abs(scale["last_two_rung_per_base_slopes"]["relative_changes"][gating]) <= 0.25

    def test_superlinear_exponent_blocks(self, monkeypatch, tmp_path):
        # current wall 1.0s vs prior 0.01s over 2.5x objects -> exponent > 1.3 -> NO_GO
        monkeypatch.setattr(R, "_timed_qc_wall", lambda *a, **k: 0.01)
        scale, go = R.compute_scale_trend(
            REPO, 1.0, 115_343_360, 64_000, 28, self._fake_allocations(), 7200.0,
        )
        assert go == "NO_GO"
        assert scale["time_exponent"]["current_n100_to_n250"] > 1.3


# --------------------------------------------------------------------------- #
# Predecessor identity (frozen)
# --------------------------------------------------------------------------- #
class TestPredecessorIdentity:
    def test_canonical_250_release_json_sha(self):
        assert R.sha_file(REPO / "manifests/canonical-cohort-250-v1/release.json") == R.CANONICAL_COHORT_RELEASE_JSON_SHA

    def test_prior_integrated_release_json_sha(self):
        p = R.PRIOR_INTEGRATED_EXTERNAL_ROOT / R.PRIOR_INTEGRATED_RELEASE_ID / "release.json"
        assert R.sha_file(p) == R.PRIOR_INTEGRATED_RELEASE_JSON_SHA

    def test_prophage_semantics_release_json_sha(self):
        p = R.PROPHAGE_SEMANTICS_EXTERNAL_ROOT / R.PROPHAGE_SEMANTICS_RELEASE_ID / "release.json"
        assert R.sha_file(p) == R.PROPHAGE_SEMANTICS_RELEASE_JSON_SHA


# --------------------------------------------------------------------------- #
# Integration: real build -> PASS / GO_500
# --------------------------------------------------------------------------- #
class TestIntegrationBuild:
    def _common_args(self, durable, scratch, run_id, extra=None):
        common = [
            "run", "--repo", str(REPO), "--durable-root", str(durable),
            "--scratch-root", str(scratch), "--run-id", run_id,
            "--assigned-ram-bytes", "68719476736", "--durable-allocation-bytes", "214748364800",
            "--scratch-allocation-bytes", "2199023255552", "--inode-allocation", "1000000",
            "--predicted-durable-peak-bytes", "107374182400", "--predicted-scratch-peak-bytes", "1099511627776",
            "--predicted-files", "10000", "--unfinished-write-bytes", "10737418240",
        ]
        if extra:
            common += extra
        return R.parser().parse_args(common)

    def test_build_passes_and_authorizes_go_500(self, tmp_path, nvme_scratch):
        durable = tmp_path / "durable"
        durable.mkdir()
        result = R.build(self._common_args(durable, nvme_scratch, "pilot-250-001"))
        assert result["verdict"] == "PASS"
        assert result["go_500"] == "GO_500"
        release_dir = next(durable.glob("integrated-pilot-250-v1-*"))
        assert (release_dir / "COMPLETE").is_file()
        assert (release_dir / "scale_trend.json").is_file()

    def test_idempotent_rerun_validates_existing(self, tmp_path, nvme_scratch):
        durable = tmp_path / "durable"
        durable.mkdir()
        first = R.build(self._common_args(durable, nvme_scratch, "pilot-250-002"))
        rid = first["release_id"]
        second = R.build(self._common_args(durable, nvme_scratch, "pilot-250-002"))
        assert second["release_id"] == rid
        assert second.get("deterministic_rerun") == "EXISTING_IMMUTABLE_RELEASE_VALIDATED"

    def test_injected_build_query_stop_then_resume(self, tmp_path, nvme_scratch):
        durable = tmp_path / "durable"
        durable.mkdir()
        # First pass: inject interruption at the new build/query stage (after SYNG)
        with pytest.raises(R.InjectedInterruption):
            R.build(self._common_args(durable, nvme_scratch, "pilot-250-003", ["--inject-stop-at-build-query"]))
        # No COMPLETE / no final release published yet (partial never exposed)
        assert not any(durable.glob("integrated-pilot-250-v1-*"))
        assert list(durable.glob(".*.staging")), "staging dir must exist after injected stop"
        # Second pass: resume -> static units SHA-validated, no mixed publication, reaches COMPLETE
        result = R.build(self._common_args(durable, nvme_scratch, "pilot-250-003"))
        assert result["verdict"] == "PASS"
        assert result["go_500"] == "GO_500"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
