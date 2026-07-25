import json
import unittest
from pathlib import Path

from workflow.prophage_semantics import pinned_caller_gate as gate
from workflow.prophage_semantics import release


ROOT = Path(__file__).resolve().parents[2]
PREDECESSOR_EXTERNAL_ROOT = gate.DEFAULT_EXTERNAL_ROOT


class PinnedCallerGateDecisionTests(unittest.TestCase):
    """Unit-test the pure decision rule in isolation (no live predecessor)."""

    def test_decisive_and_sound_authorizes_extraction(self):
        v, c = gate.decide_extraction("DECISIVE", True, "GO")
        self.assertEqual((v, c), ("EXTRACTION_GO", "ALLOW"))

    def test_decisive_and_sound_still_go_even_if_modern_pilot_is_no_go(self):
        # the modern v2.4 pilot is strictly separate and must not gate historical
        v, c = gate.decide_extraction("DECISIVE", True, "NO_GO")
        self.assertEqual((v, c), ("EXTRACTION_GO", "ALLOW"))

    def test_modern_go_alone_does_not_authorize_historical(self):
        v, c = gate.decide_extraction("NON_DECISIVE", True, "GO")
        self.assertEqual((v, c), ("EXTRACTION_BLOCKED", "REJECT"))

    def test_decisive_but_not_independently_sound_stays_blocked(self):
        v, c = gate.decide_extraction("DECISIVE", False, "GO")
        self.assertEqual((v, c), ("EXTRACTION_BLOCKED", "REJECT"))

    def test_failure_result_is_hard_no_go(self):
        result = gate.failure_result(
            ROOT, PREDECESSOR_EXTERNAL_ROOT,
            release.GateError("predecessor release is absent"),
        )
        self.assertEqual(result["verdict"], "NO_GO")
        self.assertTrue(result["hard_stop"])
        self.assertEqual(result["historical_csv_extraction"], "EXTRACTION_BLOCKED")
        self.assertEqual(result["consumer_action"], "REJECT")
        self.assertEqual(result["modern_v2_4_pilot"], "UNAVAILABLE")
        self.assertFalse(result["decisive_evidence_independently_sound"])


class PinnedCallerGateLiveTests(unittest.TestCase):
    """Integration test against the immutable predecessor external release.

    Bounded + read-only: reads only the predecessor's released native outputs
    and the immutable CSV.  No genome reads, no downloads.
    """

    def setUp(self):
        if not (PREDECESSOR_EXTERNAL_ROOT / gate.PREDECESSOR_RELEASE_ID / "COMPLETE").is_file():
            self.skipTest("predecessor pinned-caller release not present in this environment")

    def test_live_predecessor_passes_and_authorizes_extraction(self):
        result = gate.validate_pinned_caller_release(ROOT, PREDECESSOR_EXTERNAL_ROOT)
        self.assertEqual(result["verdict"], "PASS")
        self.assertEqual(result["release_id"], gate.PREDECESSOR_RELEASE_ID)
        self.assertEqual(result["historical_csv_attribution"], "DECISIVE")
        self.assertEqual(result["modern_v2_4_pilot"], "GO")
        self.assertTrue(result["modern_v2_4_pilot_separate"])
        self.assertTrue(result["decisive_evidence_independently_sound"])
        self.assertEqual(result["historical_csv_extraction"], "EXTRACTION_GO")
        self.assertEqual(result["consumer_action"], "ALLOW")
        # independent re-verification must show 56/56 exact and the +1 signature
        ind = result["independent_reverification"]
        self.assertEqual(ind["all_fields_exact_count"], 56)
        self.assertEqual(ind["csv_rows_for_cohort"], 56)
        self.assertEqual(ind["boundary_signature"]["begin_delta_v23_minus_v24_dist"], {1: 56})
        self.assertEqual(ind["boundary_signature"]["end_delta_v23_minus_v24_dist"], {1: 56})

    def test_live_predecessor_integrity_pinned(self):
        result = gate.validate_pinned_caller_release(ROOT, PREDECESSOR_EXTERNAL_ROOT)
        # SHA-256 inventory must round-trip; pin key predecessor digests
        self.assertEqual(result["inventory_rows"], 84)
        self.assertEqual(len(result["complete_sha256"]), 64)
        self.assertEqual(len(result["sha256sums_sha256"]), 64)
        self.assertEqual(len(result["release_json_sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
