#!/usr/bin/env python3
"""Offline tests for ntm/v2/pilot/evaluate_gates.py (preregistered Stage-1 gates).

Fabricates run directories (manifest + normalized TSVs + ledger) on disk and
asserts verdicts. No network access anywhere.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
LOGAN = HERE.parent
REPO = LOGAN.parent
EVAL = REPO / "ntm" / "v2" / "pilot" / "evaluate_gates.py"

D29 = "POSCON_D29_AF022214_2_MID1200"
L5 = "POSCON_L5_NC_001335_1_MID1200"
NEG1 = "NTMBAIT_SHUF_CONTROL_01_t70"
NEG2 = "NTMBAIT_SHUF_CONTROL_02_t70"
BAIT1 = "NTMBAIT_0_0000_INTERIOR_MODULE_02_t70"


def write_tsv(run: Path, bait: str, rows) -> None:
    (run / "normalized").mkdir(parents=True, exist_ok=True)
    cols = ["acc", "kmer_coverage"]
    lines = ["\t".join(cols)]
    for acc, c in rows:
        lines.append(f"{acc}\t{c}")
    (run / "normalized" / f"{bait}.normalized.tsv").write_text("\n".join(lines) + "\n")


def make_run(tmp_path: Path, *, negatives="empty", d29="own", l5="own",
             orphan_session=False, statuses=None) -> Path:
    run = tmp_path / "run"
    run.mkdir()
    rows = []
    spec = [
        (BAIT1, [("SRR123", 0.81)]),
        (D29, [("AF022214.2", 0.93)] if d29 == "own" else [("SRR999", 0.42)]),
        (L5, [("NC_001335.1", 0.88)] if l5 == "own" else []),
        (NEG1, [] if negatives == "empty" else [("SRR555", 0.95), ("SRR556", 0.72)]),
        (NEG2, [] if negatives == "empty" else [("SRR557", 0.71)]),
    ]
    for bait, hits in spec:
        row = {"bait_id": bait, "status": "complete", "error": None,
               "session_id": f"kmviz-{abs(hash(bait)) % 10**12:012d}-0000-4000-8000-000000000000",
               "zip_sha256": "a" * 64, "n_hits": len(hits)}
        if statuses and bait in statuses:
            row["status"] = statuses[bait]
        rows.append(row)
        if row["status"] == "complete":
            write_tsv(run, bait, hits)
    (run / "manifest.json").write_text(json.dumps(
        {"run_id": "test", "rows": rows}, indent=2))
    ledger = [{"event": "submission_recorded", "session": r["session_id"],
               "bait_id": r["bait_id"]} for r in rows if r.get("session_id")]
    if orphan_session:
        ledger.append({"event": "submission_recorded",
                       "session": "kmviz-orphan-0000", "bait_id": "ghost"})
    (run / "ledger.jsonl").write_text(
        "\n".join(json.dumps(e) for e in ledger) + "\n")
    return run


def run_eval(run: Path):
    return subprocess.run([sys.executable, str(EVAL), "--run-dir", str(run)],
                          capture_output=True, text=True)


def test_all_gates_pass_go(tmp_path):
    run = make_run(tmp_path)
    p = run_eval(run)
    assert p.returncode == 0, p.stdout + p.stderr
    assert "verdict=GO" in p.stdout
    rep = json.loads((run / "gates.json").read_text())
    assert all(g["pass"] for g in rep["gates"].values())


def test_negative_control_hit_at_09_is_nogo(tmp_path):
    run = make_run(tmp_path, negatives="hit")
    p = run_eval(run)
    assert p.returncode == 3
    rep = json.loads((run / "gates.json").read_text())
    assert rep["verdict"] == "NO-GO"
    assert not rep["gates"]["G2"]["pass"]
    assert rep["gates"]["G2"]["detail"][0].startswith(NEG1)
    assert rep["gates"]["G2"]["detail"][1].startswith(NEG2)  # 0.71 max -> PASS


def test_positive_control_missing_own_reference_fails_g1(tmp_path):
    run = make_run(tmp_path, d29="wrong")
    p = run_eval(run)
    assert p.returncode == 3
    rep = json.loads((run / "gates.json").read_text())
    assert not rep["gates"]["G1"]["pass"]
    assert "D29" in rep["gates"]["G1"]["detail"][0]


def test_orphan_ledger_session_fails_reconcile(tmp_path):
    run = make_run(tmp_path, orphan_session=True)
    p = run_eval(run)
    assert p.returncode == 3
    rep = json.loads((run / "gates.json").read_text())
    assert not rep["gates"]["G4"]["pass"]
    assert any("orphan" in s for s in rep["gates"]["G4"]["problems"])


def test_threshold_suffix_ids_recognized(tmp_path):
    """Real manifest ids carry a _t70 suffix; roles must still resolve."""
    run = make_run(tmp_path, d29="wrong", negatives="hit")
    m = json.loads((run / "manifest.json").read_text())
    for r in m["rows"]:
        r["bait_id"] = r["bait_id"] + "_t70" \
            if not r["bait_id"].endswith("_t70") else r["bait_id"]
    (run / "manifest.json").write_text(json.dumps(m, indent=2))
    for r in m["rows"]:
        src = run / "normalized" / (r["bait_id"].removesuffix("_t70") + ".normalized.tsv")
        if src.exists():
            src.rename(run / "normalized" / (r["bait_id"] + ".normalized.tsv"))
    p = run_eval(run)
    assert p.returncode == 3
    rep = json.loads((run / "gates.json").read_text())
    roles = {r["bait_id"]: r["role"] for r in rep["rows"]}
    assert roles["POSCON_D29_AF022214_2_MID1200_t70"] == "positive_control"
    assert roles["NTMBAIT_SHUF_CONTROL_01_t70"] == "negative_control"
    assert not rep["gates"]["G1"]["pass"] and not rep["gates"]["G2"]["pass"]


def test_pending_rows_report_incomplete(tmp_path):
    run = make_run(tmp_path, statuses={BAIT1: "pending", NEG2: "submitted"})
    p = run_eval(run)
    assert p.returncode == 2
    rep = json.loads((run / "gates.json").read_text())
    assert rep["verdict"] == "INCOMPLETE"
