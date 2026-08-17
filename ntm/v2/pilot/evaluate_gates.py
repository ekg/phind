#!/usr/bin/env python3
"""Evaluate the preregistered Stage-1 gates (PILOT_PLAN.md section 4).

Deterministic, offline: reads only the run directory (manifest.json,
ledger.jsonl, normalized/*.normalized.tsv) — no network, no inference beyond
the preregistered rules. Emits a per-gate report and a machine-readable
gates.json; exit codes: 0 = GO, 3 = NO-GO (a gate failed),
2 = INCOMPLETE (results not all fetched yet — cannot evaluate).

Gates (preregistered, fixed — no threshold fishing):
  G1 positive controls: D29 (AF022214.2) and L5 (NC_001335.1) slices must
     each hit their own reference accession with kmer_coverage >= 0.7.
  G2 negative controls: shuffled controls must produce 0 hits with
     kmer_coverage >= 0.9 in GenBank_RefSeq (submission threshold was 0.7;
     sub-0.9 hits are reported but are not a gate failure per the letter of
     the preregistration — they are flagged for interpretation).
  G3 fetch/parse: >= 90% of submitted sessions fetched and parsed without
     malformed-response failures.
  G4 reconciliation: every submitted bait has session id + status; complete
     rows have zip sha256 + normalized TSV; no duplicate or orphan sessions
     (ledger submissions not present in the manifest).
"""
import argparse
import csv
import glob
import json
import os
import re
import sys

POSITIVE_EXPECT = {
    "POSCON_D29_AF022214_2_MID1200": ("AF022214", "D29"),
    "POSCON_L5_NC_001335_1_MID1200": ("NC_001335", "L5"),
}
NEGATIVE_PREFIX = "NTMBAIT_SHUF_CONTROL"
THR_POS = 0.7
THR_NEG = 0.9
MIN_FETCH_RATE = 0.9


def base_bait_id(bait_id: str) -> str:
    """Strip the threshold suffix the prep step appends (e.g. `_t70`)."""
    return re.sub(r"_t\d+$", "", bait_id)


def role_of(bait_id: str) -> str:
    base = base_bait_id(bait_id)
    if base in POSITIVE_EXPECT:
        return "positive_control"
    if base.startswith(NEGATIVE_PREFIX):
        return "negative_control"
    return "ntm_bait"


def read_tsv_rows(path: str):
    with open(path, newline="") as fh:
        return [{k: v for k, v in r.items()} for r in csv.DictReader(fh, delimiter="\t")]


def cov(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def find_tsv(norm_dir: str, bait_id: str):
    cands = sorted(glob.glob(os.path.join(norm_dir, f"{bait_id}*.normalized.tsv")))
    return cands[0] if cands else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--out-json", default=None)
    args = ap.parse_args()

    mpath = os.path.join(args.run_dir, "manifest.json")
    manifest = json.load(open(mpath))
    rows = manifest["rows"]
    norm_dir = os.path.join(args.run_dir, "normalized")
    report = {"run_id": manifest.get("run_id"), "gates": {}, "rows": []}

    # ---- G4 reconciliation -------------------------------------------------
    g4 = {"name": "ledger/manifest reconcile", "pass": True, "problems": []}
    seen_sessions = {}
    sessions = set()
    for r in rows:
        if r["status"] in ("submitted", "complete"):
            if not r.get("session_id"):
                g4["pass"] = False
                g4["problems"].append(f"{r['bait_id']}: status {r['status']} but no session_id")
            else:
                sid = r["session_id"]
                if sid in seen_sessions:
                    g4["pass"] = False
                    g4["problems"].append(
                        f"duplicate session {sid} on {seen_sessions[sid]} and {r['bait_id']}")
                seen_sessions[sid] = r["bait_id"]
                sessions.add(sid)
        if r["status"] == "complete":
            if not r.get("zip_sha256"):
                g4["pass"] = False
                g4["problems"].append(f"{r['bait_id']}: complete but no zip_sha256")
            if not find_tsv(norm_dir, r["bait_id"]):
                g4["pass"] = False
                g4["problems"].append(f"{r['bait_id']}: complete but no normalized TSV")
    lpath = os.path.join(args.run_dir, "ledger.jsonl")
    ledger_sessions = set()
    if os.path.exists(lpath):
        for line in open(lpath):
            try:
                e = json.loads(line)
            except ValueError:
                continue
            if e.get("event") == "submission_recorded" and e.get("session"):
                ledger_sessions.add(e["session"])
    for s in sorted(ledger_sessions - sessions):
        g4["pass"] = False
        g4["problems"].append(f"orphan ledger session not in manifest: {s}")
    report["gates"]["G4"] = g4

    # ---- per-row roll-up ---------------------------------------------------
    n_sub = n_complete = n_malformed = 0
    g1 = {"name": "positive controls hit own reference >= 0.7", "pass": True, "detail": []}
    g2 = {"name": "negative controls: 0 hits at 0.9", "pass": True, "detail": []}
    for r in rows:
        role = role_of(r["bait_id"])
        entry = {"bait_id": r["bait_id"], "role": role, "status": r["status"],
                 "session": r.get("session_id"), "n_hits": r.get("n_hits")}
        if r["status"] == "submitted" and r.get("session_id"):
            n_sub += 1
        if r["status"] == "complete":
            n_sub += 1
            n_complete += 1
            tsv = find_tsv(norm_dir, r["bait_id"])
            hits = read_tsv_rows(tsv) if tsv else []
            entry["tsv_rows"] = len(hits)
            if role == "positive_control":
                prefix, label = POSITIVE_EXPECT[base_bait_id(r["bait_id"])]
                own = [h for h in hits if str(h.get("acc", "")).startswith(prefix)]
                best = max((cov(h.get("kmer_coverage")) or 0.0) for h in own) if own else None
                ok = best is not None and best >= THR_POS
                entry["own_reference_best_cov"] = best
                entry["gate1"] = "PASS" if ok else "FAIL"
                if not ok:
                    g1["pass"] = False
                g1["detail"].append(
                    f"{label} ({prefix}*): best own-reference kmer_coverage="
                    f"{best} -> {'PASS' if ok else 'FAIL'}")
            elif role == "negative_control":
                hi = [h for h in hits if (cov(h.get("kmer_coverage")) or 0.0) >= THR_NEG]
                best = max((cov(h.get("kmer_coverage")) or 0.0) for h in hits) if hits else 0.0
                ok = not hi
                entry["hits_ge_0.9"] = len(hi)
                entry["max_cov"] = best
                entry["gate2"] = "PASS" if ok else "FAIL"
                if not ok:
                    g2["pass"] = False
                g2["detail"].append(
                    f"{r['bait_id']}: total_hits={len(hits)} hits>=0.9={len(hi)} "
                    f"max_cov={best} -> {'PASS' if ok else 'FAIL'}")
        if r.get("error") and "malformed" in str(r["error"]).lower():
            n_malformed += 1
            entry["malformed"] = True
        report["rows"].append(entry)

    report["gates"]["G1"] = g1
    report["gates"]["G2"] = g2

    # control gates are only evaluable once every control row is complete
    for key, gate, roles in (("G1", g1, ("positive_control",)),
                             ("G2", g2, ("negative_control",))):
        ctrl = [r for r in report["rows"] if r["role"] in roles]
        if ctrl and not all(r["status"] == "complete" for r in ctrl):
            gate["pass"] = None
            gate["evaluable"] = False
            gate["detail"].append("not evaluable yet: control rows not all complete")
        else:
            gate["evaluable"] = True

    # ---- G3 fetch/parse ----------------------------------------------------
    rate = (n_complete / n_sub) if n_sub else 0.0
    g3 = {"name": ">= 90% sessions fetched+parsed, no malformed",
          "submitted": n_sub, "complete": n_complete, "malformed": n_malformed,
          "rate": round(rate, 4)}
    g3["evaluable"] = n_sub > 0
    g3["pass"] = (rate >= MIN_FETCH_RATE and n_malformed == 0) if n_sub > 0 else None
    report["gates"]["G3"] = g3

    pending = [r for r in rows if r["status"] in ("pending", "submitted")]
    all_gates = all(g["pass"] for g in report["gates"].values())
    if pending:
        verdict = "INCOMPLETE"
        rc = 2
    elif all_gates:
        verdict = "GO"
        rc = 0
    else:
        verdict = "NO-GO"
        rc = 3
    report["verdict"] = verdict

    out = args.out_json or os.path.join(args.run_dir, "gates.json")
    with open(out, "w") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)

    print(f"run: {manifest.get('run_id')}  rows={len(rows)}  verdict={verdict}")
    for k in ("G1", "G2", "G3", "G4"):
        g = report["gates"][k]
        if g.get("pass") is None:
            state = "NOT EVALUABLE"
        else:
            state = "PASS" if g["pass"] else "FAIL"
        print(f"  {k} {state} — {g['name']}")
        for d in g.get("detail", []) or g.get("problems", []):
            print(f"      {d}")
    print(f"gates json -> {out}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
