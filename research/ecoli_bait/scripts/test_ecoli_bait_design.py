#!/usr/bin/env python3
"""Tests for research/ecoli_bait/scripts/design_ecoli_baits.py.

Covers the E. coli-specific preregistered rules (DESIGN.md):
  - functional-QC tier derivation (A/B/ineligible)
  - traversal budget formula
  - purely-reconstructed junction selection (observed_adjacency_count == 0)
    with co-occurrence ranking
  - traversal verification gate (mismatch -> zero junction baits)
  - public-reference positive controls (exact mid-genome slices)
  - Markov negative controls (count, composition match, shared-kmer gate)
  - deterministic host-blind selection over a synthetic triangle
  - self-audit failure modes for junction/positive-control rows

Run:  python3 -m pytest research/ecoli_bait/scripts/test_ecoli_bait_design.py
      (or: python3 research/ecoli_bait/scripts/test_ecoli_bait_design.py)
"""

from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "..", "..", "..", "ntm", "v2",
                                "bait", "scripts"))

import baitlib as bl  # noqa: E402
import design_ecoli_baits as de  # noqa: E402


# ---------------------------------------------------------------------------
# Tier derivation
# ---------------------------------------------------------------------------

def test_tier_A_is_flag_ok():
    assert de.derive_tier("ok") == "A"
    assert de.derive_tier("") == "A"


def test_tier_B_allows_only_host_genes_and_low_completeness():
    assert de.derive_tier("host_genes") == "B"
    assert de.derive_tier("low_completeness") == "B"
    assert de.derive_tier("low_completeness;host_genes") == "B"
    assert de.derive_tier("host_genes;low_completeness") == "B"


def test_tier_ineligible_flags():
    assert de.derive_tier("missing_core_structural") is None
    assert de.derive_tier("missing_core_structural;low_completeness") is None
    assert de.derive_tier("many_truncated") is None
    assert de.derive_tier("contamination;host_genes") is None
    assert de.derive_tier("low_completeness;many_truncated") is None


# ---------------------------------------------------------------------------
# Traversal budget formula
# ---------------------------------------------------------------------------

def test_traversal_budget_formula():
    class FakeInp:
        pass

    import json, tempfile
    with tempfile.TemporaryDirectory() as td:
        d = os.path.join(td, "c")
        os.makedirs(d)
        with open(os.path.join(d, "manifest.json"), "w") as fh:
            json.dump({"sequence_lengths": {"max": 66686}}, fh)
        # 1.6 * 66686 = 106697.6 -> 106697 > 30k floor
        assert de.traversal_budget(d) == 106697
        with open(os.path.join(d, "manifest.json"), "w") as fh:
            json.dump({"sequence_lengths": {"max": 2625}}, fh)
        # floor applies
        assert de.traversal_budget(d) == 30000
        with open(os.path.join(d, "manifest.json"), "w") as fh:
            json.dump({"sequence_lengths": {"max": 18750}}, fh)
        # exactly at the floor boundary: int(1.6*18750)=30000
        assert de.traversal_budget(d) == 30000


# ---------------------------------------------------------------------------
# Junction candidates: purely-reconstructed joins only + ranking
# ---------------------------------------------------------------------------

def _write_trav(td, pids, occs, lens, has_seq=True):
    import json
    parts = {str(p): {"occurrence": occs[p], "representative_len": lens[p],
                      "has_seq": has_seq} for p in pids}
    obj = {"partitions": parts, "genome": {"pids": list(pids),
                                           "length_bp": sum(lens.values())}}
    pj = os.path.join(td, "ml.traversal.json")
    with open(pj, "w") as fh:
        json.dump(obj, fh)
    return pj


def _write_bed(td, rows):
    bed = os.path.join(td, "partitions.bed")
    with open(bed, "w") as fh:
        for member, s, e, pid in rows:
            fh.write(f"{member}\t{s}\t{e}\t{pid}\n")
    return bed


def test_junction_candidates_purely_reconstructed_only_and_ranking():
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        # genome pid order 1..5; joins: (1,2) observed-contiguous in M1,
        # (3,4) carried by M1+M2 but never contiguous (purely reconstructed,
        # co-occurrence 2), (4,5) contiguous in M2.
        pj = _write_trav(td, [1, 2, 3, 4, 5],
                         {1: 9, 2: 8, 3: 7, 4: 6, 5: 5},
                         {1: 500, 2: 500, 3: 500, 4: 500, 5: 500})
        bed = _write_bed(td, [
            ("M1", 0, 500, 1), ("M1", 500, 1000, 2),
            ("M1", 5000, 5500, 3), ("M1", 9000, 9500, 4),
            ("M2", 100, 600, 3), ("M2", 4000, 4500, 4),
            ("M2", 4500, 5000, 5),
            ("M3", 0, 500, 5),
        ])
        got = de.junction_candidates(bed, pj)
        # purely-reconstructed joins: (3,4) co-occ 2 ranked first, then
        # (2,3) co-occ 1; observed joins (1,2) and (4,5) are excluded.
        assert [(c["partition_a"], c["partition_b"]) for c in got] == \
            [(3, 4), (2, 3)], got
        jc = got[0]
        assert jc["observed_adjacency_count"] == 0
        assert jc["co_occurrence_count"] == 2
        assert jc["breakpoint"] == 1501  # cumulative len of pids 1,2 + 1
        assert got[1]["breakpoint"] == 1001


def test_junction_candidates_ranking_prefers_higher_co_occurrence():
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        # two purely-reconstructed joins: (1,2) co-occ 3, (3,4) co-occ 1
        pj = _write_trav(td, [1, 2, 3, 4],
                         {1: 9, 2: 9, 3: 9, 4: 9},
                         {1: 500, 2: 500, 3: 500, 4: 500})
        rows = []
        for m in ("A", "B", "C"):
            rows.append((m, 1000, 1500, 1))
            rows.append((m, 9000, 9500, 2))   # co-occ 3, never contiguous
        rows += [("D", 2000, 2500, 3), ("D", 8000, 8500, 4)]  # co-occ 1
        bed = _write_bed(td, rows)
        got = de.junction_candidates(bed, pj)
        # (1,2) co-occ 3; (3,4) co-occ 1; (2,3) co-occ 0 (still purely
        # reconstructed; ranked last)
        assert [(c["partition_a"], c["partition_b"]) for c in got] == \
            [(1, 2), (3, 4), (2, 3)], got
        assert got[0]["co_occurrence_count"] == 3
        assert got[1]["co_occurrence_count"] == 1
        assert got[2]["co_occurrence_count"] == 0
        assert all(c["observed_adjacency_count"] == 0 for c in got)


# ---------------------------------------------------------------------------
# Traversal verification gate
# ---------------------------------------------------------------------------

def test_verify_traversal_match_and_mismatch():
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "x.ml.fa")
        with open(p, "w") as fh:
            fh.write(">c\nACGTACGTAC\n")
        ok, rlen, elen = de.verify_traversal(p, "ACGTACGTAC")
        assert ok and rlen == 10 and elen == 10
        ok, rlen, elen = de.verify_traversal(p, "ACGTACGTAA")
        assert not ok


# ---------------------------------------------------------------------------
# Positive controls: exact mid-genome slices
# ---------------------------------------------------------------------------

def test_public_reference_baits_midgenome_slices():
    import tempfile

    class FakeInp:
        public_refs = {}

    with tempfile.TemporaryDirectory() as td:
        seq = ("ACGT" * 1000)  # 4000 bp synthetic reference
        for acc, (label, exp) in {"NC_000999.9": ("X", 4000)}.items():
            p = os.path.join(td, f"{acc}.fa")
            with open(p, "w") as fh:
                fh.write(f">{acc} synthetic\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")
            FakeInp.public_refs[acc] = p
        saved = dict(de.PUBLIC_REFS)
        try:
            de.PUBLIC_REFS.clear()
            de.PUBLIC_REFS["NC_000999.9"] = ("X", 4000)
            baits = de.make_public_reference_baits(FakeInp())
            assert len(baits) == 1
            b = baits[0]
            assert b.bait_class == "control_public_reference"
            assert b.expected_behavior == "own_reference_hit"
            s, e = de.bl.centered_window(4000, 2001, de.POSCON_W)
            assert (b.start, b.end) == (s, e)
            assert b.length() == de.POSCON_W
            assert b.seq == seq[s - 1:e]
        finally:
            de.PUBLIC_REFS.clear()
            de.PUBLIC_REFS.update(saved)


# ---------------------------------------------------------------------------
# Markov negative controls: count + composition + shared-kmer gate
# ---------------------------------------------------------------------------

def test_shuffled_controls_count_and_composition():
    def mk(i):
        return de.new_bait(
            bait_id=f"ECBAIT_{i}_0000_INTERIOR_MODULE_01",
            bait_class="interior_module", clade_id=f"{i}_0000",
            genome_id=f"clade_{i}_0000_ML", contig=f"clade_{i}_0000_ML",
            start=1, end=1200, seq=bl.revcomp(
                "ACGTTTAGGCACGATCGATCGGATTACAGCATTGCA" * 37)[:1200],
            module="terminase", dominated=False, priority=(0, 1, 0))

    baits = [mk(0), mk(1)]  # two accepted module baits, distinct genomes
    for b in baits:
        de.compute_masks(b)
    host_union = np.empty(0, dtype=np.uint64)
    out = de.make_shuffled_controls(baits, host_union)
    assert len(out) == de.N_SHUFFLED_CONTROLS
    for nb in out:
        assert nb.bait_class == "control_shuffled_negative"
        assert nb.expected_behavior == "no_hits"
        assert len(nb.seq) == 1200
        dev = bl.composition_deviation(nb.seq, baits[0].seq)
        assert dev["mono_l1"] < 0.10   # composition-matched in expectation
        assert dev["dinuc_l1"] < 0.35
        # gated: negligible shared canonical 31-mers with the panel
        panel = np.unique(np.concatenate(
            [bl.canonical_codes(b.seq) for b in baits]))
        shared = np.intersect1d(np.unique(bl.canonical_codes(nb.seq)),
                                panel).size
        assert shared <= de.SHARED_KMER_CONTROL_MAX
        assert f"mono_l1={dev['mono_l1']:.4f}" in nb.notes


# ---------------------------------------------------------------------------
# Deterministic host-blind selection over a synthetic triangle
# ---------------------------------------------------------------------------

class _SynthTriangle:
    """Tiny in-memory distance oracle for selection tests."""

    def __init__(self, d):
        self.d = d
        self.index = {k: i for i, k in enumerate(sorted(
            {a for a, _ in d} | {b for _, b in d}))}

    def dist(self, a, b):
        if a == b:
            return 0.0
        return self.d.get((a, b), self.d.get((b, a), 1.0))


def test_selection_deterministic_and_diverse():
    import types
    # 6 clades: 0..2 nearly identical (medoids within 0.01), 3..5 far apart
    members = {f"c{i}": [f"p{i}a", f"p{i}b"] for i in range(6)}
    d = {("p0a", "p1a"): 0.01, ("p1a", "p2a"): 0.01, ("p0a", "p2a"): 0.01}
    for i in range(6):
        for j in range(6):
            if abs(i - j) >= 2:
                a, b = f"p{i}a", f"p{j}a"
                d[(a, b)] = 0.4 + 0.1 * abs(i - j)
    tri = _SynthTriangle(d)

    saved = de.SELECT_N_A, de.SELECT_N_B
    de.SELECT_N_A, de.SELECT_N_B = 3, 1
    try:
        pools = {"A": ["c0", "c1", "c3", "c4"], "B": ["c2", "c5"]}
        medoids = {c: members[c][0] for c in
                   set(pools["A"]) | set(pools["B"])}
        cand = sorted(set(pools["A"]) | set(pools["B"]))
        dist = {(a, b): tri.dist(medoids[a], medoids[b])
                for a in cand for b in cand if a != b}
        sel = bl.greedy_farthest_point(sorted(pools["A"]), dist,
                                       sorted(pools["A"])[0], 3)
        sel += bl.greedy_farthest_point(sorted(pools["B"]), dist,
                                        sorted(pools["B"])[0],
                                        4 - len(sel))
        assert sel == ["c0", "c3", "c4", "c5"] or sel[0] == "c0"
        # near-identical clades c0/c1 never both selected
        assert not ("c0" in sel and "c1" in sel)
        # deterministic: same input -> same output
        sel2 = bl.greedy_farthest_point(sorted(pools["A"]), dist,
                                        sorted(pools["A"])[0], 3)
        assert sel2 == sel[:3]
    finally:
        de.SELECT_N_A, de.SELECT_N_B = saved


# ---------------------------------------------------------------------------
# Self-audit failure modes (E. coli rules)
# ---------------------------------------------------------------------------

def _row(bait_id, cls, contig, start, end, sha, obs="", adj="", pa="", pb="",
         bp=""):
    return {"bait_id": bait_id, "bait_class": cls, "contig": contig,
            "start": start, "end": end, "length_bp": end - start + 1,
            "seq_sha256": sha, "observed_adjacency_count": obs,
            "adjacency_observed": adj, "partition_a": pa, "partition_b": pb,
            "breakpoint": bp, "n_kmers_total": 1000, "n_usable_kmers": 900,
            "gc_frac": 0.5}


def test_self_audit_rejects_observed_adjacency_junction():
    seq = "ACGT" * 300
    row = _row("ECBAIT_0_0000_JUNCTION_01", "junction", "clade_0_0000_ML",
               1, 1200, bl.sha256_seq(seq), obs="3", adj="yes",
               pa="10", pb="20", bp="600")
    try:
        de.self_audit([("ECBAIT_0_0000_JUNCTION_01", seq)], [row], {}, {}, {})
        assert False, "should have raised"
    except AssertionError as exc:
        assert "purely reconstructed" in str(exc)


def test_self_audit_accepts_purely_reconstructed_junction():
    seq = "ACGT" * 300
    row = _row("ECBAIT_0_0000_JUNCTION_01", "junction", "clade_0_0000_ML",
               1, 1200, bl.sha256_seq(seq[0:1200]), obs="0", adj="no",
               pa="10", pb="20", bp="600")
    ml = {"clade_0_0000_ML": seq}
    de.self_audit([("ECBAIT_0_0000_JUNCTION_01", seq[0:1200])], [row],
                  ml, {}, {})


def test_self_audit_rejects_poscon_slice_mismatch():
    seq = "ACGT" * 300
    row = _row("ECBAIT_POSCON_LAMBDA_01", "control_public_reference",
               "NC_001416.1", 1, 1200, bl.sha256_seq(seq[0:1200]))
    ref = {"NC_001416.1": "TTTT" * 300}
    try:
        de.self_audit([("ECBAIT_POSCON_LAMBDA_01", seq[0:1200])], [row],
                      {}, {}, ref)
        assert False, "should have raised"
    except AssertionError as exc:
        assert "reference coordinates" in str(exc)


def test_self_audit_rejects_length_out_of_bounds():
    seq = "ACGT" * 100  # 400 bp < 500 min
    row = _row("ECBAIT_0_0000_INTERIOR_MODULE_01", "interior_module",
               "clade_0_0000_ML", 1, 400, bl.sha256_seq(seq))
    try:
        de.self_audit([("ECBAIT_0_0000_INTERIOR_MODULE_01", seq)], [row],
                      {"clade_0_0000_ML": seq}, {}, {})
        assert False
    except AssertionError as exc:
        assert "length" in str(exc)


def test_gates_reject_high_host_like_fraction():
    b = de.new_bait(bait_id="x", bait_class="interior_module",
                    clade_id="0", genome_id="g", contig="g",
                    start=1, end=1200, seq="ACGT" * 300, module="terminase",
                    dominated=False, priority=(0, 0, 0))
    de.compute_masks(b)
    b.stats["n_host_like"] = int(b.stats["n_kmers_total"] * 0.5)
    acc, ledger = de.apply_gates([b])
    assert not acc
    assert "host_like_frac" in ledger[0]["reason"]


def run_all():
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except Exception as exc:  # noqa: BLE001
                fails += 1
                print(f"FAIL {name}: {exc!r}")
    return fails


if __name__ == "__main__":
    sys.exit(1 if run_all() else 0)
