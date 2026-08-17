#!/usr/bin/env python3
"""Unit tests for the NTM bait panel design (task design-ntm-prophage).

Covers the preregistered validation contract: reverse complements, boundary
coordinates, duplicated 31-mers, masking rules, deterministic selection,
window clamping, module-gene targeting and Markov-control composition gates.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import baitlib as bl  # noqa: E402

INVALID = np.uint64(1) << np.uint64(62)


# ---------------------------------------------------------------------------
# Reverse complements / canonical codes
# ---------------------------------------------------------------------------

def test_revcomp_basic():
    assert bl.revcomp("AACCGGTT") == "AACCGGTT"
    assert bl.revcomp("AAAAT") == "ATTTT"
    assert bl.revcomp("ACGTN") == "NACGT"


def test_canonical_code_is_strand_symmetric():
    s = "ACGTACGTACGTACGTACGTACGTACGTACG"  # 31-mer
    rc = bl.revcomp(s)
    assert np.unique(bl.canonical_codes(s)).size == 1
    assert bl.canonical_codes(s)[0] == bl.canonical_codes(rc)[0]


def test_canonical_codes_exact_vs_bruteforce():
    """62-bit codes are bijective; verify against a direct base-4 encoding."""
    rng = np.random.default_rng(7)
    seq = "".join("ACGT"[i] for i in rng.integers(0, 4, 60))
    codes = bl.canonical_codes(seq)

    def brute(kmer):
        code = 0
        for c in kmer:
            code = code * 4 + "ACGT".index(c)
        rcc = 0
        for c in bl.revcomp(kmer):
            rcc = rcc * 4 + "ACGT".index(c)
        return min(code, rcc)

    for i in range(len(seq) - 30):
        assert codes[i] == brute(seq[i:i + 31])


def test_canonical_codes_ambiguous_marked():
    seq = "A" * 10 + "N" + "ACGT" * 5 + "A" * 11  # 42 bp, one N at index 10
    codes = bl.canonical_codes(seq)
    # windows [0..10] contain the N → 11 invalid positions, all others valid
    assert (codes == INVALID).sum() == 11
    assert (codes[:11] == INVALID).all()
    assert (codes[11:] != INVALID).all()


# ---------------------------------------------------------------------------
# Boundary coordinates / window clamping
# ---------------------------------------------------------------------------

def test_centered_window_interior():
    # start = center - width//2 (documented convention; 1-based inclusive)
    assert bl.centered_window(10000, 5000, 1000) == (4500, 5499)


def test_centered_window_clamps_at_start():
    assert bl.centered_window(10000, 5, 1000) == (1, 1000)


def test_centered_window_clamps_at_end_and_preserves_width():
    s, e = bl.centered_window(10000, 9999, 1000)
    assert (s, e) == (9001, 10000)
    assert e - s + 1 == 1000


def test_centered_window_short_source():
    assert bl.centered_window(300, 150, 1000) == (1, 300)


def test_junction_window_spans_breakpoint():
    glen, bp = 20000, 13742
    s, e = bl.centered_window(glen, bp, bl.K * 2)
    assert s <= bp <= e + 1
    assert e - s + 1 == 62


def test_junction_window_at_genome_edge_still_spans():
    glen, bp = 20000, 3
    s, e = bl.centered_window(glen, bp, 1000)
    assert s <= bp <= e + 1  # clamped inward, join still covered


# ---------------------------------------------------------------------------
# Duplicated 31-mers
# ---------------------------------------------------------------------------

def test_internal_dup_mask_keeps_first_occurrence():
    rng = np.random.default_rng(31)
    dup = "".join("ACGT"[i] for i in rng.integers(0, 4, 31))  # aperiodic
    seq = dup + "TTTTTTTTTTTTTTTTTTTTTTTTTTTT" + dup  # 90 bp
    canon = bl.canonical_codes(seq)
    assert canon.size == 60
    dupmask = bl.internal_dup_mask(canon[canon != INVALID])
    assert dupmask.sum() == 1          # exactly one masked (the 2nd copy)
    # the second copy starts at seq[59:90] → final window index is masked,
    # the first occurrence (window 0) is kept
    assert dupmask.argmax() == canon.size - 1
    assert not dupmask[0]


def test_internal_dup_mask_reverse_complement_repeat():
    """A k-mer and its reverse complement share a canonical code — the
    second occurrence must be masked (strand-symmetric uniqueness)."""
    rng = np.random.default_rng(32)
    a = "".join("ACGT"[i] for i in rng.integers(0, 4, 31))
    seq = a + "G" * 31 + bl.revcomp(a)
    canon = bl.canonical_codes(seq)
    dupmask = bl.internal_dup_mask(canon[canon != INVALID])
    assert dupmask.sum() >= 1
    assert not dupmask[0]   # first occurrence kept
    assert dupmask[-1]      # the rc copy is the duplicate occurrence


def test_internal_dup_mask_no_dups():
    rng = np.random.default_rng(3)
    seq = "".join("ACGT"[i] for i in rng.integers(0, 4, 200))
    canon = bl.canonical_codes(seq)
    assert bl.internal_dup_mask(canon).sum() == 0


def test_cross_bait_dup_detection_by_intersection():
    rng = np.random.default_rng(11)
    a = bl.canonical_codes("".join("ACGT"[i] for i in rng.integers(0, 4, 1200)))
    b = np.append(a[:600], bl.canonical_codes(
        "".join("ACGT"[i] for i in rng.integers(0, 4, 600))))
    frac = np.intersect1d(np.unique(a), np.unique(b)).size / np.unique(a).size
    assert frac > 0.4  # heavily shared — must trip the cross-bait gate


# ---------------------------------------------------------------------------
# Masking rules
# ---------------------------------------------------------------------------

def test_ambiguous_mask_counts_windows_with_n():
    seq = "A" * 40 + "N" * 5 + "A" * 40
    amb = bl.ambiguous_mask(seq)
    assert amb.sum() == 35  # windows [i, i+31) covering any of the 5 Ns


def test_low_complexity_homopolymer():
    seq = "ACGT" * 10 + "A" * 40 + "ACGT" * 10
    lcx = bl.low_complexity_mask(seq)
    assert lcx.sum() > 0
    # homopolymer run is [40, 81): the run merges with the leading A of the
    # trailing ACGT block. Every masked window must overlap that run.
    run_start, run_end = 40, 81
    for i in np.flatnonzero(lcx):
        assert i <= run_end - 1 and i + 30 >= run_start
    # and windows clear of the run are unmasked (entropy 2 bits, no runs)
    assert not lcx[:5].any()


def test_low_complexity_entropy():
    seq = "AC" * 300  # alternating dinucleotide: low entropy despite 50/50 GC
    lcx = bl.low_complexity_mask(seq)
    assert lcx.all()


def test_low_complexity_random_dna_unmasked():
    rng = np.random.default_rng(5)
    seq = "".join("ACGT"[i] for i in rng.integers(0, 4, 1000))
    assert bl.low_complexity_mask(seq).sum() == 0


def test_host_like_scan_finds_shared_kmers():
    rng = np.random.default_rng(9)
    bait = "".join("ACGT"[i] for i in rng.integers(0, 4, 1200))
    host = "".join("ACGT"[i] for i in rng.integers(0, 4, 3000)) \
        + bait[100:900]
    bait_canon = np.unique(bl.canonical_codes(bait))
    hit = np.isin(bl.canonical_codes(host), bait_canon)
    assert hit.sum() >= 700  # the embedded bait segment is detected
    # and the scan counter agrees
    assert bl.scan_kmers_against(host, np.sort(bait_canon)) >= 700


def test_host_scan_respects_strand():
    """A reverse-complemented bait segment must still be flagged host-like."""
    rng = np.random.default_rng(13)
    bait = "".join("ACGT"[i] for i in rng.integers(0, 4, 1200))
    host = bl.revcomp(bait[200:1000])
    bait_canon = np.sort(np.unique(bl.canonical_codes(bait)))
    canon_host = bl.canonical_codes(host)
    shared = np.intersect1d(canon_host, np.unique(bl.canonical_codes(bait)))
    assert shared.size == canon_host.size


# ---------------------------------------------------------------------------
# Region-level annotation gates
# ---------------------------------------------------------------------------

def _cds(start, stop, annot, category="tail"):
    return {"start": start, "stop": stop, "annot": annot,
            "category": category, "vfdb_hit": "None", "CARD_hit": "None"}


def test_domination_rejects_transposase_window():
    rows = [_cds(1, 600, "transposase"), _cds(601, 1200, "transposase OrfB")]
    dominated, mfrac, cfrac = bl.cds_composition_flags(rows, 1, 1200)
    assert dominated and mfrac == 1.0 and cfrac == 1.0


def test_domination_allows_module_genes():
    rows = [_cds(1, 900, "terminase large subunit", "head and packaging"),
            _cds(901, 1200, "hypothetical protein", "unknown function")]
    dominated, _, _ = bl.cds_composition_flags(rows, 1, 1200)
    assert not dominated


def test_domination_requires_cds_coverage():
    rows = [_cds(1, 300, "transposase")]
    dominated, _, cfrac = bl.cds_composition_flags(rows, 1, 1200)
    assert not dominated and cfrac < bl.CDS_COVER_FRAC_MIN


def test_housekeeping_category_counts_as_marker():
    rows = [_cds(1, 700, "DNA polymerase", "DNA, RNA and nucleotide metabolism"),
            _cds(701, 1200, "ribonucleotide reductase",
                 "DNA, RNA and nucleotide metabolism")]
    dominated, mfrac, _ = bl.cds_composition_flags(rows, 1, 1200)
    assert dominated and mfrac == 1.0


def test_pick_module_genes_priority():
    rows = [
        _cds(5000, 6200, "tail length tape measure protein", "tail"),
        _cds(100, 1300, "terminase large subunit", "head and packaging"),
        _cds(9000, 9800, "portal protein", "head and packaging"),
        _cds(2000, 2600, "major head protein", "head and packaging"),
    ]
    picked = bl.pick_module_genes(rows)
    assert [p["module"] for p in picked] == \
        ["terminase", "portal", "capsid", "tail_tape_measure"]


# ---------------------------------------------------------------------------
# Deterministic selection
# ---------------------------------------------------------------------------

def test_greedy_farthest_point_deterministic_and_diverse():
    cands = [f"c{i:02d}" for i in range(10)]
    # hand-built distances: a star with two far clusters
    dist = {}
    for i, a in enumerate(cands):
        for j, b in enumerate(cands):
            if a == b:
                continue
            d = 0.1 if (i < 5) == (j < 5) else 0.9
            dist[(a, b)] = d
    r1 = bl.greedy_farthest_point(cands, dist, "c00", 4)
    r2 = bl.greedy_farthest_point(cands, dist, "c00", 4)
    assert r1 == r2                      # deterministic
    assert r1[0] == "c00"
    assert any(c >= "c05" for c in r1)   # diversity-aware: crosses clusters


def test_greedy_farthest_point_tie_breaks_by_id():
    cands = ["b", "a", "c"]
    dist = {("a", "b"): 0.5, ("b", "a"): 0.5,
            ("a", "c"): 0.5, ("c", "a"): 0.5,
            ("b", "c"): 0.5, ("c", "b"): 0.5}
    assert bl.greedy_farthest_point(cands, dist, "a", 2) == ["a", "b"]


def test_clade_medoid_prefers_central_member():
    # three members: m2 is central (d(m1,m2)=0.1, d(m2,m3)=0.1, d(m1,m3)=0.9)
    members = ["m1", "m2", "m3"]
    sub_index = {m: i for i, m in enumerate(members)}
    D = np.full((3, 3), np.nan)
    D[0, 1] = D[1, 0] = 0.1
    D[1, 2] = D[2, 1] = 0.1
    D[0, 2] = D[2, 0] = 0.9
    assert bl.clade_medoid(members, sub_index, D) == "m2"


def test_clade_medoid_nan_distances_count_as_far():
    members = ["m1", "m2", "m3"]
    sub_index = {m: i for i, m in enumerate(members)}
    D = np.full((3, 3), np.nan)
    D[0, 1] = D[1, 0] = 0.2
    D[0, 2] = D[2, 0] = 0.2
    # m2/m3 NaN to each other → 1.0
    assert bl.clade_medoid(members, sub_index, D) == "m1"


# ---------------------------------------------------------------------------
# Markov negative controls
# ---------------------------------------------------------------------------

def test_markov_control_same_length_and_composition():
    rng = np.random.default_rng(21)
    src = "".join("ACGT"[i] for i in rng.integers(0, 4, 1200))
    a = bl.markov_control_sequence(src, seed=42)
    b = bl.markov_control_sequence(src, seed=42)
    assert a == b                                   # deterministic
    assert len(a) == len(src)
    dev = bl.composition_deviation(a, src)
    assert dev["mono_l1"] < 0.10                    # composition matched
    assert dev["dinuc_l1"] < 0.60
    shared = np.intersect1d(np.unique(bl.canonical_codes(a)),
                            np.unique(bl.canonical_codes(src))).size
    assert shared == 0                              # 31-mers destroyed


def test_markov_control_gc_tracks_source():
    src = "GC" * 600  # 100% GC source
    ctrl = bl.markov_control_sequence(src, seed=1)
    gc = ctrl.count("G") + ctrl.count("C")
    assert abs(gc - len(ctrl)) < len(ctrl) * 0.05  # stays ~100% GC


# ---------------------------------------------------------------------------
# FASTA round-trip and checksums
# ---------------------------------------------------------------------------

def test_fasta_roundtrip_and_sha(tmp_path):
    recs = [("bait1", "ACGT" * 25), ("bait2", "TTTT" * 25)]
    p = tmp_path / "t.fa"
    bl.write_fasta(recs, str(p), width=20)
    back = bl._parse_fasta(open(p))
    assert back == {k: v for k, v in recs}
    assert bl.sha256_seq("acgt" * 10) == bl.sha256_seq("ACGT" * 10)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
