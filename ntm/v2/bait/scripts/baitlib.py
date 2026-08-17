#!/usr/bin/env python3
"""baitlib — core library for the NTM prophage bait & junction panel design.

Preregistered design: ntm/v2/bait/DESIGN.md (task design-ntm-prophage).

All 31-mer operations use an exact 62-bit 2-bits-per-base encoding
(A=0, C=1, G=2, T=3) — no hashing collisions. The canonical form of a
k-mer is min(code, revcomp(code)), so masking is strand-symmetric.
"""

from __future__ import annotations

import hashlib
import re
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Preregistered constants (DESIGN.md Stage 3)
# ---------------------------------------------------------------------------

K = 31
BAIT_MIN_BP = 500
BAIT_MAX_BP = 2500

ENTROPY_MIN_BITS = 1.2      # low-complexity if k-mer base entropy below this
HOMOPOLYMER_RUN = 7          # low-complexity if k-mer overlaps >= this run
USABLE_KMER_MIN = 200        # min unmasked unique 31-mers to accept a bait
AMBIG_FRAC_MAX = 0.10
HOSTLIKE_FRAC_MAX = 0.10
LOWCX_FRAC_MAX = 0.30
CROSS_BAIT_DUP_FRAC = 0.50   # reject if > this fraction shared with kept bait
MARKER_COVER_FRAC = 0.50     # non-phage-marker domination threshold
CDS_COVER_FRAC_MIN = 0.30    # region gate needs this much CDS coverage

BASES = "ACGT"
_BASE_CODE = {b: i for i, b in enumerate(BASES)}  # A0 C1 G2 T3
_COMPLEMENT = {ord("A"): ord("T"), ord("T"): ord("A"),
               ord("C"): ord("G"), ord("G"): ord("C")}
for _b in "ACGT":
    _COMPLEMENT[ord(_b.lower())] = _COMPLEMENT[ord(_b)]

# Non-phage marker regex for region-level rejection (DESIGN.md Stage 3).
NONPHAGE_MARKER_RE = re.compile(
    r"transposase|insertion|plasmid|resistance|ribosomal|rRNA|tRNA",
    re.IGNORECASE,
)
# PHROG categories counted as housekeeping-like for domination purposes.
HOUSEKEEPING_CATEGORIES = {"DNA, RNA and nucleotide metabolism"}

# Module-gene targeting regexes, priority order (mirrors functional-QC module
# definitions; see ntm/v2/scripts/build_annotation_report_v2.py).
MODULE_PRIORITIES: List[Tuple[str, re.Pattern]] = [
    ("terminase", re.compile(r"terminase large", re.IGNORECASE)),
    ("portal", re.compile(r"portal", re.IGNORECASE)),
    ("capsid", re.compile(r"major (head|capsid)|capsid protein", re.IGNORECASE)),
    ("tail_tape_measure",
     re.compile(r"tape measure", re.IGNORECASE)),
    ("lysis", re.compile(r"endolysin|lysin|spanin|holin", re.IGNORECASE)),
    ("integrase", re.compile(r"integrase", re.IGNORECASE)),
]


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

def revcomp(seq: str) -> str:
    """Reverse complement (uppercase ACGT; other chars complemented via N)."""
    tbl = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(tbl)[::-1]


def sha256_seq(seq: str) -> str:
    """Stable sequence checksum: uppercase, no whitespace."""
    return hashlib.sha256(seq.upper().encode("ascii")).hexdigest()


def seq_to_arr(seq: str) -> np.ndarray:
    """uint8 codes 0..3 for ACGT (case-insensitive), 255 for anything else."""
    a = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8).copy()
    out = np.full(a.shape, 255, dtype=np.uint8)
    for base, code in _BASE_CODE.items():
        out[a == ord(base)] = code
    return out


def load_fasta(path: str) -> Dict[str, str]:
    """Load a (possibly bgzip/gzip) FASTA into {id: sequence} (id = first token)."""
    import gzip
    import shutil

    op = gzip.open if path.endswith((".gz", ".bgz")) else open
    with op(path, "rt") as fh:  # type: ignore[operator]
        return _parse_fasta(fh)


def _parse_fasta(fh: Iterable[str]) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    cur: Optional[str] = None
    buf: List[str] = []
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]
            buf = []
        elif cur is not None:
            buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def read_bgzf_records(path: str):
    """Yield (header_id, sequence) from a PanSN bgzip FASTA without loading all."""
    import gzip

    with gzip.open(path, "rt") as fh:
        cur = None
        buf: List[str] = []
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    yield cur, "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            elif cur is not None:
                buf.append(line.strip())
        if cur is not None:
            yield cur, "".join(buf)


# ---------------------------------------------------------------------------
# Exact 62-bit k-mer codes (numpy vectorised)
# ---------------------------------------------------------------------------

def kmer_codes(arr: np.ndarray, k: int = K) -> Tuple[np.ndarray, np.ndarray]:
    """All positional k-mer codes of a coded sequence.

    Returns (codes, valid): codes is uint64 (2 bits/base, first base in the
    high bits); valid[i] is False if the window contains any non-ACGT base.
    """
    n = arr.size
    if n < k:
        return np.empty(0, dtype=np.uint64), np.empty(0, dtype=bool)
    m = n - k + 1
    codes = np.zeros(m, dtype=np.uint64)
    bad = (arr == 255).astype(np.int64)
    csum = np.concatenate(([0], np.cumsum(bad)))
    valid = (csum[k:] - csum[:-k]) == 0
    for j in range(k):
        codes |= arr[j:j + m].astype(np.uint64) << np.uint64(2 * (k - 1 - j))
    return codes, valid


def canonical_codes(seq_or_arr, k: int = K) -> np.ndarray:
    """Canonical codes (min of forward and reverse-complement code) per position.

    Ambiguous windows get code 2**62 (an impossible ACGT code) so they never
    collide with real codes and never match host canonical codes.
    """
    arr = seq_to_arr(seq_or_arr) if isinstance(seq_or_arr, str) else seq_or_arr
    n = arr.size
    if n < k:
        return np.empty(0, dtype=np.uint64)
    fwd, valid = kmer_codes(arr, k)
    # reverse-complement sequence: complement + reverse
    comp = (arr != 255) * (3 - np.where(arr == 255, 0, arr)) + (arr == 255) * 255
    rc = comp[::-1]
    bwd, _ = kmer_codes(rc, k)
    bwd = bwd[::-1]  # align window i of fwd with its revcomp
    canon = np.minimum(fwd, bwd)
    canon[~valid] = np.uint64(1) << np.uint64(62)
    return canon


def count_kmers(canon: np.ndarray) -> Tuple[int, int]:
    """(total positional kmers, unique canonical kmers) for a bait."""
    if canon.size == 0:
        return 0, 0
    uniq = np.unique(canon).size
    return int(canon.size), int(uniq)


# ---------------------------------------------------------------------------
# Masking rules (per 31-mer)
# ---------------------------------------------------------------------------

def homopolymer_positions(arr: np.ndarray, run: int = HOMOPOLYMER_RUN) -> np.ndarray:
    """Boolean per sequence position: inside a >=`run` homopolymer (ACGT only)."""
    n = arr.size
    out = np.zeros(n, dtype=bool)
    if n == 0:
        return out
    valid = arr != 255
    change = np.empty(n, dtype=bool)
    change[0] = True
    change[1:] = (arr[1:] != arr[:-1]) | (~valid[1:]) | (~valid[:-1])
    idx = np.flatnonzero(change)
    if idx.size:
        lens = np.concatenate((np.diff(idx), [n - idx[-1]]))
        long_run = (lens >= run) & valid[idx]
        out = np.repeat(long_run, lens)  # vectorised segment expansion
    return out


def low_complexity_mask(seq_or_arr, k: int = K) -> np.ndarray:
    """Per-kmer boolean: homopolymer overlap OR entropy < ENTROPY_MIN_BITS."""
    arr = seq_to_arr(seq_or_arr) if isinstance(seq_or_arr, str) else seq_or_arr
    n = arr.size
    m = n - k + 1
    if m <= 0:
        return np.zeros(0, dtype=bool)
    hom = homopolymer_positions(arr)
    homwin = _window_sum(hom.astype(np.int64), k) > 0
    # base-composition entropy per kmer via cumulative counts
    ent = np.empty(m, dtype=np.float32)
    csums = {}
    for code in range(4):
        c = np.concatenate(([0], np.cumsum((arr == code).astype(np.int64))))
        csums[code] = c
    tot = np.stack([csums[c][k:] - csums[c][:-k] for c in range(4)])  # 4 x m
    frac = tot / k
    with np.errstate(divide="ignore", invalid="ignore"):
        logs = np.where(tot > 0, np.log2(np.where(tot > 0, frac, 1.0)), 0.0)
    ent = -(frac * logs).sum(axis=0).astype(np.float32)
    return homwin | (ent < ENTROPY_MIN_BITS)


def _window_sum(x: np.ndarray, k: int) -> np.ndarray:
    c = np.concatenate(([0], np.cumsum(x)))
    return c[k:] - c[:-k]


def internal_dup_mask(canon: np.ndarray) -> np.ndarray:
    """Per-kmer boolean: 2nd+ occurrence of a duplicated canonical kmer."""
    m = canon.size
    if m == 0:
        return np.zeros(0, dtype=bool)
    order = np.argsort(canon, kind="stable")
    s = canon[order]
    rep = np.zeros(m, dtype=bool)
    rep[1:] = s[1:] == s[:-1]
    dup = np.zeros(m, dtype=bool)
    # argsort is stable, so the first member of each sorted equal-run is the
    # leftmost original occurrence — mask every other member of the run.
    dup[order[rep]] = True
    return dup


def ambiguous_mask(seq_or_arr, k: int = K) -> np.ndarray:
    arr = seq_to_arr(seq_or_arr) if isinstance(seq_or_arr, str) else seq_or_arr
    _, valid = kmer_codes(arr, k)
    return ~valid


def scan_kmers_against(seq: str, query_canon_sorted: np.ndarray) -> int:
    """Number of `query` canonical kmers present in `seq` (host-like scan)."""
    canon = canonical_codes(seq)
    if canon.size == 0 or query_canon_sorted.size == 0:
        return 0
    hit = np.isin(canon, query_canon_sorted, assume_unique=False)
    return int(hit.sum())


def unique_intersection(a: np.ndarray, b_sorted: np.ndarray) -> np.ndarray:
    """Sorted unique values of `a` present in sorted array `b_sorted`."""
    if a.size == 0 or b_sorted.size == 0:
        return np.empty(0, dtype=a.dtype)
    return np.intersect1d(a, b_sorted, assume_unique=False)


# ---------------------------------------------------------------------------
# Window helpers (boundary coordinates)
# ---------------------------------------------------------------------------

def centered_window(genome_len: int, center: int, width: int) -> Tuple[int, int]:
    """1-based inclusive [start, end] window of `width` centered near `center`.

    Clamps to [1, genome_len] and shifts inward to preserve width when the
    genome is at least `width` long (never returns an off-genome window).
    """
    if genome_len < width:
        return 1, genome_len
    start = center - width // 2
    if start < 1:
        start = 1
    if start + width - 1 > genome_len:
        start = genome_len - width + 1
    return start, start + width - 1


def windows_overlap(a: Tuple[int, int], b: Tuple[int, int]) -> bool:
    return a[0] <= b[1] and b[0] <= a[1]


# ---------------------------------------------------------------------------
# First-order Markov negative control (composition-matched, gated)
# ---------------------------------------------------------------------------

def markov_control_sequence(seq: str, seed: int = 42) -> str:
    """First-order Markov negative control matching the source's base and
    dinucleotide composition (fitted; fixed-seed draw).

    Not an exact-count permutation: composition matches in expectation. The
    caller must gate the result on (a) negligible shared canonical 31-mers
    with the panel and (b) reported composition deviation (see DESIGN.md,
    control_shuffled_negative).
    """
    import random

    rng = random.Random(seed)
    s = seq.upper()
    n = len(s)
    bases = [c for c in s if c in "ACGT"]
    if not bases:
        return s
    mono = {b: bases.count(b) for b in "ACGT"}
    trans = {a: {b: 0 for b in "ACGT"} for a in "ACGT"}
    for a, b in zip(s, s[1:]):
        if a in trans and b in trans[a]:
            trans[a][b] += 1
    first = rng.choices("ACGT", weights=[mono[b] for b in "ACGT"])[0]
    out = [first]
    for _ in range(n - 1):
        row = trans[out[-1]]
        tot = sum(row.values())
        if tot == 0:  # dead end: resample from mono composition
            nxt = rng.choices("ACGT", weights=[mono[b] for b in "ACGT"])[0]
        else:
            nxt = rng.choices("ACGT", weights=[row[b] for b in "ACGT"])[0]
        out.append(nxt)
    return "".join(out)


def composition_deviation(a: str, b: str) -> Dict[str, float]:
    """L1 deviations of mono- and dinucleotide frequencies between two seqs."""
    def freqs(s, k):
        d: Dict[str, int] = {}
        tot = 0
        for i in range(len(s) - k + 1):
            key = s[i:i + k]
            d[key] = d.get(key, 0) + 1
            tot += 1
        return {key: v / tot for key, v in d.items()} if tot else {}

    out: Dict[str, float] = {}
    for k, tag in ((1, "mono"), (2, "dinuc")):
        fa, fb = freqs(a.upper(), k), freqs(b.upper(), k)
        keys = sorted(set(fa) | set(fb))
        out[f"{tag}_l1"] = sum(abs(fa.get(key, 0.0) - fb.get(key, 0.0))
                                for key in keys)
    return out


# ---------------------------------------------------------------------------
# Diversity-aware deterministic selection (farthest point, host-blind)
# ---------------------------------------------------------------------------

def greedy_farthest_point(
    candidates: Sequence[str],
    dist: Dict[Tuple[str, str], float],
    start: str,
    limit: int,
) -> List[str]:
    """Greedy farthest-point sampling over precomputed pairwise distances.

    Deterministic: the first candidate is `start`; each next pick maximises the
    minimum distance to all selected representatives (ties broken by candidate
    id). Returns at most `limit` ids in selection order.
    """
    if limit <= 0 or not candidates:
        return []
    selected = [start]
    remaining = [c for c in candidates if c != start]
    min_dist = {c: dist.get((c, start), dist.get((start, c), 1.0))
                for c in remaining}
    while remaining and len(selected) < limit:
        best = None
        best_d = -1.0
        for c in remaining:  # iteration order of `remaining` is stable (list)
            d = min_dist[c]
            if d > best_d + 1e-12 or (abs(d - best_d) <= 1e-12 and best is not None
                                      and c < best):
                best, best_d = c, d
        selected.append(best)
        remaining.remove(best)
        for c in remaining:
            d2 = dist.get((c, best), dist.get((best, c), 1.0))
            if d2 < min_dist[c]:
                min_dist[c] = d2
    return selected


def clade_medoid(members: Sequence[str],
                 sub_index: Dict[str, int],
                 D: np.ndarray) -> str:
    """Medoid member: minimises mean distance to co-members.

    NaN distances (mash triangles without shared hashes) count as 1.0 —
    documented in DESIGN.md. Ties break by lexicographic prophage id.
    """
    best_member, best_score = None, None
    for m in sorted(members):
        i = sub_index[m]
        vals = []
        for o in members:
            if o == m:
                continue
            j = sub_index[o]
            v = D[i, j]
            vals.append(1.0 if np.isnan(v) else float(v))
        score = float(np.mean(vals)) if vals else 0.0
        if best_score is None or score < best_score - 1e-12:
            best_member, best_score = m, score
    return best_member  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# Region-level annotation checks
# ---------------------------------------------------------------------------

def cds_composition_flags(cds_rows: Sequence[dict], start: int, end: int
                          ) -> Tuple[bool, float, float]:
    """Non-phage-marker domination test for window [start, end] (1-based incl).

    cds_rows: dicts with keys start, stop (1-based; either order per strand),
    annot, category, vfdb_hit, CARD_hit (pharokka merged table columns).
    Returns (dominated, marker_covered_frac_of_cds, cds_cover_frac_of_window).
    """
    win_len = end - start + 1
    cds_bases = 0
    marker_bases = 0
    for row in cds_rows:
        lo, hi = sorted((int(row["start"]), int(row["stop"])))
        ov = min(hi, end) - max(lo, start) + 1
        if ov <= 0:
            continue
        cds_bases += ov
        is_marker = bool(
            NONPHAGE_MARKER_RE.search(row.get("annot") or "")
            or (row.get("category") in HOUSEKEEPING_CATEGORIES)
            or ((row.get("vfdb_hit") or "None") not in ("None", "", "No hit"))
            or ((row.get("CARD_hit") or "None") not in ("None", "", "No hit"))
        )
        if is_marker:
            marker_bases += ov
    cds_cover_frac = cds_bases / win_len if win_len else 0.0
    marker_frac = (marker_bases / cds_bases) if cds_bases else 0.0
    dominated = (cds_cover_frac >= CDS_COVER_FRAC_MIN
                 and marker_frac >= MARKER_COVER_FRAC)
    return dominated, marker_frac, cds_cover_frac


def pick_module_genes(cds_rows: Sequence[dict]) -> List[dict]:
    """Key-module CDS per module, priority order, one per module.

    Rows are Pharokka merged-table dicts (start/stop/annot/category/phrog).
    A row qualifies for module `name` if its annot matches the module regex.
    Deterministic: first match by (priority, gene id sort).
    """
    picked: List[dict] = []
    for name, rx in MODULE_PRIORITIES:
        hits = [r for r in cds_rows
                if rx.search(r.get("annot") or "")]
        if not hits:
            continue
        hits.sort(key=lambda r: (min(int(r["start"]), int(r["stop"])),
                                 max(int(r["start"]), int(r["stop"])),
                                 str(r.get("gene") or "")))
        row = dict(hits[0])
        row["module"] = name
        picked.append(row)
    return picked


# ---------------------------------------------------------------------------
# FASTA writing
# ---------------------------------------------------------------------------

def write_fasta(records: Sequence[Tuple[str, str]], path: str,
                width: int = 80) -> None:
    with open(path, "w") as fh:
        for hid, seq in records:
            fh.write(f">{hid}\n")
            s = seq.upper()
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")
