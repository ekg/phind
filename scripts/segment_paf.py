#!/usr/bin/env python3
"""
segment_paf.py — CIGAR-aware PAF segmentation.

Splits each PAF alignment block into short chunks so that a downstream
partition step (impg partition -w 500) yields small intervals.

Rules (window W, default 500; max chunk span 2*W on BOTH query and target):
  * Walk the CIGAR op-by-op (query consumption: M/=/X/I/S; target
    consumption: M/=/X/D/N).
  * Any single op longer than W on its consuming axis is pre-split into
    pieces of size <= W (so chunk boundaries can land at ~W increments).
  * A chunk accumulates ops while both its query span and target span stay
    <= W; when the next op would push either span past W, the chunk is
    finalized and a new chunk starts with that op.
  * Consequence: every emitted chunk has query span <= W + last_piece and
    target span <= W + last_piece, i.e. <= 2*W = max_chunk_span. Chunks are
    typically ~W bp.
  * Chunks with no aligned content (zero query span or zero target span or
    no M/=/X ops) are dropped — they carry no usable alignment evidence.

Output: PAF with the same fields/order; columns 3-4 (qstart,qend), 8-9
(tstart,tend), 10 (nmatches), 11 (alnlen) and the cg:Z: CIGAR are rewritten
per chunk. All original non-cg tags are preserved.

Usage:
  segment_paf.py -i in.paf -o out.paf [-w 500] [--max-span 1000] [--gap 1]

A small inter-chunk gap (default 1 bp) prevents downstream range merging
(impg partition -d 0 merges only overlapping/touching ranges), keeping
partition intervals at ~window size. The gap drops that many aligned bp
from the PAF (0.2% at gap=1/window=500).

Reproducible command log is appended to commands.log in the output directory
if --log is given.
"""

import argparse
import re
import sys
from pathlib import Path

# PAF column indexes (0-based)
COL_QS = 2   # query start
COL_QE = 3   # query end
COL_TS = 7   # target start
COL_TE = 8   # target end
COL_NM = 9   # number of matching bases
COL_AL = 10  # alignment block length
COL_CG = 11  # (tag columns start here; cg:Z: found by scanning)

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

# CIGAR ops and their consumption
Q_OPS = {"M", "=", "X", "I", "S"}          # consume query
T_OPS = {"M", "=", "X", "D", "N"}          # consume target (reference)


def parse_cigar(cigar_str):
    """Return list of (op, length) in order."""
    out = []
    for m in CIGAR_RE.finditer(cigar_str):
        out.append((m.group(2), int(m.group(1))))
    return out


def op_consumption(op):
    return (op in Q_OPS, op in T_OPS)


def split_long_ops(ops, W):
    """Split ops whose consuming length exceeds W into <= W pieces.
    Returns a normalized stream of atomic pieces."""
    pieces = []
    for op, n in ops:
        q, t = op_consumption(op)
        axis = n if q else (n if t else 0)
        if axis <= W:
            pieces.append((op, n))
            continue
        # split into pieces of size <= W on the consuming axis
        # for M-like ops (consume both), pieces are equal length on both axes
        full, rem = divmod(n, W)
        for _ in range(full):
            pieces.append((op, W))
        if rem:
            pieces.append((op, rem))
    return pieces


def segment_block(cols, W, max_span, gap):
    """Segment one PAF block (list of 12+ strings). Returns list of block
    lines (each a list of strings), or [] if nothing usable."""
    qname = cols[0]
    qlen = int(cols[1])
    qstart = int(cols[2])
    qend = int(cols[3])
    strand = cols[4]
    tname = cols[5]
    tlen = int(cols[6])
    tstart = int(cols[7])
    tend = int(cols[8])

    tags = cols[12:]
    cigar = None
    for t in tags:
        if t.startswith("cg:Z:"):
            cigar = t[5:]
            break
    if not cigar:
        # no CIGAR: cannot segment safely; keep block as-is
        return [cols]

    ops = parse_cigar(cigar)
    pieces = split_long_ops(ops, W)

    # allwave emits '-' strand alignments with ASCENDING target coordinates
    # (strand records the query orientation, not a coordinate inversion), so
    # both query and target coordinates increase along the CIGAR walk.
    target_dir = 1

    qpos = qstart
    tpos = tstart
    chunks = []          # list of (qstart,qend,tstart,tend,pieces,nmatch,alnlen)
    cur_q = qpos
    cur_t = tpos
    cur_pieces = []
    cur_qs = 0
    cur_ts = 0
    cur_nm = 0
    cur_al = 0

    def finalize():
        nonlocal cur_q, cur_t, cur_pieces, cur_qs, cur_ts, cur_nm, cur_al
        if cur_pieces and cur_qs > 0 and cur_ts > 0 and cur_nm > 0:
            chunks.append((cur_q, cur_q + cur_qs, cur_t, cur_t + cur_ts,
                           list(cur_pieces), cur_nm, cur_al))
        cur_q = qpos
        cur_t = tpos
        cur_pieces = []
        cur_qs = 0
        cur_ts = 0
        cur_nm = 0
        cur_al = 0

    gap_remain = 0  # bp of alignment content to skip before next chunk

    for op, n in pieces:
        qadd = n if op in Q_OPS else 0
        tadd = n if op in T_OPS else 0
        # drop gap bp of alignment content after a chunk cut, so consecutive
        # chunks are not touching (keeps -d 0 partition ranges separate)
        if gap_remain:
            drop = min(gap_remain, n)
            gap_remain -= drop
            if drop == n:
                qpos += qadd
                tpos += tadd
                continue
            # partial drop inside this op: split it
            keep = n - drop
            # emulate: q/t consumed proportionally for M-like ops; for
            # single-axis ops only that axis advances
            qpos += drop if qadd else 0
            tpos += drop if tadd else 0
            qadd = keep if qadd else 0
            tadd = keep if tadd else 0
            n = keep
        if cur_qs > 0 and (cur_qs + qadd > W or cur_ts + tadd > W):
            finalize()
            gap_remain = gap
            if gap_remain:
                drop = min(gap_remain, n)
                gap_remain -= drop
                if drop == n:
                    qpos += qadd
                    tpos += tadd
                    continue
                keep = n - drop
                qpos += drop if qadd else 0
                tpos += drop if tadd else 0
                qadd = keep if qadd else 0
                tadd = keep if tadd else 0
                n = keep
        cur_pieces.append((op, n))
        cur_qs += qadd
        cur_ts += tadd
        if op in ("=", "M"):
            cur_nm += n
        if op in ("M", "=", "X", "I", "D"):
            cur_al += n
        qpos += qadd
        tpos += tadd * target_dir
    finalize()

    out = []
    for qs, qe, ts, te, cpieces, nm, al in chunks:
        # enforce the hard cap (defensive; should already hold)
        if qe - qs > max_span or abs(te - ts) > max_span:
            continue
        newcigar = "".join(f"{n}{op}" for op, n in cpieces)
        newcols = list(cols)
        newcols[COL_QS] = str(qs)
        newcols[COL_QE] = str(qe)
        newcols[COL_TS] = str(ts)
        newcols[COL_TE] = str(te)
        newcols[COL_NM] = str(nm)
        newcols[COL_AL] = str(al)
        newtags = [t for t in tags if not t.startswith("cg:Z:")]
        newtags.append(f"cg:Z:{newcigar}")
        newcols[12:] = newtags
        out.append(newcols)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-w", "--window", type=int, default=500)
    ap.add_argument("--max-span", type=int, default=None,
                    help="hard cap on chunk span (default 2*window)")
    ap.add_argument("--gap", type=int, default=1,
                    help="aligned bp dropped between consecutive chunks "
                         "(default 1; keeps ranges from merging at -d 0)")
    ap.add_argument("--log", default=None, help="commands.log path to append to")
    args = ap.parse_args()

    W = args.window
    max_span = args.max_span or 2 * W

    n_in = 0
    n_out = 0
    with open(args.input) as fin, open(args.output, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            n_in += 1
            for blk in segment_block(cols, W, max_span, args.gap):
                fout.write("\t".join(blk) + "\n")
                n_out += 1

    if args.log:
        with open(args.log, "a") as f:
            f.write(
                f"# segment_paf.py -i {args.input} -o {args.output} "
                f"-w {W} --max-span {max_span}  "
                f"[{n_in} blocks -> {n_out} chunks]\n"
            )

    print(f"segment_paf: {n_in} input blocks -> {n_out} chunks "
          f"(window {W}, max span {max_span})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
