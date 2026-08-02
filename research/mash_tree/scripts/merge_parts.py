#!/usr/bin/env python3
"""merge_parts.py — assemble per-job part files into the float32 triangle.

Part files (data/parts/part_III_JJJ.bin) each hold float32 distances in
triangle row order for chunk pair (III, JJJ) (III <= JJJ):
  - offdiag (III<JJJ): (a-major, b-minor) for a in chunk_III, b in chunk_JJJ
  - diag    (III==JJJ): upper triangle row-major within the chunk
The merge writes them into full_prophages_mash.dist at the canonical offsets,
row by row (contiguous ~20 KB writes).

Reads: ids.txt, data/parts/*.bin
Writes: full_prophages_mash.dist (little-endian float32, n*(n-1)/2 values)
"""
import json, os, sys, time
import numpy as np

def off(a, b, n):
    # upper-triangle offset for a<b (row-major, no diagonal)
    return a * (2 * n - a - 1) // 2 + (b - a - 1)

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ids", required=True)
    ap.add_argument("--chunk-size", type=int, default=5000)
    a = ap.parse_args()
    chunk_size = a.chunk_size
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ids_path = a.ids
    parts_dir = os.path.join(a.workdir, "parts")
    out_path = a.out

    with open(ids_path) as f:
        ids = [l.strip() for l in f if l.strip()]
    n = len(ids)
    nchunks = (n + chunk_size - 1) // chunk_size
    sizes = [min(chunk_size, n - c * chunk_size) for c in range(nchunks)]
    starts = [sum(sizes[:c]) for c in range(nchunks)]
    tri_len = n * (n - 1) // 2
    print(f"[merge] n={n} chunks={nchunks} tri_len={tri_len} ({tri_len*4/1e9:.1f} GB)", flush=True)

    tri = None
    with open(out_path, "wb") as f:
        f.truncate(tri_len * 4)
    tri = np.memmap(out_path, dtype="<f4", mode="r+", shape=(tri_len,))
    t0 = time.time()
    total_written = 0
    for c in range(nchunks):
        for d in range(c, nchunks):
            pf = os.path.join(parts_dir, f"part_{c:03d}_{d:03d}.bin")
            if not os.path.exists(pf):
                sys.exit(f"[merge] missing part file {pf}")
            vals = np.fromfile(pf, dtype="<f4")
            if c < d:
                mi, mj = sizes[c], sizes[d]
                expect = mi * mj
                if len(vals) != expect:
                    sys.exit(f"[merge] {pf}: expected {expect} values, got {len(vals)}")
                sj = starts[d]
                for x in range(mi):
                    a = starts[c] + x
                    row_off = off(a, sj, n)
                    tri[row_off:row_off + mj] = vals[x * mj:(x + 1) * mj]
                total_written += expect
            else:
                m = sizes[c]
                expect = m * (m - 1) // 2
                if len(vals) != expect:
                    sys.exit(f"[merge] {pf}: expected {expect} values, got {len(vals)}")
                pos = 0
                for x in range(m):
                    a = starts[c] + x
                    row_off = off(a, a + 1, n)
                    cnt = m - 1 - x
                    tri[row_off:row_off + cnt] = vals[pos:pos + cnt]
                    pos += cnt
                total_written += expect
            if (c * nchunks + d) % 54 == 0:
                print(f"[merge] {c* nchunks + d}/{nchunks*nchunks} block pairs "
                      f"({(c*nchunks+d)*100//(nchunks*nchunks)}%) "
                      f"{time.time()-t0:.0f}s", flush=True)
    tri.flush()
    del tri
    size_actual = os.path.getsize(out_path)
    size_expected = tri_len * 4
    print(f"[merge] wrote {total_written} values in {time.time()-t0:.0f}s", flush=True)
    print(f"[merge] file size {size_actual} expected {size_expected} "
          f"match={size_actual == size_expected}", flush=True)

if __name__ == "__main__":
    main()
