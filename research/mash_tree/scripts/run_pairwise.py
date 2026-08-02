#!/usr/bin/env python3
"""run_pairwise.py — chunked all-vs-all mash dist -> binary float32 triangle.

Pipeline:
  1. Split the full FASTA into CHUNK_SIZE chunks (by sequence, preserving order).
  2. mash sketch -i each chunk (k=21, s=1000, seed 42 — must match the full sketch).
  3. For every chunk pair (i<=j): `mash dist -p THREADS ci.msh cj.msh | triangle_writer`.
     - i<j  -> offdiag mode (every output line is a unique upper-triangle pair)
     - i==j -> diag mode (keep only a<b)
  4. Writer appends float32 values at the canonical triangle offsets in one file.

Concurrency: NPROC parallel worker processes (mash dist does not scale past ~32
threads; process-level parallelism is what scales).

Output: full_prophages_mash.dist (little-endian float32, length n*(n-1)/2,
row-major upper triangle without diagonal). Side outputs: ids.txt (sequence
order), per-chunk logs, and a manifest of chunk indices.
"""
import argparse, os, subprocess, sys, time, random
from concurrent.futures import ThreadPoolExecutor

def run(cmd, logf, cwd):
    with open(logf, "ab") as lf:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=cwd)
        for blk in iter(lambda: p.stdout.read(1 << 16), b""):
            lf.write(blk)
        p.wait()
        return p.returncode

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ids", required=True)
    ap.add_argument("--chunk-size", type=int, default=5000)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--sketch-only", action="store_true")
    ap.add_argument("--dist-only", action="store_true")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    args.workdir = os.path.abspath(args.workdir)
    args.out = os.path.abspath(args.out)
    args.ids = os.path.abspath(args.ids)
    args.fasta = os.path.abspath(args.fasta)

    os.makedirs(args.workdir, exist_ok=True)
    os.makedirs(os.path.join(args.workdir, "chunks"), exist_ok=True)
    os.makedirs(os.path.join(args.workdir, "logs"), exist_ok=True)

    # 1. split fasta
    if not args.dist_only:
        t0 = time.time()
        chunks = {}   # ci -> [(header, [seqlines]), ...]
        buff = []
        with open(args.fasta) as f:
            for line in f:
                if line.startswith(">"):
                    if len(buff) >= args.chunk_size:
                        chunks[len(chunks)] = buff
                        buff = []
                    buff.append([line.strip()[1:].split()[0], []])
                else:
                    if buff:
                        buff[-1][1].append(line.strip())
        if buff:
            chunks[len(chunks)] = buff
        nchunks = len(chunks)
        # write chunk fastas + ids.txt
        with open(args.ids, "w") as f:
            for ci in sorted(chunks):
                with open(os.path.join(args.workdir, "chunks", f"chunk_{ci:03d}.fa"), "w") as cf:
                    for hid, seqlines in chunks[ci]:
                        cf.write(f">{hid}\n" + "".join(seqlines) + "\n")
                        f.write(hid + "\n")
        print(f"[split] {nchunks} chunks in {time.time()-t0:.1f}s", flush=True)
        # verify: ids length == expected
        with open(args.ids) as f:
            nids = sum(1 for _ in f)
        print(f"[split] total ids: {nids}", flush=True)

        # 2. sketch chunks
        t0 = time.time()
        for ci in sorted(chunks):
            msh = os.path.join(args.workdir, "chunks", f"chunk_{ci:03d}.msh")
            if os.path.exists(msh):
                continue
            r = run(["mash", "sketch", "-i", "-k", "21", "-s", "1000",
                     "-p", str(args.threads), "-o", msh,
                     os.path.join(args.workdir, "chunks", f"chunk_{ci:03d}.fa")],
                    os.path.join(args.workdir, "logs", f"sketch_{ci:03d}.log"),
                    args.workdir)
            if r != 0:
                sys.exit(f"sketch chunk {ci} failed (see log)")
        print(f"[sketch] {nchunks} chunk sketches in {time.time()-t0:.1f}s", flush=True)

    if args.sketch_only:
        return

    # 3. pairwise -> per-job part files (separate files => no shared-file
    #    write contention, which we measured as ~8x pathological slowdown)
    t0 = time.time()
    nchunks = len([x for x in os.listdir(os.path.join(args.workdir, "chunks"))
                   if x.endswith(".msh")])
    with open(args.ids) as _f:
        n_total = sum(1 for _ in _f)
    sizes = [min(args.chunk_size, n_total - c * args.chunk_size) for c in range(nchunks)]
    starts = [sum(sizes[:c]) for c in range(nchunks)]
    parts_dir = os.path.join(args.workdir, "parts")
    os.makedirs(parts_dir, exist_ok=True)
    jobs = []
    for i in range(nchunks):
        for j in range(i, nchunks):
            jobs.append((i, j))
    random.Random(1234).shuffle(jobs)
    writer = os.path.join(os.path.dirname(os.path.abspath(__file__)), "part_writer")

    def do_job(pair):
        i, j = pair
        out = os.path.join(args.workdir, "logs", f"part_{i:03d}_{j:03d}.log")
        msh_i = os.path.join(args.workdir, "chunks", f"chunk_{i:03d}.msh")
        msh_j = os.path.join(args.workdir, "chunks", f"chunk_{j:03d}.msh")
        part = os.path.join(parts_dir, f"part_{i:03d}_{j:03d}.bin")
        with open(part, "wb") as pf:
            if i < j:
                # ref = later chunk (b-side), query = earlier (a-side) so the
                # value stream is (a-major, b-minor) = triangle row order
                p1 = subprocess.Popen(["mash", "dist", "-p", str(args.threads),
                                       msh_j, msh_i], stdout=subprocess.PIPE)
                p2 = subprocess.Popen([writer, args.ids, "offdiag", "0", "0"],
                                      stdin=p1.stdout, stdout=pf,
                                      stderr=open(out, "ab"))
            else:
                p1 = subprocess.Popen(["mash", "dist", "-p", str(args.threads),
                                       msh_i, msh_i], stdout=subprocess.PIPE)
                p2 = subprocess.Popen([writer, args.ids, "diag",
                                       str(starts[i]), str(sizes[i])],
                                      stdin=p1.stdout, stdout=pf,
                                      stderr=open(out, "ab"))
            p1.stdout.close()
            rc = p2.wait()
            p1.wait()
        if rc != 0 or p1.returncode != 0:
            print(f"  FAILED chunk pair {i},{j} (mash rc={p1.returncode} writer rc={rc})",
                  flush=True)
        return (i, j, p1.returncode, rc)

    with ThreadPoolExecutor(max_workers=args.procs) as ex:
        for k, res in enumerate(ex.map(do_job, jobs), 1):
            if k % 27 == 0 or k == len(jobs):
                elapsed = time.time() - t0
                print(f"[dist] {k}/{len(jobs)} jobs done in {elapsed:.0f}s "
                      f"({elapsed/k*(len(jobs)-k)/60:.1f} min remaining)", flush=True)
    print(f"[dist] all {len(jobs)} chunk-pair jobs in {time.time()-t0:.1f}s", flush=True)

    # 4. merge part files into the final triangle
    t0 = time.time()
    m = subprocess.run([sys.executable, os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "merge_parts.py"),
        "--workdir", args.workdir, "--out", args.out, "--ids", args.ids,
        "--chunk-size", str(args.chunk_size)],
        capture_output=True, text=True)
    print(m.stdout, flush=True)
    print(m.stderr, flush=True)
    print(f"[merge] finished in {time.time()-t0:.1f}s", flush=True)

if __name__ == "__main__":
    main()
