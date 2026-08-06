#!/usr/bin/env python3
"""host_pipeline.py — sketch pre-built chunk fastas, compute float32 triangle,
merge. Mirrors the prophage run_pairwise.py but starts from chunks produced by
build_chunks.py (one concatenated genome per record) and uses k=21 s=10000.

Steps:
  1. mash sketch -i each chunk (k=21, s=SKETCH_SIZE, seed 42)
  2. pairwise chunk dist -> per-job part files (part_writer: offdiag/diag)
  3. merge parts -> OUT (float32 upper triangle, n*(n-1)/2 values)
  4. also write full hosts.msh (mash paste of chunk sketches)
"""
import argparse, os, subprocess, sys, time, random
from concurrent.futures import ThreadPoolExecutor

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ids", required=True)
    ap.add_argument("--msh", required=True, help="output full hosts.msh")
    ap.add_argument("--chunk-size", type=int, default=5000)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--sketch-size", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--mash", default="mash")
    ap.add_argument("--skip-sketch", action="store_true")
    ap.add_argument("--skip-dist", action="store_true")
    args = ap.parse_args()

    workdir = os.path.abspath(args.workdir)
    chunks_dir = os.path.join(workdir, "chunks")
    parts_dir = os.path.join(workdir, "parts")
    logs_dir = os.path.join(workdir, "logs")
    os.makedirs(parts_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    with open(args.ids) as f:
        n_total = sum(1 for _ in f)
    chunk_fas = sorted([x for x in os.listdir(chunks_dir) if x.endswith(".fa")])
    nchunks = len(chunk_fas)
    sizes = [min(args.chunk_size, n_total - c * args.chunk_size)
             for c in range(nchunks)]
    starts = [sum(sizes[:c]) for c in range(nchunks)]
    tri_len = n_total * (n_total - 1) // 2
    print(f"[pipe] n={n_total} chunks={nchunks} tri_len={tri_len} "
          f"({tri_len*4/1e9:.2f} GB)", flush=True)

    # 1. sketch chunks
    if not args.skip_sketch:
        t0 = time.time()
        for ci in range(nchunks):
            fa = os.path.join(chunks_dir, f"chunk_{ci:03d}.fa")
            msh = os.path.join(chunks_dir, f"chunk_{ci:03d}.msh")
            if os.path.exists(msh):
                continue
            r = subprocess.run([args.mash, "sketch", "-i", "-k", "21",
                                "-s", str(args.sketch_size), "-p", str(args.threads),
                                "-o", msh, fa],
                               capture_output=True, text=True)
            if r.returncode != 0:
                sys.exit(f"sketch chunk {ci} failed: {r.stderr}")
        print(f"[pipe] sketched {nchunks} chunks in {time.time()-t0:.1f}s", flush=True)

    # 2. pairwise -> per-job part files (separate files => no shared write contention)
    jobs = [(i, j) for i in range(nchunks) for j in range(i, nchunks)]
    random.Random(1234).shuffle(jobs)
    writer = os.path.join(os.path.dirname(os.path.abspath(__file__)), "part_writer")

    def do_job(pair):
        i, j = pair
        logf = os.path.join(logs_dir, f"part_{i:03d}_{j:03d}.log")
        msh_i = os.path.join(chunks_dir, f"chunk_{i:03d}.msh")
        msh_j = os.path.join(chunks_dir, f"chunk_{j:03d}.msh")
        part = os.path.join(parts_dir, f"part_{i:03d}_{j:03d}.bin")
        with open(part, "wb") as pf:
            if i < j:
                p1 = subprocess.Popen([args.mash, "dist", "-p", str(args.threads),
                                       msh_j, msh_i], stdout=subprocess.PIPE)
                p2 = subprocess.Popen([writer, args.ids, "offdiag", "0", "0"],
                                      stdin=p1.stdout, stdout=pf,
                                      stderr=open(logf, "ab"))
            else:
                p1 = subprocess.Popen([args.mash, "dist", "-p", str(args.threads),
                                       msh_i, msh_i], stdout=subprocess.PIPE)
                p2 = subprocess.Popen([writer, args.ids, "diag",
                                       str(starts[i]), str(sizes[i])],
                                      stdin=p1.stdout, stdout=pf,
                                      stderr=open(logf, "ab"))
            p1.stdout.close()
            rc = p2.wait(); p1.wait()
        if rc != 0 or p1.returncode != 0:
            print(f"  FAILED chunk pair {i},{j} (mash={p1.returncode} writer={rc})", flush=True)
        return (i, j, p1.returncode, rc)

    if not args.skip_dist:
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=args.procs) as ex:
            for k, res in enumerate(ex.map(do_job, jobs), 1):
                if k % max(1, len(jobs)//20) == 0 or k == len(jobs):
                    el = time.time() - t0
                    print(f"[pipe] {k}/{len(jobs)} jobs in {el:.0f}s "
                          f"({el/k*(len(jobs)-k)/60:.1f} min remaining)", flush=True)
        print(f"[pipe] all {len(jobs)} chunk-pair jobs in {time.time()-t0:.1f}s", flush=True)
        # 3. merge
        t0 = time.time()
        m = subprocess.run([sys.executable,
                            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "merge_parts.py"),
                            "--workdir", workdir, "--out", args.out, "--ids", args.ids,
                            "--chunk-size", str(args.chunk_size)],
                           capture_output=True, text=True)
        print(m.stdout, flush=True); print(m.stderr, flush=True)
        print(f"[pipe] merge finished in {time.time()-t0:.1f}s", flush=True)

    # 4. full hosts.msh via mash paste (concatenate chunk sketches)
    if not os.path.exists(args.msh):
        chunk_mshs = [os.path.join(chunks_dir, f"chunk_{c:03d}.msh")
                      for c in range(nchunks)]
        with open(os.path.join(workdir, "chunks", "list.txt"), "w") as f:
            for p in chunk_mshs:
                f.write(p + "\n")
        r = subprocess.run([args.mash, "paste", "-l", args.msh,
                            os.path.join(workdir, "chunks", "list.txt")],
                           capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"mash paste failed: {r.stderr}")
        print(f"[pipe] wrote {args.msh}", flush=True)

if __name__ == "__main__":
    main()