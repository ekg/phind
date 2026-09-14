#!/usr/bin/env python3
"""
rederive_sorted_clades.py — re-derive corrected NTM v1/v2 tight-clade counts.

Background (found during ntm-v3-prophage; see ntm/v3/mash_clades_report.md
"Mechanics note" and task verify-v2-clade-order): scripts/build_tight_clades.py
sorts community members with sorted(set(members)) and read_community_matrix only
fills D[i][j] when the member-list order agrees with triangle row order. v1/v2
wrote ids.txt in FASTA header order, so the string-sorted member list disagreed
with the triangle row order for 35.2% of v2's within-community pairs (measured:
12,730,123 of 36,137,751 unique pairs NaN in v2/clades/0/distances.npz). NaN
compares False against every threshold, so those pairs silently behaved as
maximally distant: missed joins -> fragmented (inflated) clade counts.

v3 fixed this upstream by writing ids.txt pre-sorted. For v1/v2 the triangle
bytes are frozen, and a bare sorted ids.txt would MISINDEX them (triangle row
k is the k-th FASTA-order id, not the k-th sorted id). This driver therefore
re-derives the counterfactual "v2 as if ids.txt had been written pre-sorted":

  1. read the frozen ids.txt (== triangle row order) and the frozen triangle;
  2. write a NEW mash_clades_sorted/ with ids.txt sorted(), labels.csv in the
     same order, and the triangle PERMUTED into sorted-id row order — a
     lossless reordering of the same float32 values (MASH pairwise distances
     do not depend on sketch input order, so this is exactly the triangle mash
     would have produced over a pre-sorted FASTA; v3's approach);
  3. run scripts/build_tight_clades.py UNMODIFIED with the original v1/v2
     flags (--threshold 0.25 --max-size 100 --communities 0) into a NEW
     clades_sorted/ outdir;
  4. validate:
       - 0 NaN in the new community matrix;
       - exact (bitwise) agreement with the OLD matrix wherever the old one
         was filled — proves the permuted triangle reuses the frozen bytes;
       - TSV spot checks against the frozen mash dist text output;
       - clade invariants: sum of sizes == n, ids unique, size <= max-size,
         every non-singleton internal median <= threshold (recomputed
         independently from the matrix, not trusted from the run's json);
  5. emit clade_summary.tsv / alignable_clades.tsv / singletons.tsv (same
     format as the original v2 driver) + sorted_vs_published.json comparing
     against the published (fragmented) numbers.

The frozen v1/v2 outputs are only ever READ (new outdirs only).

Usage:
  python3 ntm/scripts/rederive_sorted_clades.py --version v2
  python3 ntm/scripts/rederive_sorted_clades.py --version v1
  python3 ntm/scripts/rederive_sorted_clades.py --selftest
"""
import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np

NVME = "/mnt/nvme3n1/erikg/phind-genome-work/ntm"
THRESHOLD = 0.25
MAX_SIZE = 100
COMMUNITY = "0"


def repo_root():
    return os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))


def sha256(path, buf=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def triangle_offset(a, b, n):
    return a * (2 * n - a - 1) // 2 + (b - a - 1)


def selftest():
    """Tiny synthetic check of the permutation math: 40 ids in scrambled
    order, random triangle; permute; brute-force compare every pair."""
    rng = np.random.default_rng(7)
    n = 40
    base = [f"id_{i:03d}" for i in range(n)]
    scrambled = list(base)
    rng.shuffle(scrambled)
    # triangle in scrambled order
    old = rng.random(n * (n - 1) // 2).astype(np.float32)
    pos = {s: i for i, s in enumerate(scrambled)}
    order = sorted(range(n), key=lambda i: scrambled[i])  # orig idx, sorted ids
    new = np.empty_like(old)
    for k in range(n):
        a = order[k]
        oj = np.array(order[k + 1:], dtype=np.int64)
        if not len(oj):
            continue
        lo = np.minimum(a, oj)
        hi = np.maximum(a, oj)
        src = lo * (2 * n - lo - 1) // 2 + (hi - lo - 1)
        s = triangle_offset(k, k + 1, n)
        new[s: s + len(oj)] = old[src]
    # brute force: value for sorted pair (x, y) must equal old triangle
    # value for (pos[x], pos[y]) with x < y in sorted order
    for x in range(n):
        for y in range(x + 1, n):
            a, b = pos[base[x]], pos[base[y]]
            lo, hi = min(a, b), max(a, b)
            v_old = old[triangle_offset(lo, hi, n)]
            v_new = new[triangle_offset(x, y, n)]
            assert v_old == v_new, (x, y, v_old, v_new)
    print("selftest PASS (40x40 permutation exact)")


def permute_triangle(old_tri, new_tri, ids):
    n = len(ids)
    npairs = n * (n - 1) // 2
    assert os.path.getsize(old_tri) == npairs * 4, "old triangle size mismatch"
    order = sorted(range(n), key=lambda i: ids[i])   # orig row idx, sorted ids
    mm = np.memmap(old_tri, dtype="<f4", mode="r", shape=(npairs,))
    t0 = time.time()
    with open(new_tri + ".tmp", "wb") as f:
        for k in range(n):
            a = order[k]
            oj = np.array(order[k + 1:], dtype=np.int64)
            if not len(oj):
                continue
            lo = np.minimum(a, oj)
            hi = np.maximum(a, oj)
            src = lo * (2 * n - lo - 1) // 2 + (hi - lo - 1)
            f.write(mm[src].tobytes())
    os.replace(new_tri + ".tmp", new_tri)
    assert os.path.getsize(new_tri) == npairs * 4, "new triangle size mismatch"
    return order, time.time() - t0


def tsv_spotcheck(dist_tsv, ids, new_tri, seed=0, k=5):
    """Independent source: pick k id pairs, find their distance in the frozen
    mash dist TEXT output (one awk pass), compare to the permuted triangle."""
    n = len(ids)
    rng = np.random.default_rng(seed)
    pairs = []
    while len(pairs) < k:
        i, j = sorted(rng.integers(0, n, 2))
        if i < j:
            pairs.append((int(i), int(j)))
    conds = []
    for i, j in pairs:
        a, b = ids[i], ids[j]
        conds.append(f'($1=="{a}" && $2=="{b}")')
        conds.append(f'($1=="{b}" && $2=="{a}")')
    awk = 'awk \'BEGIN{FS="\\t"} ' + " || ".join(conds) + \
          ' {print $1"\\t"$2"\\t"$3}\''
    env = dict(os.environ, LC_ALL="C")
    out = subprocess.run(["bash", "-c", awk + " " + dist_tsv], env=env,
                         capture_output=True, text=True, check=True).stdout
    pairset = {frozenset((ids[i], ids[j])) for i, j in pairs}
    found = {}
    for line in out.strip().splitlines():
        a, b, d = line.split("\t")
        key = frozenset((a, b))
        if key in pairset:
            found[key] = float(d)
    assert len(found) >= k, f"spot check: only {len(found)}/{k} pairs found in TSV"
    mm = np.memmap(new_tri, dtype="<f4", mode="r", shape=(n * (n - 1) // 2,))
    worst = 0.0
    for i, j in pairs:
        key = frozenset((ids[i], ids[j]))
        d_tsv = found[key]
        d_new = float(mm[triangle_offset(i, j, n)])
        worst = max(worst, abs(d_tsv - d_new))
    assert worst <= 2e-7, f"spot check mismatch: worst |delta| {worst}"
    return len(found), worst


def frozen_files(src_mash, old_clades):
    """Frozen inputs this driver reads (never writes): hashed before and
    after the run to prove they are untouched."""
    return [
        f"{src_mash}/ids.txt",
        f"{src_mash}/labels.csv",
        f"{src_mash}/prophages_mash.dist",
        f"{src_mash}/prophages.dist.tsv",
        f"{old_clades}/tight_clades_summary.json",
        f"{old_clades}/0/tight_clades.json",
        f"{old_clades}/0/clade_similarity.json",
        f"{old_clades}/0/members.json",
        f"{old_clades}/0/distances.npz",
    ]


def hash_files(paths):
    """sha256 per file (dist.tsv is multi-GB; hashing it takes ~20-30 s —
    acceptable for a one-shot verification run)."""
    out = {}
    for p in paths:
        if os.path.exists(p):
            out[p] = sha256(p)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--version", choices=["v1", "v2"], default="v2")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--resume", action="store_true",
                    help="reuse an existing mash_clades_sorted/ + clades_sorted/ "
                         "(e.g. after a crashed validation step); only ever "
                         "writes inside those new outdirs")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return

    v = args.version
    src_mash = f"{NVME}/{v}/mash_clades"
    old_clades = f"{NVME}/{v}/clades"
    new_mash = f"{NVME}/{v}/mash_clades_sorted"
    new_clades = f"{NVME}/{v}/clades_sorted"
    resume = args.resume
    if not resume:
        for p in (new_mash, new_clades):
            if os.path.exists(p):
                raise SystemExit(f"REFUSING to overwrite existing {p}")
    build_tight = os.path.join(repo_root(), "scripts", "build_tight_clades.py")
    assert os.path.exists(build_tight), build_tight

    # frozen-input hashes BEFORE any work (also asserted again at the end)
    frozen_paths = frozen_files(src_mash, old_clades)
    frozen_before = hash_files(frozen_paths)
    print(f"[{v}] hashed {len(frozen_before)} frozen inputs before work",
          flush=True)

    # ---- 1. frozen inputs -------------------------------------------------
    ids = [l.strip() for l in open(f"{src_mash}/ids.txt")]
    n = len(ids)
    assert len(set(ids)) == n, "duplicate ids in frozen ids.txt"
    assert ids != sorted(ids), "ids.txt already sorted?!"
    old_tri = f"{src_mash}/prophages_mash.dist"
    labels = {}
    with open(f"{src_mash}/labels.csv") as f:
        for row in csv.DictReader(f):
            labels[row["sequence"]] = row["community"]
    assert set(labels) == set(ids), "labels.csv does not cover ids exactly"
    print(f"[{v}] frozen inputs: {n} ids (unsorted, FASTA order), "
          f"triangle {os.path.getsize(old_tri)} bytes", flush=True)

    # ---- 2. sorted copies + permuted triangle -----------------------------
    new_tri = f"{new_mash}/prophages_mash_sorted.dist"
    if resume and os.path.exists(new_tri) and \
            os.path.exists(f"{new_mash}/ids.txt") and \
            os.path.exists(f"{new_mash}/labels.csv"):
        sorted_ids = [l.strip() for l in open(f"{new_mash}/ids.txt")]
        assert sorted_ids == sorted(ids), "resumed ids.txt is not sorted(ids)"
        assert os.path.getsize(new_tri) == n * (n - 1) // 2 * 4, \
            "resumed triangle size mismatch"
        order = sorted(range(n), key=lambda i: ids[i])
        print(f"[{v}] RESUME: reusing {new_mash} (verified vs frozen ids.txt)",
              flush=True)
    else:
        os.makedirs(new_mash, exist_ok=True)
        sorted_ids = sorted(ids)
        with open(f"{new_mash}/ids.txt", "w") as f:
            f.write("\n".join(sorted_ids) + "\n")
        with open(f"{new_mash}/labels.csv", "w") as f:
            f.write("sequence,community\n")
            for s in sorted_ids:
                f.write(f"{s},{labels[s]}\n")
        order, t_perm = permute_triangle(old_tri, new_tri, ids)
        print(f"[{v}] permuted triangle -> {new_tri} ({t_perm:.1f}s)", flush=True)

    provenance = {
        "frozen_ids_sha256": sha256(f"{src_mash}/ids.txt"),
        "frozen_triangle_sha256": sha256(old_tri),
        "sorted_triangle_sha256": sha256(new_tri),
        "frozen_triangle_bytes": os.path.getsize(old_tri),
        "n": n,
    }
    with open(f"{new_mash}/provenance.json", "w") as f:
        json.dump(provenance, f, indent=1)

    # ---- 3. run build_tight_clades.py unmodified --------------------------
    if resume and os.path.exists(f"{new_clades}/0/tight_clades.json"):
        print(f"[{v}] RESUME: reusing existing build output in {new_clades}",
              flush=True)
    else:
        os.makedirs(new_clades, exist_ok=True)
        cmd = [sys.executable, build_tight,
               "--threshold", str(THRESHOLD), "--max-size", str(MAX_SIZE),
               "--communities", COMMUNITY,
               "--outdir", new_clades,
               "--ids-file", f"{new_mash}/ids.txt",
               "--triangle", new_tri,
               "--labels-csv", f"{new_mash}/labels.csv"]
        with open(f"{new_clades}/commands.log", "a") as f:
            f.write("# " + " ".join(cmd) + "\n")
        print(f"[{v}] running build_tight_clades.py unmodified ...", flush=True)
        t0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        with open(f"{new_clades}/build_tight_clades.stdout.log", "w") as f:
            f.write(r.stdout + "\n--- stderr ---\n" + r.stderr)
        if r.returncode != 0:
            raise SystemExit(f"build_tight_clades.py FAILED rc={r.returncode}; "
                             f"see {new_clades}/build_tight_clades.stdout.log")
        print(f"[{v}] build_tight_clades.py done in {time.time()-t0:.0f}s",
              flush=True)

    # ---- 4. validation -----------------------------------------------------
    checks = {}
    d_new = np.load(f"{new_clades}/0/distances.npz", allow_pickle=True)
    D_new, mem_new = d_new["D"], [str(x) for x in d_new["members"]]
    d_old = np.load(f"{old_clades}/0/distances.npz", allow_pickle=True)
    D_old, mem_old = d_old["D"], [str(x) for x in d_old["members"]]
    offdiag = ~np.eye(n, dtype=bool)
    nan_old = int(np.isnan(D_old[offdiag]).sum()) // 2
    nan_new = int(np.isnan(D_new[offdiag]).sum()) // 2
    total_pairs = n * (n - 1) // 2
    assert mem_old == mem_new == sorted_ids, "member order mismatch"
    assert nan_new == 0, f"new matrix has {nan_new} NaN pairs"
    filled = ~np.isnan(D_old)
    mism = int((D_new[filled] != D_old[filled]).sum())
    assert mism == 0, f"{mism} filled entries disagree with old matrix"
    checks["old_nan_pairs"] = nan_old
    checks["old_nan_fraction"] = nan_old / total_pairs
    checks["new_nan_pairs"] = nan_new
    checks["filled_entries_exact_match"] = mism == 0
    print(f"[{v}] matrix: old NaN {nan_old}/{total_pairs} "
          f"({100*nan_old/total_pairs:.2f}%), new NaN {nan_new}, "
          f"filled-entry exact match: {mism == 0}", flush=True)

    n_found, worst = tsv_spotcheck(f"{src_mash}/prophages.dist.tsv",
                                   sorted_ids, new_tri)
    checks["tsv_spotcheck_pairs"] = n_found
    checks["tsv_spotcheck_worst_abs_delta"] = worst
    print(f"[{v}] TSV spot check: {n_found} pairs, worst |delta| {worst:.1e}",
          flush=True)

    tc = json.load(open(f"{new_clades}/0/tight_clades.json"))
    pos = {s: i for i, s in enumerate(sorted_ids)}
    sizes = [len(m) for m in tc.values()]
    allmem = [m for ms in tc.values() for m in ms]
    assert sum(sizes) == n, f"sum of clade sizes {sum(sizes)} != {n}"
    assert len(set(allmem)) == n, "duplicate ids across clades"
    assert max(sizes) <= MAX_SIZE, "clade exceeds max-size"
    worst_med = 0.0
    sim_med_mism = 0
    sim = json.load(open(f"{new_clades}/0/clade_similarity.json"))["per_clade"]
    for cid, ms in tc.items():
        if len(ms) < 2:
            continue
        idx = np.array([pos[s] for s in ms])
        vals = D_new[np.ix_(idx, idx)][np.triu_indices(len(ms), 1)]
        med = float(np.median(vals))
        worst_med = max(worst_med, med)
        if sim[cid]["median"] is None or abs(sim[cid]["median"] - med) > 1e-9:
            sim_med_mism += 1
    assert worst_med <= THRESHOLD + 1e-6, \
        f"non-singleton median {worst_med} > {THRESHOLD}"
    assert sim_med_mism == 0, "clade_similarity.json medians disagree"
    checks["sum_sizes_eq_n"] = True
    checks["ids_unique"] = True
    checks["max_clade_size"] = max(sizes)
    checks["worst_non_singleton_median"] = worst_med
    checks["similarity_json_medians_recomputed_ok"] = True
    print(f"[{v}] clades: sum sizes == {n}, unique, max size {max(sizes)}, "
          f"worst non-singleton median {worst_med:.6f} <= {THRESHOLD}",
          flush=True)

    # ---- 5. summaries + comparison ----------------------------------------
    rows = []
    for cid, ms in tc.items():
        s = sim.get(cid, {})
        rows.append((cid, len(ms), s.get("median"), s.get("min"), s.get("max")))
    rows.sort()

    def g(x):
        return "" if x is None else f"{x:.6f}"

    with open(f"{new_clades}/clade_summary.tsv", "w") as f:
        f.write("clade_id\tn_members\tmedian_mash\tmin_mash\tmax_mash\n")
        for cid, nm, med, mn, mx in rows:
            f.write(f"{cid}\t{nm}\t{g(med)}\t{g(mn)}\t{g(mx)}\n")
    with open(f"{new_clades}/alignable_clades.tsv", "w") as f:
        f.write("clade_id\tn_members\tmedian_mash\tmin_mash\tmax_mash\n")
        for cid, nm, med, mn, mx in rows:
            if nm >= 2:
                f.write(f"{cid}\t{nm}\t{g(med)}\t{g(mn)}\t{g(mx)}\n")
    with open(f"{new_clades}/singletons.tsv", "w") as f:
        f.write("clade_id\tmember_id\n")
        for cid, nm, med, mn, mx in rows:
            if nm == 1:
                f.write(f"{cid}\t{tc[cid][0]}\n")

    # ---- frozen-outputs-untouched gate: re-hash and compare -------------
    frozen_after = hash_files(frozen_paths)
    assert frozen_before == frozen_after, \
        "FROZEN INPUTS CHANGED during the run: " + \
        str({k for k in frozen_before if frozen_before[k] != frozen_after.get(k)})
    checks["frozen_inputs_hashed_before"] = len(frozen_before)
    checks["frozen_inputs_unchanged_before_vs_after"] = True
    print(f"[{v}] frozen inputs unchanged before vs after "
          f"({len(frozen_after)} files sha256-verified)", flush=True)

    new_sum = json.load(open(f"{new_clades}/tight_clades_summary.json"))["0"]
    old_sum = json.load(open(f"{old_clades}/tight_clades_summary.json"))["0"]
    med_of_meds = new_sum["internal_similarity"]["median_median"]
    comparison = {
        "version": v,
        "threshold": THRESHOLD,
        "max_size": MAX_SIZE,
        "n_prophages": n,
        "published_fragmented": {
            "n_clades": old_sum["n_clades"],
            "n_alignable": old_sum["n_clades"] - old_sum["n_singletons"],
            "n_singletons": old_sum["n_singletons"],
            "median_median": old_sum["internal_similarity"]["median_median"],
            "max_clade_size": old_sum["max_clade_size"],
            "clade_size_distribution": old_sum["clade_size_distribution"],
        },
        "corrected_sorted": {
            "n_clades": new_sum["n_clades"],
            "n_alignable": new_sum["n_clades"] - new_sum["n_singletons"],
            "n_singletons": new_sum["n_singletons"],
            "median_median": med_of_meds,
            "max_clade_size": new_sum["max_clade_size"],
            "clade_size_distribution": new_sum["clade_size_distribution"],
        },
        "delta_clades": new_sum["n_clades"] - old_sum["n_clades"],
        "checks": checks,
        "provenance": provenance,
        "frozen_input_sha256_before": frozen_before,
        "new_outdirs": {"mash": new_mash, "clades": new_clades},
        "frozen_outputs_touched": False,
    }
    with open(f"{new_clades}/sorted_vs_published.json", "w") as f:
        json.dump(comparison, f, indent=1)

    print(f"[{v}] RESULT: published {old_sum['n_clades']} clades "
          f"({old_sum['n_clades']-old_sum['n_singletons']} alignable + "
          f"{old_sum['n_singletons']} singletons, median "
          f"{old_sum['internal_similarity']['median_median']:.4f}) -> "
          f"corrected {new_sum['n_clades']} clades "
          f"({new_sum['n_clades']-new_sum['n_singletons']} alignable + "
          f"{new_sum['n_singletons']} singletons, median {med_of_meds:.4f}); "
          f"delta {new_sum['n_clades']-old_sum['n_clades']}", flush=True)
    print(f"[{v}] comparison -> {new_clades}/sorted_vs_published.json", flush=True)
    print(f"[{v}] DONE", flush=True)


if __name__ == "__main__":
    main()
