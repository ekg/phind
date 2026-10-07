#!/usr/bin/env python3
"""Tests for the structure-driven clustering mode of
scripts/build_tight_clades.py (`--split-mode medoids`).

Run:  python3 -m pytest scripts/test_medoid_split.py
      (or: python3 scripts/test_medoid_split.py)

Covers:
  * unit tests on small synthetic distance matrices: medoid / radius
    helpers, k-medoids (k=2) split shape, determinism, non-empty disjoint
    halves.
  * `split_until_radius`: a clade that already satisfies the radius is left
    alone (no cap -> one cluster); a spread clade is split until every
    member is within `--medoid-radius` of its medoid, with a split trace
    whose parent/child ids resolve against the emitted tight_clades.json.
  * end-to-end on a synthetic triangle: `--split-mode medoids` emits
    `split_trace.json` and ignores `--max-size`; the default invocation is
    byte-identical to `--split-mode cap` and emits no `split_trace.json`.
"""

import json
import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import build_tight_clades as btc  # noqa: E402

SCRIPT = os.path.join(HERE, "build_tight_clades.py")
THRESHOLD = 0.25
RADIUS = 0.10


# ─── helpers / fixtures ────────────────────────────────────────────────────

def line_distances(coords):
    """Pairwise MASH-like distances from 1-D coordinates: d(i, j) = |xi - xj|."""
    return [[abs(a - b) for b in coords] for a in coords]


def write_synthetic(root, coords, community=0):
    """Write ids.txt / labels.csv / triangle for a synthetic community whose
    members are the points in `coords` (triangle = |xi - xj| float32)."""
    n = len(coords)
    names = ["s%04d" % i for i in range(n)]
    dist = line_distances(coords)

    ids_path = os.path.join(root, "ids.txt")
    with open(ids_path, "w") as f:
        f.write("\n".join(names) + "\n")

    labels_path = os.path.join(root, "labels.csv")
    with open(labels_path, "w") as f:
        f.write("sequence,community\n")
        for s in names:
            f.write("%s,%d\n" % (s, community))

    tri_path = os.path.join(root, "synthetic.dist")
    tri = np.zeros(n * (n - 1) // 2, dtype="<f4")
    for i in range(n):
        for j in range(i + 1, n):
            tri[btc.triangle_offset(i, j, n)] = dist[i][j]
    tri.tofile(tri_path)
    return ids_path, labels_path, tri_path, names, dist


def run_build(root, tri_path, ids_path, labels_path, extra=(), outdir=None):
    outdir = outdir or os.path.join(root, "clades")
    cmd = [sys.executable, SCRIPT,
           "--triangle", tri_path,
           "--ids-file", ids_path,
           "--labels-csv", labels_path,
           "--communities", "0",
           "--threshold", str(THRESHOLD),
           "--outdir", outdir]
    cmd.extend(extra)
    subprocess.run(cmd, check=True, capture_output=True)
    return outdir


def read_json(path):
    with open(path) as f:
        return json.load(f)


def radius_of(clade_ids, names, dist):
    """Independent radius of a clade, straight from the raw distance matrix."""
    idx = [names.index(s) for s in clade_ids]
    if len(idx) <= 1:
        return 0.0
    sums = [sum(dist[i][j] for j in idx) for i in idx]
    med = idx[int(np.argmin(sums))]
    return max(dist[med][i] for i in idx)


# ─── unit tests on synthetic matrices ──────────────────────────────────────

class TestMedoidHelpers(unittest.TestCase):

    def test_medoid_and_radius_on_square(self):
        D = np.array([
            [0.0, 0.1, 0.1, 0.4],
            [0.1, 0.0, 0.2, 0.5],
            [0.1, 0.2, 0.0, 0.3],
            [0.4, 0.5, 0.3, 0.0],
        ], dtype=np.float32)
        # total distances: 0.6, 0.8, 0.6, 1.2 -> ties resolved to lowest index
        self.assertEqual(btc.medoid_of(D, [0, 1, 2, 3]), 0)
        self.assertAlmostEqual(btc.cluster_radius(D, [0, 1, 2, 3]), 0.4, places=6)
        self.assertEqual(btc.cluster_radius(D, [3]), 0.0)

    def test_kmedoids_split_is_disjoint_and_deterministic(self):
        coords = [-0.11, -0.10, 0.0, 0.10, 0.11]
        D = np.array(line_distances(coords), dtype=np.float32)
        idx = list(range(len(coords)))
        left, right, iters = btc.kmedoids_split(D, idx)
        self.assertTrue(left)
        self.assertTrue(right)
        self.assertEqual(sorted(left + right), idx)
        self.assertGreaterEqual(iters, 1)
        # deterministic on repeat
        self.assertEqual(btc.kmedoids_split(D, idx), (left, right, iters))
        # both halves are strictly smaller than the parent
        self.assertLess(len(left), len(idx))
        self.assertLess(len(right), len(idx))


class TestSplitUntilRadius(unittest.TestCase):

    def test_cluster_within_radius_is_untouched(self):
        coords = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05]
        D = np.array(line_distances(coords), dtype=np.float32)
        names = ["0_0000"]
        out, trace = btc.split_until_radius(D, [list(range(len(coords)))],
                                            names, RADIUS)
        self.assertEqual(len(out), 1)
        self.assertEqual(trace, [])
        self.assertEqual(sorted(out[0][1]), list(range(len(coords))))

    def test_spread_cluster_splits_until_within_radius(self):
        coords = [-0.11, -0.10, 0.0, 0.10, 0.11]   # all pairs <= threshold
        D = np.array(line_distances(coords), dtype=np.float32)
        root = list(range(len(coords)))
        out, trace = btc.split_until_radius(D, [root], ["0_0000"], RADIUS)
        self.assertGreater(len(out), 1)
        self.assertGreaterEqual(len(trace), 1)
        # every final clade satisfies the radius criterion
        for _, members in out:
            self.assertLessEqual(
                btc.cluster_radius(D, members), RADIUS + 1e-9,
                "clade %s exceeds radius %s" % (members, RADIUS))
        # members are partitioned
        flat = sorted(m for _, members in out for m in members)
        self.assertEqual(flat, root)
        # trace bookkeeping: sizes multiply out, radii shrink or hold
        for rec in trace:
            self.assertEqual(sum(rec["sizes_after"]), rec["size_before"])
            self.assertGreaterEqual(rec["iterations"], 1)
            for r in rec["radii_after"]:
                self.assertLessEqual(r, rec["radius_before"])

    def test_child_ids_resolve_against_final_clades(self):
        coords = [-0.11, -0.10, 0.0, 0.10, 0.11]
        D = np.array(line_distances(coords), dtype=np.float32)
        out, trace = btc.split_until_radius(
            D, [list(range(len(coords)))], ["0_0000"], RADIUS)
        final_ids = [n for n, _ in out]
        parents = {rec["parent_clade_id"] for rec in trace}
        children = {c for rec in trace for c in
                    (rec["parent_clade_id"] + "_s1", rec["parent_clade_id"] + "_s2")}
        # every parent is either an initial clade or a child of an earlier split
        for p in parents:
            self.assertTrue(p == "0_0000" or p in children,
                            "%s is neither initial nor a split child" % p)
        # every child is either emitted or split further, never a stale parent
        for cid in children:
            self.assertIn(cid, final_ids + sorted(parents),
                          "%s is neither a final clade nor split further" % cid)
        for cid in final_ids:
            self.assertNotIn(cid, parents, "%s was split but still emitted" % cid)

    def test_singleton_never_split(self):
        D = np.zeros((1, 1), dtype=np.float32)
        out, trace = btc.split_until_radius(D, [[0]], ["0_0000"], RADIUS)
        self.assertEqual(out, [("0_0000", [0])])
        self.assertEqual(trace, [])

    def test_missing_distance_is_split_out(self):
        """A NaN (missing) MASH distance cannot be certified within the radius,
        so the clade is split until no clade contains one; its radius is
        reported as null rather than leaking NaN into the trace."""
        coords = [-0.02, 0.0, 0.02, 0.04]
        D = np.array(line_distances(coords), dtype=np.float32)
        D[0, 1] = D[1, 0] = np.nan          # one missing distance inside
        self.assertTrue(np.isinf(btc.cluster_radius(D, [0, 1, 2, 3])))
        out, trace = btc.split_until_radius(D, [[0, 1, 2, 3]], ["0_0000"], RADIUS)
        self.assertGreaterEqual(len(trace), 1)
        self.assertIsNone(trace[0]["radius_before"])
        for _, members in out:
            self.assertFalse(np.isnan(D[np.ix_(members, members)]).any())
            self.assertLessEqual(btc.cluster_radius(D, members), RADIUS + 1e-9)
        self.assertEqual(sorted(m for _, ms in out for m in ms), [0, 1, 2, 3])


# ─── end-to-end on a synthetic triangle ────────────────────────────────────

class TestEndToEndSynthetic(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="medoid_split_")
        self.addCleanup(lambda: __import__("shutil").rmtree(self.tmp, True))

    def test_medoids_mode_no_cap_single_cluster(self):
        """All members within the radius -> one clade, even with --max-size 1."""
        coords = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
        ids, labels, tri, names, dist = write_synthetic(self.tmp, coords)
        out = run_build(self.tmp, tri, ids, labels,
                        extra=["--split-mode", "medoids",
                               "--medoid-radius", str(RADIUS),
                               "--max-size", "1"])
        tc = read_json(os.path.join(out, "0", "tight_clades.json"))
        self.assertEqual(len(tc), 1)
        self.assertEqual(sorted(next(iter(tc.values()))), sorted(names))
        trace = read_json(os.path.join(out, "0", "split_trace.json"))
        self.assertEqual(trace["n_splits"], 0)
        self.assertEqual(trace["medoid_radius"], RADIUS)

    def test_medoids_mode_splits_spread_cluster(self):
        coords = [-0.11, -0.10, 0.0, 0.10, 0.11]
        ids, labels, tri, names, dist = write_synthetic(self.tmp, coords)
        out = run_build(self.tmp, tri, ids, labels,
                        extra=["--split-mode", "medoids",
                               "--medoid-radius", str(RADIUS),
                               "--max-size", "1"])
        tc = read_json(os.path.join(out, "0", "tight_clades.json"))
        self.assertGreater(len(tc), 1)
        members = sorted(s for v in tc.values() for s in v)
        self.assertEqual(members, sorted(names))
        for cid, clade in tc.items():
            self.assertLessEqual(radius_of(clade, names, dist), RADIUS + 1e-6,
                                 "%s radius > %s" % (cid, RADIUS))
        trace = read_json(os.path.join(out, "0", "split_trace.json"))
        self.assertEqual(trace["n_splits"], len(trace["splits"]))
        self.assertGreaterEqual(trace["n_splits"], 1)
        for rec in trace["splits"]:
            self.assertEqual(set(rec.keys()),
                             {"parent_clade_id", "size_before", "sizes_after",
                              "radius_before", "radii_after", "iterations"})

    def test_default_mode_is_cap_mode_and_unchanged(self):
        """Default invocation == explicit --split-mode cap, byte for byte, and
        emits no split_trace.json."""
        coords = [-0.11, -0.10, 0.0, 0.10, 0.11]
        ids, labels, tri, names, dist = write_synthetic(self.tmp, coords)
        out_default = run_build(self.tmp, tri, ids, labels,
                                extra=["--max-size", "2"],
                                outdir=os.path.join(self.tmp, "default"))
        out_cap = run_build(self.tmp, tri, ids, labels,
                            extra=["--max-size", "2", "--split-mode", "cap"],
                            outdir=os.path.join(self.tmp, "cap"))
        for fname in ("tight_clades.json", "clade_similarity.json",
                      "members.json"):
            with open(os.path.join(out_default, "0", fname), "rb") as f:
                a = f.read()
            with open(os.path.join(out_cap, "0", fname), "rb") as f:
                b = f.read()
            self.assertEqual(a, b, "%s differs between default and cap" % fname)
        for out in (out_default, out_cap):
            self.assertFalse(os.path.exists(os.path.join(out, "0", "split_trace.json")))
            tc = read_json(os.path.join(out, "0", "tight_clades.json"))
            self.assertTrue(all(len(v) <= 2 for v in tc.values()))
            report = read_json(os.path.join(out, "0", "clade_similarity.json"))["report"]
            self.assertEqual(report["max_size"], 2)

    def test_cap_mode_under_cap_matches_v3_shape(self):
        """With no cap pressure the cap path still emits one clade id 0_0000."""
        coords = [0.0, 0.01, 0.02, 0.03]
        ids, labels, tri, names, dist = write_synthetic(self.tmp, coords)
        out = run_build(self.tmp, tri, ids, labels, extra=["--max-size", "800"])
        tc = read_json(os.path.join(out, "0", "tight_clades.json"))
        self.assertEqual(list(tc.keys()), ["0_0000"])
        self.assertEqual(sorted(tc["0_0000"]), sorted(names))
        self.assertFalse(os.path.exists(os.path.join(out, "0", "split_trace.json")))

    def test_medoids_mode_output_json_is_strict(self):
        """split_trace.json must be valid strict JSON (no NaN/Infinity tokens)."""
        coords = [-0.11, -0.10, 0.0, 0.10, 0.11]
        ids, labels, tri, names, dist = write_synthetic(self.tmp, coords)
        out = run_build(self.tmp, tri, ids, labels,
                        extra=["--split-mode", "medoids"])
        with open(os.path.join(out, "0", "split_trace.json")) as f:
            raw = f.read()
        self.assertNotIn("NaN", raw)
        self.assertNotIn("Infinity", raw)
        json.loads(raw, parse_constant=lambda c: (_ for _ in ()).throw(
            AssertionError("non-strict JSON constant: %s" % c)))


if __name__ == "__main__":
    unittest.main()
