#!/usr/bin/env python3
"""E. coli prophage bait & junction panel driver (preregistered: ../DESIGN.md).

Deterministic end-to-end design mirroring ntm/v2/bait/scripts/
design_ntm_baits.py and reusing the shared baitlib (implement-bait-design):

  pilot selection (diversity-aware, host-blind, mash-triangle medoids)
  -> traversal regeneration + byte-identity verification (junction provenance)
  -> bait generation (interior module / purely-reconstructed junction /
     member / public-reference positive / shuffled + host negatives)
  -> 31-mer masking with acceptance gates
  -> committed panel outputs + external ledgers.

Run:
  python3 design_ecoli_baits.py \
      --repo-root <repo> --external-root /mnt/nvme3n1/erikg/phind-genome-work/ecoli_bait \
      --out-dir <repo>/research/ecoli_bait/panel

Committed outputs are byte-deterministic (no wall-clock inside them).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "..", "ntm", "v2",
                                "bait", "scripts"))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "..", "scripts"))
import baitlib as bl  # noqa: E402

# ---------------------------------------------------------------------------
# Preregistered parameters (DESIGN.md)
# ---------------------------------------------------------------------------

SELECT_N_A = 18
SELECT_N_B = 6
MIN_GENOME_BP = 5000
MEMBER_MIN_BP = 600
MODULE_W = 1200
JUNCTION_W = 1000
MEMBER_W = 1200
HOSTCTL_W = 1000
POSCON_W = 1200
MODULE_CAP = 1                 # one interior bait per clade (DESIGN.md)
JUNCTION_CAP = 3
N_SHUFFLED_CONTROLS = 2
N_HOST_CONTROLS = 2
SHUFFLED_SEED = 42
SHARED_KMER_CONTROL_MAX = 5
TIERS = ("A", "B")
TIER_B_ALLOWED_FLAGS = {"host_genes", "low_completeness"}
RELEASE_FLAG_EXCLUDE = ("too_short", "long_n_run")
HOST_PANEL_MAX = 60
K12_ACCESSION = "GCF_000005845.2"   # E. coli K-12 MG1655
PUBLIC_REFS = {                    # accession -> (label, expected length)
    "NC_001416.1": ("LAMBDA", 48502),
    "NC_002167.1": ("HK97", 39732),
}
TRAV_N_SAMPLES = 5
TRAV_SEED = 42
TRAV_MIN_BUDGET = 30000
TRAV_MEMBER_MULT = 1.6
INVALID = np.uint64(1) << np.uint64(62)

BAIT_CLASSES = (
    "interior_module", "interior_generic", "junction", "member_interior",
    "control_shuffled_negative", "control_host_negative",
    "control_public_reference",
)

COHORT = "ecoli_ml_reconstructed"


def log(msg: str) -> None:
    print(msg, flush=True)


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: str) -> List[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: str, rows: List[dict], cols: List[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            line = "\t".join(str(r.get(c, "")) for c in cols).rstrip("\t")
            fh.write(line + "\n")  # no trailing whitespace anywhere


def iter_fasta(path: str):
    """Stream (id, sequence) from a plain or gzipped FASTA."""
    op = gzip_open if path.endswith((".gz", ".bgz")) else open
    with op(path, "rt") as fh:  # type: ignore[operator]
        cur, buf = None, []
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    yield cur, "".join(buf)
                cur, buf = line[1:].split()[0], []
            elif cur is not None:
                buf.append(line.strip())
        if cur is not None:
            yield cur, "".join(buf)


def gzip_open(path: str, mode: str = "rt"):
    import gzip
    return gzip.open(path, mode)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

class Inputs:
    def __init__(self, repo: str, clades_source: str, prophages_fa: str,
                 external_root: str, host_genomes_dir: str, pharokka_cds: str):
        self.repo = repo
        self.clades_source = clades_source          # main-tree research/clades
        self.prophages_fa = prophages_fa            # main-tree full_prophages.fa
        self.external_root = external_root
        self.ml_dir = os.path.join(repo, "research", "ml_phage_genomes")
        self.ml_fasta = os.path.join(self.ml_dir, "all_ml_phage_genomes.fa.gz")
        self.release_manifest = os.path.join(self.ml_dir, "release_manifest.tsv")
        self.release_qc = os.path.join(self.ml_dir, "qc_table.tsv")
        self.functional_qc = os.path.join(repo, "research", "phage_annotation",
                                          "per_genome_annotation_qc.tsv")
        self.clade_genomes = os.path.join(self.ml_dir, "clade_genomes.json")
        self.mash_dir = os.path.join(repo, "research", "mash_tree")
        self.triangle = os.path.join(self.mash_dir, "full_prophages_mash.dist")
        self.triangle_ids = os.path.join(self.mash_dir, "data", "ids.txt")
        self.prophage_calls = os.path.join(repo, "26k_prophage1.csv")
        self.pharokka_cds = pharokka_cds
        self.host_genomes_dir = host_genomes_dir
        self.traverse_script = os.path.join(repo, "scripts",
                                            "traverse_partitions.py")
        self.public_refs = {acc: os.path.join(repo, "research", "ecoli_bait",
                                              "panel", "public_refs",
                                              f"{acc}.fa")
                            for acc in PUBLIC_REFS}
        self.work = os.path.join(external_root, "work")
        self.staged_clades = os.path.join(self.work, "clades")
        self.trav_dir = os.path.join(self.work, "traversal")
        self.member_fa = os.path.join(self.work, "member_prophages.fa")

    def community_tight_clades(self) -> Dict[str, Dict[str, List[str]]]:
        out: Dict[str, List[str]] = {}
        root = os.path.join(self.repo, "research", "clades")
        for comm in sorted(os.listdir(root)):
            p = os.path.join(root, comm, "tight_clades.json")
            if os.path.isfile(p):
                with open(p) as fh:
                    out[comm] = json.load(fh)
        return out

    def clade_dir_source(self, community: str, clade: str) -> str:
        return os.path.join(self.clades_source, community, clade)

    def clade_dir_staged(self, community: str, clade: str) -> str:
        return os.path.join(self.staged_clades, f"{community}_{clade}")

    def host_fasta(self, accession: str) -> str:
        return os.path.join(self.host_genomes_dir, accession,
                            f"{accession}.pansn.fa.gz")


# ---------------------------------------------------------------------------
# Functional-QC tier derivation (preregistered, DESIGN.md stage 1)
# ---------------------------------------------------------------------------

def derive_tier(flag: str) -> Optional[str]:
    flags = set(x for x in (flag or "").split(";") if x)
    if not flags or flags == {"ok"}:
        return "A"
    if flags <= TIER_B_ALLOWED_FLAGS:
        return "B"
    return None


# ---------------------------------------------------------------------------
# Stage 1 — pilot selection (deterministic, diversity-aware, host-blind)
# ---------------------------------------------------------------------------

class Triangle:
    """Condensed float32 mash triangle with id index (research/mash_tree)."""

    def __init__(self, dist_path: str, ids_path: str):
        self.ids = [l.strip() for l in open(ids_path) if l.strip()]
        self.index = {p: i for i, p in enumerate(self.ids)}
        self.n = len(self.ids)
        expected = self.n * (self.n - 1) // 2
        actual = os.path.getsize(dist_path) // 4
        if actual != expected:
            raise ValueError(f"triangle size mismatch: {actual} vs {expected}")
        self.tri = np.memmap(dist_path, dtype="<f4", mode="r",
                             shape=(expected,))

    def cond_index(self, i: int, j: int) -> int:
        if i > j:
            i, j = j, i
        return self.n * i - i * (i + 1) // 2 + (j - i - 1)

    def clade_medoid_fast(self, member_ids: Sequence[str]) -> str:
        """Medoid minimising mean mash distance over co-members.

        Uses contiguous row-slice gathers from the condensed triangle (row i
        is contiguous in j), NaN -> 1.0, lexicographic-id ties.
        """
        ids = sorted(member_ids)
        arr = np.asarray(sorted(self.index[m] for m in ids), dtype=np.int64)
        k = arr.size
        D = np.zeros((k, k), dtype=np.float64)
        for r in range(k - 1):
            i = int(arr[r])
            j_rest = arr[r + 1:]
            base = self.cond_index(i, int(j_rest[0]))
            end = self.cond_index(i, int(j_rest[-1])) + 1
            offs = j_rest - i - 1
            vals = np.asarray(self.tri[base:end], dtype=np.float64)[offs]
            vals[np.isnan(vals)] = 1.0
            D[r, r + 1:] = vals
            D[r + 1:, r] = vals
        scores = D.mean(axis=1)
        best = float(scores.min())
        tied = [ids[r] for r in range(k) if scores[r] <= best + 1e-12]
        return min(tied)

    def dist(self, a: str, b: str) -> float:
        if a == b:
            return 0.0
        if a not in self.index or b not in self.index:
            return 1.0
        v = float(self.tri[self.cond_index(self.index[a], self.index[b])])
        return 1.0 if np.isnan(v) else v


def select_pilot(inp: Inputs, tri: Triangle
                 ) -> Tuple[List[dict], List[dict], List[dict], Dict[str, str]]:
    """Returns (selected genomes, member reps, eligibility ledger, medoids)."""
    qc = {r["clade_id"]: r for r in read_tsv(inp.functional_qc)
          if r["source"] == "ecoli_ml"}
    rel = {r["clade_id"]: r for r in read_tsv(inp.release_manifest)}
    relq = {r["clade_id"]: r for r in read_tsv(inp.release_qc)}
    tight_all = inp.community_tight_clades()

    ledger: List[dict] = []
    pools: Dict[str, List[str]] = {t: [] for t in TIERS}
    tier_by_clade: Dict[str, str] = {}
    members_by_clade: Dict[str, List[str]] = {}

    for clade in sorted(rel):
        genome = f"clade_{clade}_ML"
        q = qc.get(genome)
        r, rq = rel[clade], relq.get(genome, {})
        tier = derive_tier(q["flag"]) if q else None
        reason = None
        if q is None:
            reason = "no_functional_qc_row"
        else:
            if tier is None:
                reason = f"tier_ineligible:flags={q['flag']}"
            elif r["status"] != "ml":
                reason = f"release_status={r['status']}"
            elif int(r["length_bp"]) < MIN_GENOME_BP:
                reason = f"release_length<{MIN_GENOME_BP}"
            elif int(r["n_members"]) < 2:
                reason = "n_members<2"
            elif any(f in (rq.get("flags") or "")
                     for f in RELEASE_FLAG_EXCLUDE):
                reason = f"release_flags={rq.get('flags')}"
            else:
                comm = clade.split("_")[0]
                mem = sorted(m for m in tight_all.get(comm, {}).get(clade, [])
                             if m in tri.index)
                if len(mem) < 2:
                    reason = "fewer_than_2_triangle_members"
        ledger.append({
            "genome_id": genome, "clade_id": clade,
            "functional_qc_flag": q["flag"] if q else "",
            "tier": tier or "",
            "eligible": "yes" if reason is None else "no",
            "reason": reason or "eligible",
        })
        if reason is None:
            pools[tier].append(clade)             # type: ignore[index]
            tier_by_clade[clade] = tier            # type: ignore[index]
            members_by_clade[clade] = mem           # type: ignore[index]

    # medoids over prophage mash space (host labels never enter here)
    medoids: Dict[str, str] = {}
    for clade in sorted(members_by_clade):
        medoids[clade] = tri.clade_medoid_fast(members_by_clade[clade])

    cand_all = sorted(set(pools["A"]) | set(pools["B"]))
    dist: Dict[Tuple[str, str], float] = {}
    for a in cand_all:
        for b in cand_all:
            if a != b:
                dist[(a, b)] = tri.dist(medoids[a], medoids[b])

    a_pool = sorted(pools["A"])
    b_pool = sorted(pools["B"])
    selected: List[str] = []
    if a_pool:
        selected += bl.greedy_farthest_point(a_pool, dist, a_pool[0],
                                             SELECT_N_A)
    quota_b = SELECT_N_A + SELECT_N_B - len(selected)
    if quota_b > 0 and b_pool:
        selected += bl.greedy_farthest_point(b_pool, dist, b_pool[0],
                                             quota_b)
    if len(selected) < SELECT_N_A + SELECT_N_B:   # pool shortfall -> fill
        rest = [c for c in a_pool + b_pool if c not in selected]
        selected += rest[:SELECT_N_A + SELECT_N_B - len(selected)]

    member_reps: List[dict] = []
    for rank, clade in enumerate(selected, 1):
        mem = members_by_clade[clade]
        medoid = medoids[clade]
        member_reps.append({
            "selection_rank": rank, "clade_id": clade,
            "member_prophage_id": medoid, "role": "medoid",
            "is_medoid": "yes",
        })
        if len(mem) >= 2:
            order = sorted(((tri.dist(medoid, m), m) for m in mem
                            if m != medoid), key=lambda t: (-t[0], t[1]))
            if order:
                member_reps.append({
                    "selection_rank": rank, "clade_id": clade,
                    "member_prophage_id": order[0][1],
                    "role": "farthest_from_medoid", "is_medoid": "no",
                })

    selected_rows = []
    for rank, clade in enumerate(selected, 1):
        r = rel[clade]
        selected_rows.append({
            "selection_rank": rank, "clade_id": clade,
            "genome_id": f"clade_{clade}_ML", "tier": tier_by_clade[clade],
            "n_members": r["n_members"], "n_partitions": r["n_partitions"],
            "length_bp": r["length_bp"], "release_flags": "",
            "functional_qc_flag": qc[genome_of(clade)]["flag"],
            "medoid_member": medoids[clade],
            "selection_basis": "greedy_farthest_point_prophage_mash_"
                               "space_host_blind",
        })
    write_tsv(os.path.join(inp.work, "eligibility_ledger.tsv"), ledger,
              ["genome_id", "clade_id", "functional_qc_flag", "tier",
               "eligible", "reason"])
    return selected_rows, member_reps, ledger, medoids


def genome_of(clade: str) -> str:
    return f"clade_{clade}_ML"


def clade_of(genome: str) -> str:
    # clade_0_0000_ML -> 0_0000 (community may be multi-digit)
    return genome[len("clade_"):-len("_ML")]


# ---------------------------------------------------------------------------
# Staging: clade partitions + deterministic traversal regeneration
# ---------------------------------------------------------------------------

def stage_clade(inp: Inputs, community: str, clade: str) -> Optional[str]:
    """Copy partitions.bed + partitions/ + manifest.json to NVMe (once)."""
    src = inp.clade_dir_source(community, clade)
    bed = os.path.join(src, "partitions.bed")
    parts = os.path.join(src, "partitions")
    if not (os.path.isfile(bed) and os.path.isdir(parts)):
        return None
    dst = inp.clade_dir_staged(community, clade)
    if not os.path.exists(os.path.join(dst, "partitions.bed")):
        os.makedirs(dst, exist_ok=True)
        shutil.copy2(bed, os.path.join(dst, "partitions.bed"))
        if not os.path.isdir(os.path.join(dst, "partitions")):
            shutil.copytree(parts, os.path.join(dst, "partitions"))
        mf = os.path.join(src, "manifest.json")
        if os.path.isfile(mf):
            shutil.copy2(mf, os.path.join(dst, "manifest.json"))
    return dst


def traversal_budget(staged_dir: str) -> int:
    m = json.load(open(os.path.join(staged_dir, "manifest.json")))
    mx = int(m["sequence_lengths"]["max"])
    return max(TRAV_MIN_BUDGET, int(TRAV_MEMBER_MULT * mx))


def regenerate_traversal(inp: Inputs, community: str, clade: str,
                         staged_dir: str) -> str:
    """Run the seeded traversal; returns the traversal.json path.

    Deterministic (same inputs + seed + budget). Cached per clade.
    """
    os.makedirs(inp.trav_dir, exist_ok=True)
    prefix = os.path.join(inp.trav_dir, f"{community}_{clade}")
    tjson = prefix + ".traversal.json"
    if os.path.exists(tjson):
        return tjson
    budget = traversal_budget(staged_dir)
    cmd = [sys.executable, inp.traverse_script,
           "--partitions-dir", os.path.join(staged_dir, "partitions"),
           "--bed", os.path.join(staged_dir, "partitions.bed"),
           "--output", prefix, "--mode", "ml",
           "--n-samples", str(TRAV_N_SAMPLES), "--seed", str(TRAV_SEED),
           "--max-length", str(budget)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"traversal failed for {clade}: {r.stderr[-1500:]}")
    return tjson


def verify_traversal(trav_ml_fa: str, released_seq: str) -> Tuple[bool, int, int]:
    seq = None
    for _hid, s in iter_fasta(trav_ml_fa):
        seq = s.upper()
        break
    if seq is None:
        return False, 0, len(released_seq)
    return seq == released_seq, len(seq), len(released_seq)


def load_released_ml(inp: Inputs, wanted: set) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    for hid, seq in iter_fasta(inp.ml_fasta):
        if hid in wanted:
            seqs[hid] = seq.upper()
    return seqs


# ---------------------------------------------------------------------------
# Junction candidates (purely-reconstructed joins only)
# ---------------------------------------------------------------------------

def junction_candidates(bed_path: str, trav_json: str) -> List[dict]:
    """Ranked purely-reconstructed joins: both partitions + exact breakpoint.

    E. coli rule (DESIGN.md): only joins with observed_adjacency_count == 0
    are eligible; ranking = co_occurrence desc, min(occurrence) desc, pids.
    """
    trav = json.load(open(trav_json))
    parts = trav["partitions"]
    pids = trav["genome"]["pids"]
    spans: Dict[int, Tuple[int, int]] = {}
    pos = 0
    for pid in pids:
        L = int(parts[str(pid)]["representative_len"])
        spans[pid] = (pos + 1, pos + L)
        pos += L

    per_member: Dict[str, List[Tuple[int, int, int]]] = {}
    with open(bed_path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 4:
                try:
                    pid = int(f[3])
                except ValueError:
                    continue
                if pid in spans:
                    per_member.setdefault(f[0], []).append(
                        (int(f[1]), int(f[2]), pid))

    cands: List[dict] = []
    for a, b in zip(pids, pids[1:]):
        pa, pb = parts[str(a)], parts[str(b)]
        if not (pa.get("has_seq", True) and pb.get("has_seq", True)):
            continue
        co_occ = 0
        obs_adj = 0
        for ivs in per_member.values():
            pids_here = {p for _, _, p in ivs}
            if a in pids_here and b in pids_here:
                co_occ += 1
                ivs.sort()
                for (s1, e1, p1), (s2, e2, p2) in zip(ivs, ivs[1:]):
                    if {p1, p2} == {a, b} and e1 == s2:
                        obs_adj += 1
                        break
        if obs_adj != 0:
            continue  # purely-reconstructed joins only
        cands.append({
            "partition_a": a, "partition_b": b, "breakpoint": spans[b][0],
            "occ_a": int(pa["occurrence"]), "occ_b": int(pb["occurrence"]),
            "co_occurrence_count": co_occ, "observed_adjacency_count": 0,
            "pa_coords": f"{spans[a][0]}-{spans[a][1]}",
            "pb_coords": f"{spans[b][0]}-{spans[b][1]}",
            "support": min(int(pa["occurrence"]), int(pb["occurrence"])),
        })
    cands.sort(key=lambda c: (-c["co_occurrence_count"], -c["support"],
                              c["partition_a"], c["partition_b"]))
    return cands


# ---------------------------------------------------------------------------
# Stage 2 — bait candidates
# ---------------------------------------------------------------------------

class Bait:
    __slots__ = ("bait_id", "bait_class", "clade_id", "genome_id", "member_id",
                 "contig", "start", "end", "seq", "module", "phrog", "annot",
                 "gene_lo", "gene_hi", "partition_a", "partition_b",
                 "breakpoint", "occ_a", "occ_b", "co_occ", "obs_adj",
                 "pa_coords", "pb_coords", "dominated", "priority",
                 "expected_behavior", "notes", "stats")

    def length(self) -> int:
        return self.end - self.start + 1


def new_bait(**kw) -> Bait:
    b = Bait()
    for s in Bait.__slots__:
        setattr(b, s, kw.get(s))
    b.stats = {}
    return b


def load_cds_rows(path: str, genomes: set) -> Dict[str, List[dict]]:
    """genome_id (clade_X_ML) -> pharokka rows (contig ecoli_ml|clade_X_ML)."""
    prefix = {f"ecoli_ml|{g}": g for g in genomes}
    rows: Dict[str, List[dict]] = {g: [] for g in genomes}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = prefix.get(row.get("contig") or "")
            if g is not None:
                rows[g].append(row)
    return rows


def load_prophage_intervals(path: str) -> Tuple[
        Dict[str, Dict[str, List[Tuple[int, int]]]], Dict[str, int]]:
    """26k_prophage1.csv -> (accession -> contig -> [(s,e)], prophage len).

    Also returns prophage_id -> length for member-bait provenance.
    """
    iv: Dict[str, Dict[str, List[Tuple[int, int]]]] = {}
    lens: Dict[str, int] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            pid = row.get("prophage_id") or ""
            if not pid:
                continue
            try:
                s = int(round(float(row["begin"])))
                e = int(round(float(row["end"])))
            except (ValueError, TypeError, KeyError):
                continue
            if s > e:
                s, e = e, s
            lens[pid] = e - s + 1
            acc = row["genome"]
            iv.setdefault(acc, {}).setdefault(
                row["scaffold"], []).append((s, e))
    for acc in iv:
        for contig in iv[acc]:
            iv[acc][contig].sort()
    return iv, lens


def generate_baits(inp: Inputs, selected: List[dict], member_reps: List[dict],
                   ml_seqs: Dict[str, str], cds_rows: Dict[str, List[dict]],
                   proph_seqs: Dict[str, str],
                   staged: Dict[str, str], trav_verified: Dict[str, bool]
                   ) -> Tuple[List[Bait], List[dict]]:
    baits: List[Bait] = []
    counters: Dict[Tuple[str, str], int] = {}
    rej: List[dict] = []

    def bait_id(clade: str, cls: str) -> str:
        key = (clade, cls)
        counters[key] = counters.get(key, 0) + 1
        return f"ECBAIT_{clade}_{cls.upper()}_{counters[key]:02d}"

    for row in selected:
        clade = row["clade_id"]
        genome = row["genome_id"]
        comm = clade.split("_")[0]
        seq = ml_seqs[genome]
        glen = len(seq)
        rows = cds_rows.get(genome, [])
        used: List[Tuple[int, int]] = []
        n_module = 0
        for gene in bl.pick_module_genes(rows):
            if n_module >= MODULE_CAP:
                break
            lo, hi = sorted((int(gene["start"]), int(gene["stop"])))
            s, e = bl.centered_window(glen, (lo + hi) // 2, MODULE_W)
            if any(bl.windows_overlap((s, e), u) for u in used):
                continue
            baits.append(new_bait(
                bait_id=bait_id(clade, "interior_module"),
                bait_class="interior_module", clade_id=clade, genome_id=genome,
                member_id="", contig=genome, start=s, end=e,
                seq=seq[s - 1:e], module=gene["module"],
                phrog=gene.get("phrog", ""), annot=gene.get("annot", ""),
                gene_lo=lo, gene_hi=hi,
                dominated=bl.cds_composition_flags(rows, s, e)[0],
                priority=(0, row["selection_rank"], n_module)))
            used.append((s, e))
            n_module += 1
        if n_module < MODULE_CAP:  # generic midpoint fill to 1/clade
            s, e = bl.centered_window(glen, glen // 2 + 1, MODULE_W)
            if not any(bl.windows_overlap((s, e), u) for u in used):
                baits.append(new_bait(
                    bait_id=bait_id(clade, "interior_generic"),
                    bait_class="interior_generic", clade_id=clade,
                    genome_id=genome, member_id="", contig=genome, start=s,
                    end=e, seq=seq[s - 1:e], module="generic_interior",
                    dominated=bl.cds_composition_flags(rows, s, e)[0],
                    priority=(1, row["selection_rank"], 0)))
                used.append((s, e))
        # junctions: regenerated + verified traversal only (DESIGN.md)
        nj = 0
        if staged.get(clade) and trav_verified.get(clade):
            try:
                jcands = junction_candidates(
                    os.path.join(staged[clade], "partitions.bed"),
                    os.path.join(inp.trav_dir, f"{comm}_{clade}.traversal.json"))
            except (OSError, KeyError, ValueError, json.JSONDecodeError):
                jcands = []
            for jc in jcands:
                if nj >= JUNCTION_CAP:
                    break
                bp = jc["breakpoint"]
                s, e = bl.centered_window(glen, bp, JUNCTION_W)
                if not (s <= bp <= e + 1):  # must span the join
                    continue
                baits.append(new_bait(
                    bait_id=bait_id(clade, "junction"),
                    bait_class="junction", clade_id=clade, genome_id=genome,
                    member_id="", contig=genome, start=s, end=e,
                    seq=seq[s - 1:e], module="partition_junction",
                    partition_a=jc["partition_a"],
                    partition_b=jc["partition_b"], breakpoint=bp,
                    occ_a=jc["occ_a"], occ_b=jc["occ_b"],
                    co_occ=jc["co_occurrence_count"],
                    obs_adj=jc["observed_adjacency_count"],
                    pa_coords=jc["pa_coords"], pb_coords=jc["pb_coords"],
                    dominated=bl.cds_composition_flags(rows, s, e)[0],
                    priority=(2, row["selection_rank"], nj)))
                nj += 1
        elif staged.get(clade) and not trav_verified.get(clade, False):
            rej.append({"clade_id": clade, "reason":
                        "traversal_regeneration_mismatch"})
        # member bait: medoid representative
        med = next((m for m in member_reps
                    if m["clade_id"] == clade and m["is_medoid"] == "yes"),
                   None)
        if med:
            pid = med["member_prophage_id"]
            mseq = proph_seqs.get(pid)
            if mseq:
                w = min(MEMBER_W, len(mseq))
                if w >= bl.BAIT_MIN_BP:
                    s, e = bl.centered_window(len(mseq), len(mseq) // 2 + 1, w)
                    baits.append(new_bait(
                        bait_id=bait_id(clade, "member_interior"),
                        bait_class="member_interior", clade_id=clade,
                        genome_id=genome, member_id=pid, contig=pid,
                        start=s, end=e, seq=mseq[s - 1:e],
                        module="member_midpoint", dominated=False,
                        priority=(3, row["selection_rank"], 0)))
    baits.sort(key=lambda b: b.priority)
    return baits, rej


# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------

def make_public_reference_baits(inp: Inputs) -> List[Bait]:
    """Mid-genome 1200 bp slices of lambda + HK97 (exact accessions)."""
    baits: List[Bait] = []
    for i, (acc, (label, exp_len)) in enumerate(sorted(PUBLIC_REFS.items()),
                                                1):
        path = inp.public_refs[acc]
        seq = None
        for _hid, s in iter_fasta(path):
            seq = s.upper()
            break
        if seq is None:
            raise RuntimeError(f"public reference missing: {path}")
        if len(seq) != exp_len:
            raise RuntimeError(f"{acc}: length {len(seq)} != expected "
                               f"{exp_len}")
        s, e = bl.centered_window(len(seq), len(seq) // 2 + 1, POSCON_W)
        baits.append(new_bait(
            bait_id=f"ECBAIT_POSCON_{label}_{i:02d}",
            bait_class="control_public_reference", clade_id="",
            genome_id="", member_id="", contig=acc, start=s, end=e,
            seq=seq[s - 1:e], module="public_reference_poscon",
            dominated=False, expected_behavior="own_reference_hit",
            notes=(f"positive control: mid-genome slice of RefSeq {acc} "
                   f"({label}); ref_len={len(seq)}; "
                   f"ref_sha256={sha256_file(path)}"),
            priority=(4, i, 0)))
    return baits


def make_host_control_baits(inp: Inputs, selected: List[dict],
                            member_reps: List[dict],
                            intervals: Dict[str, Dict[str, List[Tuple[int, int]]]],
                            prophage_calls: Dict[str, dict]
                            ) -> List[Bait]:
    """1-kb windows from prophage-free stretches: K-12 first, then a source
    host of the first selected clade whose rep source differs."""
    order: List[Tuple[str, str]] = [(K12_ACCESSION, "k12_reference")]
    for row in selected:
        reps = [m for m in member_reps
                if m["clade_id"] == row["clade_id"] and m["is_medoid"] == "yes"]
        for m in reps:
            acc = prophage_calls.get(m["member_prophage_id"], {}).get("genome")
            if acc and acc != K12_ACCESSION and all(a != acc for a, _ in order):
                order.append((acc, "selected_clade_source_host"))
        if len(order) >= N_HOST_CONTROLS:
            break
    baits: List[Bait] = []
    for acc, role in order[:N_HOST_CONTROLS]:
        if not os.path.exists(inp.host_fasta(acc)):
            continue
        best = None  # (stretch_len, contig, start, end, seq)
        for contig, seq in bl.read_bgzf_records(inp.host_fasta(acc)):
            iv = sorted(intervals.get(acc, {}).get(contig.split("#")[-1], []))
            clean = np.ones(len(seq), dtype=bool)
            for s, e in iv:
                s, e = max(1, s), min(len(seq), e)
                if s <= e:
                    clean[s - 1:e] = False
            if not clean.any():
                continue
            edges = np.diff(np.concatenate(
                ([False], clean, [False])).astype(np.int8))
            change = np.flatnonzero(edges)
            starts, ends = change[0::2], change[1::2]
            lens = ends - starts
            j = int(np.argmax(lens))
            if lens[j] < HOSTCTL_W or (best is not None and
                                       lens[j] <= best[0]):
                continue
            mid = int(starts[j] + lens[j] // 2) + 1
            s, e = bl.centered_window(len(seq), mid, HOSTCTL_W)
            best = (int(lens[j]), contig, s, e, seq.upper())
        if best is None:
            continue
        _, contig, s, e, seq = best
        baits.append(new_bait(
            bait_id=f"ECBAIT_HOST_CONTROL_{len(baits) + 1:02d}",
            bait_class="control_host_negative", clade_id="", genome_id="",
            member_id="", contig=contig, start=s, end=e,
            seq=seq[s - 1:e], module="host_control", dominated=False,
            expected_behavior="host_like",
            notes=(f"host-derived negative/filter-validation control "
                   f"({role} {acc}): prophage-free stretch of a source-host "
                   f"contig"),
            priority=(5, len(baits) + 1, 0)))
    return baits


def make_shuffled_controls(accepted: List[Bait],
                           host_like_union: np.ndarray) -> List[Bait]:
    """Fixed-seed Markov negatives from the first accepted module baits."""
    baits: List[Bait] = []
    panel_union = np.unique(np.concatenate(
        [b.stats["_canon"][b.stats["_valid_pos"]] for b in accepted])) \
        if accepted else np.empty(0, dtype=np.uint64)
    seen_genomes: set = set()
    for b in accepted:
        if (b.bait_class != "interior_module" or b.genome_id in seen_genomes
                or len(baits) >= N_SHUFFLED_CONTROLS):
            continue
        seed = SHUFFLED_SEED + len(baits)
        tries = 0
        while True:
            cand = bl.markov_control_sequence(b.seq, seed=seed)
            canon = np.unique(bl.canonical_codes(cand))
            shared = int(np.intersect1d(canon, panel_union).size) \
                if panel_union.size else 0
            if shared <= SHARED_KMER_CONTROL_MAX or tries >= 10:
                break
            tries += 1
            seed += 1000
        dev = bl.composition_deviation(cand, b.seq)
        hl = int(np.intersect1d(canon, host_like_union).size) \
            if host_like_union.size else 0
        seen_genomes.add(b.genome_id)
        nb = new_bait(
            bait_id=f"ECBAIT_SHUF_CONTROL_{len(baits) + 1:02d}",
            bait_class="control_shuffled_negative", clade_id=b.clade_id,
            genome_id="", member_id="", contig=b.contig, start=b.start,
            end=b.end, seq=cand, module="markov_negative", dominated=False,
            expected_behavior="no_hits",
            notes=(f"source_bait={b.bait_id};seed={seed};"
                   f"shared_panel_kmers={shared};host_like_kmers={hl};"
                   f"mono_l1={dev['mono_l1']:.4f};"
                   f"dinuc_l1={dev['dinuc_l1']:.4f};redraws={tries}"),
            priority=(6, len(baits) + 1, 0))
        compute_masks(nb)
        nb.stats["n_host_like"] = hl
        nb.stats["_hostmask"].fill(False)
        if hl:
            vals = nb.stats["_canon"][nb.stats["_valid_pos"]]
            idx = np.clip(np.searchsorted(host_like_union, vals), 0,
                          host_like_union.size - 1)
            nb.stats["_hostmask"][nb.stats["_valid_pos"][
                host_like_union[idx] == vals]] = True
        usable = ~(nb.stats["_ambig"] | nb.stats["_lcx"] |
                   nb.stats["_dup"] | nb.stats["_hostmask"])
        nb.stats["n_usable_kmers"] = int(usable.sum())
        nb.stats["cross_bait_max_shared_frac"] = round(
            shared / max(1, int(nb.stats["_valid_pos"].size)), 4)
        baits.append(nb)
    return baits


# ---------------------------------------------------------------------------
# Stage 3 — masking + gates (shared with NTM driver)
# ---------------------------------------------------------------------------

def compute_masks(b: Bait) -> None:
    arr = bl.seq_to_arr(b.seq)
    canon = bl.canonical_codes(arr)
    valid_pos = np.flatnonzero(canon != INVALID)
    dup = bl.internal_dup_mask(canon[valid_pos]) if valid_pos.size else \
        np.zeros(0, dtype=bool)
    dup_full = np.zeros(canon.size, dtype=bool)
    dup_full[valid_pos] = dup
    valid_bases = int((arr != 255).sum())
    gc = float(((arr == 1) | (arr == 2)).sum()) / max(1, valid_bases)
    b.stats.update({
        "n_kmers_total": int(canon.size),
        "n_ambiguous": int(bl.ambiguous_mask(arr).sum()),
        "n_low_complexity": int(bl.low_complexity_mask(arr).sum()),
        "n_internal_dup": int(dup_full.sum()),
        "gc_frac": round(gc, 4),
        "_canon": canon, "_valid_pos": valid_pos, "_dup": dup_full,
        "_ambig": bl.ambiguous_mask(arr), "_lcx": bl.low_complexity_mask(arr),
        "_hostmask": np.zeros(canon.size, dtype=bool),
    })


def host_scan(inp: Inputs, baits: List[Bait], panel: List[dict],
              intervals: Dict[str, Dict[str, List[Tuple[int, int]]]],
              log_path: str) -> np.ndarray:
    """Single host pass marking host-like canonical 31-mers per bait."""
    if baits:
        union = np.unique(np.concatenate(
            [b.stats["_canon"][b.stats["_valid_pos"]] for b in baits]))
    else:
        union = np.empty(0, dtype=np.uint64)
    host_like_union = np.empty(0, dtype=np.uint64)
    with open(log_path, "w") as lf:
        lf.write("accession\tcontigs\tprophage_intervals\tunion_kmers_hit\n")
        for entry in panel:
            acc = str(entry["accession"])
            iv = intervals.get(acc, {})
            ncontigs = 0
            found_all = np.empty(0, dtype=np.uint64)
            for contig, seq in bl.read_bgzf_records(inp.host_fasta(acc)):
                ncontigs += 1
                for (s, e) in iv.get(contig.split("#")[-1], []):
                    s, e = max(1, s), min(len(seq), e)
                    if s <= e:
                        seq = seq[:s - 1] + "N" * (e - s + 1) + seq[e:]
                canon = bl.canonical_codes(seq)
                if canon.size and union.size:
                    found_all = np.union1d(
                        found_all, np.intersect1d(canon, union))
            host_like_union = np.union1d(host_like_union, found_all)
            lf.write(f"{acc}\t{ncontigs}\t"
                     f"{sum(len(v) for v in iv.values())}\t"
                     f"{found_all.size}\n")
    for b in baits:
        m = b.stats["_hostmask"]
        vals = b.stats["_canon"][b.stats["_valid_pos"]]
        if vals.size and host_like_union.size:
            idx = np.clip(np.searchsorted(host_like_union, vals), 0,
                          host_like_union.size - 1)
            hit = host_like_union[idx] == vals
            m[b.stats["_valid_pos"][hit]] = True
            b.stats["n_host_like"] = int(hit.sum())
        else:
            b.stats["n_host_like"] = 0
    return host_like_union


def apply_gates(baits: List[Bait]) -> Tuple[List[Bait], List[dict]]:
    accepted: List[Bait] = []
    ledger: List[dict] = []
    seen_sha: set = set()
    accepted_kmers: List[np.ndarray] = []
    for b in baits:
        st = b.stats
        L = b.length()
        total = st["n_kmers_total"]
        ambig_f = st["n_ambiguous"] / total if total else 1.0
        host_f = st["n_host_like"] / total if total else 1.0
        lcx_f = st["n_low_complexity"] / total if total else 1.0
        reason = None
        if not (bl.BAIT_MIN_BP <= L <= bl.BAIT_MAX_BP):
            reason = f"length_out_of_bounds:{L}"
        elif b.dominated:
            reason = "dominated_nonphage_markers"
        elif ambig_f > bl.AMBIG_FRAC_MAX:
            reason = f"ambiguous_frac>{bl.AMBIG_FRAC_MAX}:{ambig_f:.3f}"
        elif (getattr(b, "expected_behavior", None) != "host_like"
              and host_f > bl.HOSTLIKE_FRAC_MAX):
            reason = f"host_like_frac>{bl.HOSTLIKE_FRAC_MAX}:{host_f:.3f}"
        elif lcx_f > bl.LOWCX_FRAC_MAX:
            reason = f"low_complexity_frac>{bl.LOWCX_FRAC_MAX}:{lcx_f:.3f}"
        else:
            masks = st["_ambig"] | st["_lcx"] | st["_dup"]
            if getattr(b, "expected_behavior", None) != "host_like":
                masks = masks | st["_hostmask"]
            st["n_usable_kmers"] = int((~masks).sum())
            if st["n_usable_kmers"] < bl.USABLE_KMER_MIN:
                reason = (f"usable_kmers<{bl.USABLE_KMER_MIN}:"
                          f"{st['n_usable_kmers']}")
            elif bl.sha256_seq(b.seq) in seen_sha:
                reason = "duplicate_sequence"
            else:
                vals = np.unique(st["_canon"][st["_valid_pos"]])
                max_shared = 0.0
                for ak in accepted_kmers:
                    frac = (np.intersect1d(vals, ak).size / max(1, vals.size))
                    max_shared = max(max_shared, frac)
                st["cross_bait_max_shared_frac"] = round(max_shared, 4)
                if max_shared > bl.CROSS_BAIT_DUP_FRAC:
                    reason = f"cross_bait_kmer_dup:{max_shared:.3f}"
        if reason is None:
            accepted.append(b)
            seen_sha.add(bl.sha256_seq(b.seq))
            accepted_kmers.append(
                np.unique(st["_canon"][st["_valid_pos"]]))
        ledger.append(_ledger_row(b, reason))
    return accepted, ledger


def _ledger_row(b: Bait, reason: Optional[str]) -> dict:
    st = b.stats
    return {
        "bait_id": b.bait_id, "bait_class": b.bait_class,
        "clade_id": b.clade_id, "genome_id": b.genome_id or "",
        "member_prophage_id": b.member_id or "", "contig": b.contig,
        "start": b.start, "end": b.end, "length_bp": b.length(),
        "n_kmers_total": st.get("n_kmers_total", 0),
        "n_ambiguous": st.get("n_ambiguous", 0),
        "n_low_complexity": st.get("n_low_complexity", 0),
        "n_host_like": st.get("n_host_like", 0),
        "n_internal_dup": st.get("n_internal_dup", 0),
        "n_usable_kmers": st.get("n_usable_kmers", ""),
        "seq_sha256": bl.sha256_seq(b.seq) if b.seq else "",
        "status": "accepted" if reason is None else "rejected",
        "reason": reason or "accepted",
    }


# ---------------------------------------------------------------------------
# Host panel
# ---------------------------------------------------------------------------

def choose_host_panel(inp: Inputs, selected: List[dict],
                      member_reps: List[dict],
                      prophage_calls: Dict[str, dict]) -> List[Dict[str, object]]:
    """K-12 first, then <=4 source accessions/clade, <=60 total."""
    clade_genomes = json.load(open(inp.clade_genomes))
    panel: List[Dict[str, object]] = []
    seen: set = set()
    if os.path.exists(inp.host_fasta(K12_ACCESSION)):
        seen.add(K12_ACCESSION)
        panel.append({"accession": K12_ACCESSION, "clade_id": "",
                      "role": "k12_reference"})
    for row in selected:  # selection-rank order
        clade = row["clade_id"]
        rep_accs = []
        for m in member_reps:
            if m["clade_id"] != clade:
                continue
            acc = prophage_calls.get(m["member_prophage_id"], {}).get("genome")
            if acc and acc not in rep_accs:
                rep_accs.append(acc)
        accs = list(rep_accs)
        # remaining member sources via clade source genomes (deterministic)
        cg = clade_genomes.get(clade, [])
        for g in sorted(cg):
            if g not in accs:
                accs.append(g)
            if len(accs) >= 4:
                break
        for acc in accs:
            if acc not in seen and os.path.exists(inp.host_fasta(acc)):
                seen.add(acc)
                panel.append({"accession": acc, "clade_id": clade,
                              "role": "source_host"})
            if len(panel) >= HOST_PANEL_MAX + 1:
                return panel[:HOST_PANEL_MAX + 1]
    return panel


# ---------------------------------------------------------------------------
# Stage 4 — outputs
# ---------------------------------------------------------------------------

MANIFEST_COLS = [
    "bait_id", "bait_class", "clade_id", "genome_id", "tier", "cohort",
    "member_prophage_id", "member_source_accession", "contig", "start", "end",
    "length_bp", "module", "phrog", "gene_annot", "gene_start", "gene_end",
    "partition_a", "partition_b", "breakpoint", "occurrence_a",
    "occurrence_b", "co_occurrence_count", "observed_adjacency_count",
    "adjacency_observed", "partition_a_coords", "partition_b_coords",
    "seq_sha256", "n_kmers_total", "n_kmers_unique_canonical",
    "n_ambiguous", "n_low_complexity", "n_host_like", "n_internal_dup",
    "cross_bait_max_shared_frac", "n_usable_kmers", "gc_frac",
    "expected_behavior", "notes",
]

LEDGER_COLS = ["bait_id", "bait_class", "clade_id", "genome_id",
               "member_prophage_id", "contig", "start", "end", "length_bp",
               "n_kmers_total", "n_ambiguous", "n_low_complexity",
               "n_host_like", "n_internal_dup", "n_usable_kmers",
               "seq_sha256", "status", "reason"]


def manifest_rows(accepted: List[Bait], selected: List[dict],
                  prophage_calls: Dict[str, dict]) -> List[dict]:
    tier = {r["genome_id"]: r["tier"] for r in selected}
    out = []
    for b in accepted:
        st = b.stats
        vals = st["_canon"][st["_valid_pos"]]
        out.append({
            "bait_id": b.bait_id, "bait_class": b.bait_class,
            "clade_id": b.clade_id, "genome_id": b.genome_id or "",
            "tier": tier.get(b.genome_id or "", ""),
            "cohort": COHORT if b.genome_id else "",
            "member_prophage_id": b.member_id or "",
            "member_source_accession": prophage_calls.get(
                b.member_id or "", {}).get("genome", ""),
            "contig": b.contig, "start": b.start, "end": b.end,
            "length_bp": b.length(), "module": b.module or "",
            "phrog": b.phrog or "", "gene_annot": b.annot or "",
            "gene_start": b.gene_lo if b.gene_lo else "",
            "gene_end": b.gene_hi if b.gene_hi else "",
            "partition_a": "" if b.partition_a is None else b.partition_a,
            "partition_b": "" if b.partition_b is None else b.partition_b,
            "breakpoint": "" if b.breakpoint is None else b.breakpoint,
            "occurrence_a": "" if b.occ_a is None else b.occ_a,
            "occurrence_b": "" if b.occ_b is None else b.occ_b,
            "co_occurrence_count": "" if b.co_occ is None else b.co_occ,
            "observed_adjacency_count": "" if b.obs_adj is None else b.obs_adj,
            "adjacency_observed": ("yes" if (b.obs_adj or 0) > 0 else "no")
            if b.bait_class == "junction" else "",
            "partition_a_coords": b.pa_coords or "",
            "partition_b_coords": b.pb_coords or "",
            "seq_sha256": bl.sha256_seq(b.seq),
            "n_kmers_total": st["n_kmers_total"],
            "n_kmers_unique_canonical": int(np.unique(vals).size),
            "n_ambiguous": st["n_ambiguous"],
            "n_low_complexity": st["n_low_complexity"],
            "n_host_like": st.get("n_host_like", 0),
            "n_internal_dup": st["n_internal_dup"],
            "cross_bait_max_shared_frac": st.get(
                "cross_bait_max_shared_frac", ""),
            "n_usable_kmers": st["n_usable_kmers"],
            "gc_frac": st["gc_frac"],
            "expected_behavior": b.expected_behavior or "",
            "notes": b.notes or "",
        })
    return out


def self_audit(records: List[Tuple[str, str]], mrows: List[dict],
               ml_seqs: Dict[str, str], proph_seqs: Dict[str, str],
               ref_seqs: Dict[str, str]) -> None:
    ids = [h for h, _ in records]
    assert len(ids) == len(set(ids)), "duplicate bait ids"
    shas = [bl.sha256_seq(s) for _, s in records]
    assert len(shas) == len(set(shas)), "duplicate bait sequences"
    for (hid, seq), row in zip(records, mrows):
        assert hid == row["bait_id"], f"id mismatch {hid} vs {row['bait_id']}"
        L = len(seq)
        assert bl.BAIT_MIN_BP <= L <= bl.BAIT_MAX_BP, f"{hid} length {L}"
        assert L == int(row["length_bp"]), f"{hid} length col mismatch"
        assert row["seq_sha256"] == bl.sha256_seq(seq), f"{hid} sha mismatch"
        assert row["bait_class"] in BAIT_CLASSES
        if row["bait_class"] == "junction":
            assert row["partition_a"] != "" and row["partition_b"] != ""
            bp = int(row["breakpoint"])
            assert int(row["start"]) <= bp <= int(row["end"]) + 1, \
                f"{hid} does not span its breakpoint"
            # E. coli rule: purely-reconstructed joins only
            assert int(row["observed_adjacency_count"]) == 0, \
                f"{hid} junction is not purely reconstructed"
            assert row["adjacency_observed"] == "no"
        if row["bait_class"] in ("interior_module", "interior_generic",
                                 "junction"):
            src = ml_seqs.get(row["contig"])
            assert src is not None, f"{hid} missing ML source seq"
            assert seq == src[int(row["start"]) - 1:int(row["end"])], \
                f"{hid} sequence does not match source coordinates"
        elif row["bait_class"] == "member_interior":
            src = proph_seqs.get(row["contig"])
            assert src is not None, f"{hid} missing member source seq"
            assert seq == src[int(row["start"]) - 1:int(row["end"])], \
                f"{hid} sequence does not match member coordinates"
        elif row["bait_class"] == "control_public_reference":
            src = ref_seqs.get(row["contig"])
            assert src is not None, f"{hid} missing public reference"
            assert seq == src[int(row["start"]) - 1:int(row["end"])], \
                f"{hid} sequence does not match reference coordinates"
        for col in ("n_kmers_total", "n_usable_kmers", "gc_frac"):
            assert row[col] != "", f"{hid} missing {col}"
        assert 0 <= float(row["gc_frac"]) <= 1
    log(f"SELF-AUDIT OK: {len(records)} baits verified")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_pipeline(args: argparse.Namespace) -> Dict[str, object]:
    repo = os.path.abspath(args.repo_root)
    inp = Inputs(repo, os.path.abspath(args.clades_source),
                 os.path.abspath(args.prophages_fa),
                 os.path.abspath(args.external_root),
                 os.path.abspath(args.host_genomes_dir),
                 os.path.abspath(args.pharokka_cds))
    os.makedirs(inp.work, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)

    log("[select] pilot selection (host-blind, mash-triangle medoids)")
    tri = Triangle(inp.triangle, inp.triangle_ids)
    selected, member_reps, ledger, medoids = select_pilot(inp, tri)
    log(f"  {len(selected)} genomes, {len(member_reps)} member reps, "
        f"{sum(1 for r in ledger if r['eligible'] == 'yes')} eligible "
        f"(A={sum(1 for r in selected if r['tier'] == 'A')}, "
        f"B={sum(1 for r in selected if r['tier'] == 'B')})")

    prophage_calls: Dict[str, dict] = {}
    with open(inp.prophage_calls, newline="") as fh:
        for row in csv.DictReader(fh):
            pid = row.get("prophage_id") or ""
            if pid:
                prophage_calls[pid] = row
    host_intervals, prophage_lens = load_prophage_intervals(inp.prophage_calls)

    log("[stage] clade partitions + deterministic traversal regeneration")
    staged: Dict[str, str] = {}
    trav_verified: Dict[str, bool] = {}
    trav_prov: List[dict] = []
    ml_seqs = load_released_ml(inp, {r["genome_id"] for r in selected})
    verification_log = os.path.join(inp.work, "traversal_verification.log")
    with open(verification_log, "w") as vf:
        vf.write("clade\tbudget\tregenerated_len\treleased_len\t"
                 "byte_identical\n")
        for row in selected:
            clade = row["clade_id"]
            comm = clade.split("_")[0]
            sdir = stage_clade(inp, comm, clade)
            if sdir is None:
                trav_verified[clade] = False
                trav_prov.append({"clade_id": clade,
                                  "status": "traversal_inputs_missing"})
                vf.write(f"{clade}\t\t\t\tinputs_missing\n")
                continue
            budget = traversal_budget(sdir)
            try:
                tjson = regenerate_traversal(inp, comm, clade, sdir)
                ml_fa = tjson.replace(".traversal.json", ".ml.fa")
                ok, rlen, elen = verify_traversal(
                    ml_fa, ml_seqs[row["genome_id"]])
            except RuntimeError as exc:
                ok, rlen, elen = False, "ERR", \
                    len(ml_seqs[row["genome_id"]])
                vf.write(f"{clade}\t{budget}\tERR\t{elen}\t{exc}\n")
            trav_verified[clade] = ok
            trav_prov.append({
                "clade_id": clade, "status": "verified" if ok else
                "traversal_regeneration_mismatch",
                "budget_bp": budget,
                "traversal_json_sha256": sha256_file(tjson) if ok else "",
                "partitions_bed_sha256": sha256_file(
                    os.path.join(sdir, "partitions.bed")),
                "regenerated_len_bp": rlen if ok else "",
                "released_len_bp": elen,
            })
            if ok:
                vf.write(f"{clade}\t{budget}\t{rlen}\t{elen}\tTrue\n")
                staged[clade] = sdir
    n_ok = sum(trav_verified.values())
    log(f"  traversal verified byte-identical for {n_ok}/{len(selected)} "
        f"clades; junctions disabled elsewhere")
    if n_ok == 0:
        log("  STOP-AND-REPORT CONDITION: no selected clade reproduced its "
            "released ML genome; emitting a junction-free panel")

    log("[design] bait candidates")
    wanted_members = {m["member_prophage_id"] for m in member_reps
                      if m["is_medoid"] == "yes"}
    proph_seqs: Dict[str, str] = {}
    if os.path.exists(inp.member_fa):
        for hid, seq in iter_fasta(inp.member_fa):
            if hid in wanted_members:
                proph_seqs[hid] = seq.upper()
    missing_members = wanted_members - set(proph_seqs)
    if missing_members:
        with open(inp.member_fa, "a") as out:
            for hid, seq in iter_fasta(inp.prophages_fa):
                if hid in missing_members:
                    out.write(f">{hid}\n")
                    s = seq.upper()
                    for i in range(0, len(s), 80):
                        out.write(s[i:i + 80] + "\n")
                    proph_seqs[hid] = s
        for hid in missing_members:
            prophage_lens.setdefault(
                hid, len(proph_seqs.get(hid, "")))
    cds_rows = load_cds_rows(inp.pharokka_cds, {r["genome_id"]
                                                for r in selected})
    baits, pre_rej = generate_baits(inp, selected, member_reps, ml_seqs,
                                    cds_rows, proph_seqs, staged,
                                    trav_verified)
    log(f"  {len(baits)} real baits; adding controls")

    ref_seqs: Dict[str, str] = {}
    for hb in make_host_control_baits(inp, selected, member_reps,
                                      host_intervals, prophage_calls):
        compute_masks(hb)
        baits.append(hb)
    for pb in make_public_reference_baits(inp):
        ref_seqs[pb.contig] = next(
            iter_fasta(inp.public_refs[pb.contig]))[1].upper()
        compute_masks(pb)
        baits.append(pb)
    baits.sort(key=lambda b: b.priority)

    log("[mask] 31-mer statistics + host scan")
    for b in baits:
        compute_masks(b)
    panel = choose_host_panel(inp, selected, member_reps, prophage_calls)
    log(f"[host-panel] {len(panel)} host assemblies (prophages masked), "
        f"K-12 first")
    host_like_union = host_scan(inp, baits, panel, host_intervals,
                                os.path.join(inp.work, "host_scan.log"))
    accepted, rej_ledger = apply_gates(baits)
    log(f"[gates] {len(accepted)} accepted, "
        f"{len(rej_ledger) - len(accepted)} rejected")

    for sb in make_shuffled_controls(accepted, host_like_union):
        accepted.append(sb)
        rej_ledger.append(_ledger_row(sb, None))

    log("[report] outputs")
    mrows = manifest_rows(accepted, selected, prophage_calls)
    records = [(b.bait_id, b.seq) for b in accepted]
    bl.write_fasta(records, os.path.join(args.out_dir, "baits.fa"))

    ps_rows = []
    for row in selected:
        clade = row["clade_id"]
        reps = [m for m in member_reps if m["clade_id"] == clade]
        accs = sorted({prophage_calls.get(m["member_prophage_id"], {})
                       .get("genome", "") for m in reps} - {""})
        ps_rows.append({**row,
                        "member_reps": ";".join(
                            m["member_prophage_id"] for m in reps),
                        "source_accessions": ";".join(accs),
                        "traversal_verified": "yes" if trav_verified[clade]
                        else "no"})
    write_tsv(os.path.join(args.out_dir, "pilot_selection.tsv"), ps_rows,
              ["selection_rank", "clade_id", "genome_id", "tier", "n_members",
               "n_partitions", "length_bp", "release_flags",
               "functional_qc_flag", "medoid_member", "member_reps",
               "source_accessions", "traversal_verified", "selection_basis"])
    write_tsv(os.path.join(args.out_dir, "bait_manifest.tsv"), mrows,
              MANIFEST_COLS)
    write_tsv(os.path.join(args.out_dir, "rejection_ledger.tsv"),
              sorted(rej_ledger, key=lambda r: (r["status"] != "rejected",
                                                r["bait_id"])), LEDGER_COLS)

    fs_rows = []
    for cls in BAIT_CLASSES:
        rows = [r for r in mrows if r["bait_class"] == cls]
        if not rows:
            continue
        usable = [int(r["n_usable_kmers"]) for r in rows]
        hl = [int(r["n_host_like"]) / max(1, int(r["n_kmers_total"]))
              for r in rows]
        fs_rows.append({
            "bait_class": cls, "n_baits": len(rows),
            "median_usable_kmers": int(np.median(usable)),
            "min_usable_kmers": int(min(usable)),
            "median_host_like_frac": round(float(np.median(hl)), 4),
            "max_host_like_frac": round(float(max(hl)), 4),
        })
    write_tsv(os.path.join(args.out_dir, "filter_stats.tsv"), fs_rows,
              ["bait_class", "n_baits", "median_usable_kmers",
               "min_usable_kmers", "median_host_like_frac",
               "max_host_like_frac"])

    # provenance — deterministic run id over input checksums + params
    driver = os.path.abspath(__file__)
    checksums = {
        "ml_fasta": sha256_file(inp.ml_fasta),
        "release_manifest": sha256_file(inp.release_manifest),
        "release_qc": sha256_file(inp.release_qc),
        "functional_qc": sha256_file(inp.functional_qc),
        "clade_genomes": sha256_file(inp.clade_genomes),
        "triangle_ids": sha256_file(inp.triangle_ids),
        "prophage_calls": sha256_file(inp.prophage_calls),
        "pharokka_cds": sha256_file(inp.pharokka_cds),
        "baitlib": sha256_file(os.path.join(
            _HERE, "..", "..", "..", "ntm", "v2", "bait", "scripts",
            "baitlib.py")),
        "traverse_script": sha256_file(inp.traverse_script),
    }
    for acc, path in sorted(inp.public_refs.items()):
        checksums[f"public_ref_{acc}"] = sha256_file(path)
    for comm_file in sorted(inp.community_tight_clades().keys()):
        p = os.path.join(repo, "research", "clades", comm_file,
                         "tight_clades.json")
        checksums[f"tight_clades_{comm_file}"] = sha256_file(p)
    for prov in trav_prov:
        if prov["partitions_bed_sha256"]:
            checksums[f"staged_{prov['clade_id']}/partitions.bed"] = \
                prov["partitions_bed_sha256"]
        if prov.get("traversal_json_sha256"):
            checksums[f"traversal_{prov['clade_id']}"] = \
                prov["traversal_json_sha256"]
    host_panel_prov = [{
        "accession": str(e["accession"]),
        "fasta": inp.host_fasta(str(e["accession"])),
        "fasta_sha256": sha256_file(inp.host_fasta(str(e["accession"]))),
        "prophage_intervals_masked": sum(
            len(v) for v in host_intervals.get(str(e["accession"]), {}).values()),
        "clade_id": e["clade_id"], "role": e["role"],
    } for e in panel]
    params = {
        "K": bl.K, "SELECT_N_A": SELECT_N_A, "SELECT_N_B": SELECT_N_B,
        "MIN_GENOME_BP": MIN_GENOME_BP, "MEMBER_MIN_BP": MEMBER_MIN_BP,
        "MODULE_W": MODULE_W, "JUNCTION_W": JUNCTION_W, "MEMBER_W": MEMBER_W,
        "HOSTCTL_W": HOSTCTL_W, "POSCON_W": POSCON_W, "MODULE_CAP": MODULE_CAP,
        "JUNCTION_CAP": JUNCTION_CAP,
        "N_SHUFFLED_CONTROLS": N_SHUFFLED_CONTROLS,
        "N_HOST_CONTROLS": N_HOST_CONTROLS, "SHUFFLED_SEED": SHUFFLED_SEED,
        "TRAV_N_SAMPLES": TRAV_N_SAMPLES, "TRAV_SEED": TRAV_SEED,
        "TRAV_MIN_BUDGET": TRAV_MIN_BUDGET,
        "TRAV_MEMBER_MULT": TRAV_MEMBER_MULT,
        "tier_B_allowed_flags": sorted(TIER_B_ALLOWED_FLAGS),
        "release_flag_exclude": list(RELEASE_FLAG_EXCLUDE),
        "host_panel_max": HOST_PANEL_MAX, "k12_accession": K12_ACCESSION,
        "entropy_min_bits": bl.ENTROPY_MIN_BITS,
        "homopolymer_run": bl.HOMOPOLYMER_RUN,
        "usable_kmer_min": bl.USABLE_KMER_MIN,
        "ambig_frac_max": bl.AMBIG_FRAC_MAX,
        "hostlike_frac_max": bl.HOSTLIKE_FRAC_MAX,
        "lowcx_frac_max": bl.LOWCX_FRAC_MAX,
        "cross_bait_dup_frac": bl.CROSS_BAIT_DUP_FRAC,
    }
    run_id = hashlib.sha256(json.dumps(
        {"inputs": checksums, "params": params,
         "script_sha256": sha256_file(driver)},
        sort_keys=True).encode()).hexdigest()[:16]
    provenance = {
        "run_id": run_id,
        "design": "research/ecoli_bait/DESIGN.md",
        "task": "design-e-coli",
        "inputs_sha256": checksums,
        "host_panel": host_panel_prov,
        "traversal_regeneration": trav_prov,
        "public_references": [
            {"accession": acc, "label": label, "length_bp": exp_len,
             "sha256": checksums[f"public_ref_{acc}"]}
            for acc, (label, exp_len) in sorted(PUBLIC_REFS.items())],
        "params": params,
        "script_sha256": sha256_file(driver),
        "counts": {
            "selected_genomes": len(selected),
            "member_reps": len(member_reps),
            "baits_accepted": len(accepted),
            "baits_rejected": len(rej_ledger) - len(accepted),
            "eligibility_ledger_rows": len(ledger),
            "host_panel_genomes": len(panel),
            "traversal_verified_clades": n_ok,
        },
        "environment": {"python": sys.version.split()[0],
                        "numpy": np.__version__},
        "junction_rule": ("purely-reconstructed joins only "
                          "(observed_adjacency_count == 0); regeneration "
                          "verified byte-identical to the frozen release "
                          "per clade; unverified clades contribute zero "
                          "junction baits"),
    }
    with open(os.path.join(args.out_dir, "provenance.json"), "w") as fh:
        json.dump(provenance, fh, indent=2, sort_keys=True)
        fh.write("\n")

    tier_counts: Dict[str, int] = {}
    for r in selected:
        tier_counts[r["tier"]] = tier_counts.get(r["tier"], 0) + 1
    cls_counts: Dict[str, int] = {}
    for r in mrows:
        cls_counts[r["bait_class"]] = cls_counts[r["bait_class"]] + 1
    n_recon = sum(1 for r in mrows if r["bait_class"] == "junction")
    with open(os.path.join(args.out_dir, "README.md"), "w") as fh:
        fh.write(f"""# E. coli prophage bait & junction panel

Generated by `../scripts/design_ecoli_baits.py` per the preregistered
`../DESIGN.md` (mirrors the NTM v1 panel; shared baitlib). Deterministic
run id: `{run_id}`.

- Pilot selection: {len(selected)} ML genomes (tiers: {tier_counts});
  host/phylogroup labels are provenance only, never selection features.
- Baits accepted: {len(mrows)} by class: {cls_counts}.
- Junction baits are **purely reconstructed** (observed_adjacency_count == 0)
  with co-occurrence ranking; traversal provenance regenerated deterministically
  and verified byte-identical to the frozen release for {n_ok}/{len(selected)}
  selected clades (unverified clades contribute no junctions).
- Positive controls: exact mid-genome slices of lambda NC_001416.1 and HK97
  NC_002167.1 (reference sha256 in provenance + `public_refs/`).
- Negative controls: {cls_counts.get('control_shuffled_negative', 0)} Markov
  composition-matched negatives (expected no hits) and
  {cls_counts.get('control_host_negative', 0)} E. coli host-derived controls
  (expected host-like; validates the host mask).

Files: `baits.fa`, `pilot_selection.tsv`, `bait_manifest.tsv`,
`rejection_ledger.tsv`, `filter_stats.tsv`, `provenance.json`,
`public_refs/`. Reruns are byte-identical (`--verify-determinism`).
""")
    with open(os.path.join(inp.work, "host_panel.json"), "w") as fh:
        json.dump(host_panel_prov, fh, indent=2, sort_keys=True)

    self_audit(records, mrows, ml_seqs, proph_seqs, ref_seqs)
    log(f"[done] panel in {args.out_dir} (run_id {run_id})")
    return {"provenance": provenance, "counts": provenance["counts"]}


def verify_determinism(args: argparse.Namespace) -> int:
    """Run the whole pipeline twice into temp dirs; compare committed bytes."""
    hashes: List[Dict[str, str]] = []
    for _ in range(2):
        with tempfile.TemporaryDirectory() as td:
            ext = os.path.join(td, "external")
            out = os.path.join(td, "panel")
            sub = argparse.Namespace(**{**vars(args), "external_root": ext,
                                        "out_dir": out})
            run_pipeline(sub)  # fresh staging + traversal regen each time
            files = [f for f in sorted(os.listdir(out))
                     if os.path.isfile(os.path.join(out, f))]
            hashes.append({f: sha256_file(os.path.join(out, f))
                           for f in files})
    if hashes[0] != hashes[1]:
        for k in sorted(set(hashes[0]) | set(hashes[1])):
            if hashes[0].get(k) != hashes[1].get(k):
                log(f"DETERMINISM FAIL: {k}")
        return 1
    log("DETERMINISM OK: both runs byte-identical")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    here = os.path.dirname(os.path.abspath(__file__))
    default_repo = os.path.abspath(os.path.join(here, "..", "..", ".."))
    ap.add_argument("--repo-root", default=default_repo)
    ap.add_argument("--clades-source",
                    default="/home/erikg/phind/research/clades",
                    help="main-tree research/clades (untracked partition data)")
    ap.add_argument("--prophages-fa",
                    default="/home/erikg/phind/prophage_homology_survey/"
                            "full_prophages.fa")
    ap.add_argument("--external-root",
                    default="/mnt/nvme3n1/erikg/phind-genome-work/ecoli_bait")
    ap.add_argument("--host-genomes-dir",
                    default="/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/"
                            "v1/26k/canonical_objects")
    ap.add_argument("--pharokka-cds",
                    default="/mnt/nvme3n1/erikg/phind-genome-work/annotation/"
                            "pharokka_out/pharokka_cds_final_merged_output.tsv")
    ap.add_argument("--out-dir",
                    default=os.path.join(default_repo, "research", "ecoli_bait",
                                         "panel"))
    ap.add_argument("--verify-determinism", action="store_true")
    args = ap.parse_args(argv)

    if args.verify_determinism:
        return verify_determinism(args)
    run_pipeline(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
