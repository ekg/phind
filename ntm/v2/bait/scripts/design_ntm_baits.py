#!/usr/bin/env python3
"""NTM prophage bait & junction panel driver (preregistered: ../DESIGN.md).

Deterministic end-to-end design: pilot selection (diversity-aware, host-blind)
→ bait generation (interior module / junction / member / controls) → 31-mer
masking with acceptance gates → committed panel outputs + external ledgers.

Run:
  python3 design_ntm_baits.py \
      --v2-root /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2 \
      --external-root /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/bait \
      --out-dir ntm/v2/bait/panel

Committed outputs are byte-deterministic (no wall-clock inside them).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import tempfile
from typing import Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import baitlib as bl  # noqa: E402

SELECT_N_A = 15
SELECT_N_B = 5
MIN_GENOME_BP = 5000
MEMBER_MIN_BP = 600
MODULE_W = 1200
JUNCTION_W = 1000
MEMBER_W = 1200
HOSTCTL_W = 1000
MODULE_CAP = 3
JUNCTION_CAP = 3
N_SHUFFLED_CONTROLS = 3
N_HOST_CONTROLS = 3
SHUFFLED_SEED = 42
SHARED_KMER_CONTROL_MAX = 5
TIERS = ("A_strong_candidate", "B_moderate_candidate")
INVALID = np.uint64(1) << np.uint64(62)

BAIT_CLASSES = (
    "interior_module", "interior_generic", "junction", "member_interior",
    "control_shuffled_negative", "control_host_negative",
    "control_public_reference",
)


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
    """Stream (id, sequence) from a plain FASTA."""
    with open(path) as fh:
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


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

class Inputs:
    def __init__(self, v2_root: str):
        self.root = v2_root
        p = lambda *a: os.path.join(v2_root, *a)  # noqa: E731
        self.functional_qc = p("annotation/report/per_genome_functional_qc.tsv")
        self.release_manifest = p("ml/release_manifest.tsv")
        self.ml_fasta = p("ml/all_ntm2_ml_phage_genomes.fa")
        self.distances = p("clades/0/distances.npz")
        self.members_json = p("clades/0/members.json")
        self.tight_clades = p("clades/0/tight_clades.json")
        self.prophages_fa = p("full_prophages.fa")
        self.master_manifest = p(
            "inputs/NTM_QC_passed_prophage_master_manifest.tsv")
        self.pharokka_cds = p(
            "annotation/pharokka_out/pharokka_cds_final_merged_output.tsv")
        self.host_clades_tsv = p("host_clades/host_clades.tsv")
        self.genomes_dir = p("genomes/canonical_objects")

    def clade_dir(self, clade_id: str) -> str:
        return os.path.join(self.root, "clades", clade_id)

    def host_fasta(self, accession: str) -> str:
        return os.path.join(self.genomes_dir, accession,
                            f"{accession}.pansn.fa.gz")

    def checksums(self, paths: Dict[str, str]) -> Dict[str, str]:
        return {k: sha256_file(v) for k, v in sorted(paths.items())}


def load_master_prophage_rows(path: str) -> Dict[str, dict]:
    """prophage_id -> row (master manifest; only rows with a prophage call)."""
    out: Dict[str, dict] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pid = row.get("prophage_id") or ""
            if pid:
                out[pid] = row
    return out


def load_host_intervals(master: Dict[str, dict]
                        ) -> Dict[str, Dict[str, List[Tuple[int, int]]]]:
    """accession -> contig -> sorted [(start, end)] prophage intervals (1-based)."""
    out: Dict[str, Dict[str, List[Tuple[int, int]]]] = {}
    for row in master.values():
        try:
            s = int(round(float(row["prophage_start"])))
            e = int(round(float(row["prophage_end"])))
        except (ValueError, TypeError):
            continue
        out.setdefault(row["accession"], {}).setdefault(
            row["prophage_contig"], []).append((s, e))
    for acc in out:
        for contig in out[acc]:
            out[acc][contig].sort()
    return out


def load_host_clades(path: str) -> Dict[str, Tuple[str, str]]:
    """accession -> (host_clade_id, species) — provenance only, never features."""
    out: Dict[str, Tuple[str, str]] = {}
    for row in read_tsv(path):
        out[row["accession"]] = (row["host_clade_id"], row.get("species", ""))
    return out


def load_selected_ml_sequences(path: str, wanted: set) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    for hid, seq in iter_fasta(path):
        if hid in wanted:
            seqs[hid] = seq.upper()
    return seqs


def load_cds_rows(path: str, contigs: set) -> Dict[str, List[dict]]:
    rows: Dict[str, List[dict]] = {c: [] for c in contigs}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            c = row.get("contig")
            if c in rows:
                rows[c].append(row)
    return rows


# ---------------------------------------------------------------------------
# Stage 1 — pilot selection (deterministic, diversity-aware, host-blind)
# ---------------------------------------------------------------------------

def select_pilot(inp: Inputs, work: str
                 ) -> Tuple[List[dict], List[dict], List[dict], Dict[str, str]]:
    """Returns (selected genomes, member reps, eligibility ledger, medoids)."""
    qc = read_tsv(inp.functional_qc)
    rel = {r["clade_id"]: r for r in read_tsv(inp.release_manifest)}
    tight = json.load(open(inp.tight_clades))
    z = np.load(inp.distances, allow_pickle=True)
    D = z["D"]
    sub_index = {m: i for i, m in enumerate(list(z["members"]))}

    ledger: List[dict] = []
    pools: Dict[str, List[str]] = {t: [] for t in TIERS}
    tier_by_clade: Dict[str, str] = {}

    for row in qc:
        clade = row["clade_id"]
        reason = None
        if row["cohort"] != "ntm2_ml_reconstructed":
            reason = f"cohort={row['cohort']}"
        elif row["candidate_functionality"] not in TIERS:
            reason = f"tier={row['candidate_functionality']}"
        elif row["flag"] != "ok":
            reason = f"functional_qc_flag={row['flag']}"
        elif clade not in rel:
            reason = "not_in_release_manifest"
        else:
            r = rel[clade]
            if r["status"] != "ml":
                reason = f"release_status={r['status']}"
            elif int(r["length_bp"]) < MIN_GENOME_BP:
                reason = f"release_length<{MIN_GENOME_BP}"
            elif "too_short" in (r["flags"] or ""):
                reason = f"release_flags={r['flags']}"
            elif not os.path.exists(os.path.join(
                    inp.clade_dir(clade), "ml.traversal.json")):
                reason = "missing_ml_traversal"
            elif int(r["n_members"]) < 2:
                reason = "n_members<2"
        ledger.append({
            "genome_id": row["genome_id"], "clade_id": clade,
            "cohort": row["cohort"], "tier": row["candidate_functionality"],
            "functional_qc_flag": row["flag"],
            "eligible": "yes" if reason is None else "no",
            "reason": reason or "eligible",
        })
        if reason is None:
            pools[row["candidate_functionality"]].append(clade)
            tier_by_clade[clade] = row["candidate_functionality"]

    # medoids over prophage mash space (host labels never enter here)
    medoids: Dict[str, str] = {}
    for tier in TIERS:
        for clade in pools[tier]:
            mem = sorted(m for m in tight[clade] if m in sub_index)
            if len(mem) < 2:
                medoids[clade] = mem[0] if mem else ""
                continue
            idx = np.array([sub_index[m] for m in mem])
            sub = D[np.ix_(idx, idx)].astype(np.float64)
            sub[np.isnan(sub)] = 1.0
            medoids[clade] = mem[int(np.argmin(sub.mean(axis=1)))]

    cand_all = sorted(set(pools[TIERS[0]]) | set(pools[TIERS[1]]))
    dist: Dict[Tuple[str, str], float] = {}
    for a in cand_all:
        for b in cand_all:
            if a != b:
                v = D[sub_index[medoids[a]], sub_index[medoids[b]]]
                dist[(a, b)] = 1.0 if np.isnan(v) else float(v)

    a_pool = sorted(pools[TIERS[0]])
    b_pool = sorted(pools[TIERS[1]])
    selected: List[str] = []
    if a_pool:
        selected += bl.greedy_farthest_point(a_pool, dist, a_pool[0],
                                             SELECT_N_A)
    quota_b = SELECT_N_A + SELECT_N_B - len(selected)
    if quota_b > 0 and b_pool:
        selected += bl.greedy_farthest_point(b_pool, dist, b_pool[0],
                                             quota_b)
    if len(selected) < SELECT_N_A + SELECT_N_B:
        rest = [c for c in a_pool + b_pool if c not in selected]
        selected += rest[:SELECT_N_A + SELECT_N_B - len(selected)]

    # member representatives: medoid + farthest-from-medoid (>= MEMBER_MIN_BP)
    proph_rows = load_master_prophage_rows(inp.master_manifest)

    def plen(pid: str) -> int:
        try:
            return int(round(float(proph_rows[pid]["prophage_length_bp"])))
        except (KeyError, ValueError, TypeError):
            return 0

    member_reps: List[dict] = []
    for rank, clade in enumerate(selected, 1):
        mem = sorted(m for m in tight[clade] if m in sub_index)
        if not mem:
            continue
        medoid = medoids[clade]
        reps: List[Tuple[str, str]] = []
        if medoid and plen(medoid) >= MEMBER_MIN_BP:
            reps.append((medoid, "medoid"))
        if len(mem) >= 2:
            i = sub_index[medoid]
            order = sorted(
                ((1.0 if np.isnan(D[i, sub_index[m]]) else
                  float(D[i, sub_index[m]]), m) for m in mem if m != medoid),
                key=lambda t: (-t[0], t[1]))
            for _, m in order:
                if len(reps) >= 2:
                    break
                if plen(m) >= MEMBER_MIN_BP and m not in (r[0] for r in reps):
                    reps.append((m, "farthest_from_medoid"))
                    break
        if not reps:  # fallback: longest member
            best = max(mem, key=lambda m: (plen(m), m))
            if plen(best) >= MEMBER_MIN_BP:
                reps.append((best, "longest_member_fallback"))
        for rep, role in reps:
            member_reps.append({
                "selection_rank": rank, "clade_id": clade,
                "member_prophage_id": rep, "member_length_bp": plen(rep),
                "is_medoid": "yes" if rep == medoid else "no", "role": role,
            })

    selected_rows = []
    for rank, clade in enumerate(selected, 1):
        r = rel[clade]
        selected_rows.append({
            "selection_rank": rank, "clade_id": clade,
            "genome_id": f"ntm2_{clade}_ML", "tier": tier_by_clade[clade],
            "n_members": r["n_members"], "length_bp": r["length_bp"],
            "release_flags": r["flags"] or "",
            "medoid_member": medoids[clade],
            "selection_basis": "greedy_farthest_point_prophage_mash_"
                               "space_host_blind",
        })
    write_tsv(os.path.join(work, "eligibility_ledger.tsv"), ledger,
              ["genome_id", "clade_id", "cohort", "tier",
               "functional_qc_flag", "eligible", "reason"])
    return selected_rows, member_reps, ledger, medoids


def choose_host_panel(inp: Inputs, selected: List[dict],
                      member_reps: List[dict]
                      ) -> List[Dict[str, object]]:
    """Deterministic host panel: rep/medoid sources, <=4/clade, <=60 total."""
    proph_rows = load_master_prophage_rows(inp.master_manifest)
    by_clade: Dict[str, List[str]] = {}
    for mr in member_reps:
        acc = proph_rows.get(mr["member_prophage_id"], {}).get("accession")
        if acc:
            by_clade.setdefault(mr["clade_id"], []).append(acc)
    tight = json.load(open(inp.tight_clades))
    panel: List[Dict[str, object]] = []
    seen: set = set()
    for row in selected:  # selection-rank order
        clade = row["clade_id"]
        accs = list(dict.fromkeys(by_clade.get(clade, [])))  # dedupe, keep order
        for pid in sorted(tight[clade]):
            acc = proph_rows.get(pid, {}).get("accession")
            if acc and acc not in accs:
                accs.append(acc)
            if len(accs) >= 4:
                break
        for acc in accs:
            if acc not in seen and os.path.exists(inp.host_fasta(acc)):
                seen.add(acc)
                panel.append({"accession": acc, "clade_id": clade,
                              "role": "source_host"})
            if len(panel) >= 60:
                return panel
    return panel


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


def junction_candidates(inp: Inputs, clade: str) -> List[dict]:
    """Ranked junction candidates: both partitions + exact breakpoint."""
    cdir = inp.clade_dir(clade)
    trav = json.load(open(os.path.join(cdir, "ml.traversal.json")))
    parts = trav["partitions"]
    pids = trav["genome"]["pids"]
    spans: Dict[int, Tuple[int, int]] = {}
    pos = 0
    for pid in pids:
        L = int(parts[str(pid)]["representative_len"])
        spans[pid] = (pos + 1, pos + L)
        pos += L

    per_member: Dict[str, List[Tuple[int, int, int]]] = {}
    with open(os.path.join(cdir, "partitions.bed")) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 4:
                pid = int(f[3])
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
        cands.append({
            "partition_a": a, "partition_b": b, "breakpoint": spans[b][0],
            "occ_a": int(pa["occurrence"]), "occ_b": int(pb["occurrence"]),
            "co_occurrence_count": co_occ, "observed_adjacency_count": obs_adj,
            "pa_coords": f"{spans[a][0]}-{spans[a][1]}",
            "pb_coords": f"{spans[b][0]}-{spans[b][1]}",
            "support": min(int(pa["occurrence"]), int(pb["occurrence"])),
        })
    cands.sort(key=lambda c: (-c["support"], c["partition_a"],
                              c["partition_b"]))
    return cands


def generate_baits(inp: Inputs, selected: List[dict], member_reps: List[dict],
                   ml_seqs: Dict[str, str], cds_rows: Dict[str, List[dict]],
                   proph_seqs: Dict[str, str]
                   ) -> List[Bait]:
    baits: List[Bait] = []
    counters: Dict[Tuple[str, str], int] = {}

    def bait_id(clade: str, cls: str) -> str:
        key = (clade, cls)
        counters[key] = counters.get(key, 0) + 1
        return f"NTMBAIT_{clade}_{cls.upper()}_{counters[key]:02d}"

    for row in selected:
        clade = row["clade_id"]
        genome = row["genome_id"]
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
        gi = 0
        for frac in (0.25, 0.55, 0.85):  # generic fill to 3 interior baits
            if n_module >= MODULE_CAP:
                break
            s, e = bl.centered_window(glen, int(glen * frac) + 1, MODULE_W)
            if any(bl.windows_overlap((s, e), u) for u in used):
                continue
            baits.append(new_bait(
                bait_id=bait_id(clade, "interior_generic"),
                bait_class="interior_generic", clade_id=clade, genome_id=genome,
                member_id="", contig=genome, start=s, end=e, seq=seq[s - 1:e],
                module="generic_interior",
                dominated=bl.cds_composition_flags(rows, s, e)[0],
                priority=(1, row["selection_rank"], gi)))
            used.append((s, e))
            n_module += 1
            gi += 1
        nj = 0
        try:
            trav = json.load(open(os.path.join(
                inp.clade_dir(clade), "ml.traversal.json")))
            trav_len = int(trav["genome"]["length_bp"])
            jcands = junction_candidates(inp, clade)
        except (OSError, KeyError, ValueError, json.JSONDecodeError):
            trav_len, jcands = None, []
        if trav_len == glen:
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
        for mr in [m for m in member_reps if m["clade_id"] == clade]:
            pid = mr["member_prophage_id"]
            mseq = proph_seqs.get(pid)
            if not mseq:
                continue
            w = min(MEMBER_W, len(mseq))
            if w < bl.BAIT_MIN_BP:
                continue
            s, e = bl.centered_window(len(mseq), len(mseq) // 2 + 1, w)
            baits.append(new_bait(
                bait_id=bait_id(clade, "member_interior"),
                bait_class="member_interior", clade_id=clade,
                genome_id=genome, member_id=pid, contig=pid, start=s, end=e,
                seq=mseq[s - 1:e], module="member_midpoint", dominated=False,
                priority=(3, row["selection_rank"],
                          0 if mr["is_medoid"] == "yes" else 1)))
    baits.sort(key=lambda b: b.priority)
    return baits


# ---------------------------------------------------------------------------
# Stage 3 — masking + gates
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
    """Single host pass marking host-like canonical 31-mers per bait position.

    Returns the sorted union of host-like canonical kmers (for the Markov
    negative controls, which are generated after this scan).
    """
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
            # usable kmers: unmasked, unique. For host-derived negative
            # controls the host mask is exempt — being host-like is their
            # purpose (they validate the mask); hostmask is excluded here.
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
                vals = np.unique(b.stats["_canon"][b.stats["_valid_pos"]])
                max_shared = 0.0
                for ak in accepted_kmers:
                    frac = (np.intersect1d(vals, ak).size /
                            max(1, vals.size))
                    max_shared = max(max_shared, frac)
                st["cross_bait_max_shared_frac"] = round(max_shared, 4)
                if max_shared > bl.CROSS_BAIT_DUP_FRAC:
                    reason = f"cross_bait_kmer_dup:{max_shared:.3f}"
        if reason is None:
            accepted.append(b)
            seen_sha.add(bl.sha256_seq(b.seq))
            accepted_kmers.append(
                np.unique(b.stats["_canon"][b.stats["_valid_pos"]]))
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
# Controls
# ---------------------------------------------------------------------------

def make_host_control_baits(
        inp: Inputs, selected: List[dict], member_reps: List[dict],
        proph_rows: Dict[str, dict],
        intervals: Dict[str, Dict[str, List[Tuple[int, int]]]]) -> List[Bait]:
    """1-kb windows from prophage-free stretches of source-host contigs."""
    baits: List[Bait] = []
    used_acc: set = set()
    for row in selected:
        clade = row["clade_id"]
        reps = [m for m in member_reps if m["clade_id"] == clade]
        if not reps:
            continue
        acc = proph_rows.get(reps[0]["member_prophage_id"], {}).get("accession")
        if not acc or acc in used_acc or not os.path.exists(inp.host_fasta(acc)):
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
            bait_id=f"NTMBAIT_HOST_CONTROL_{len(baits) + 1:02d}",
            bait_class="control_host_negative", clade_id=clade, genome_id="",
            member_id=reps[0]["member_prophage_id"], contig=contig,
            start=s, end=e, seq=seq[s - 1:e], module="host_control",
            dominated=False, expected_behavior="host_like",
            notes=("host-derived negative/filter-validation control: "
                   "prophage-free stretch of a source-host contig"),
            priority=(4, len(selected), len(baits) + 1)))
        used_acc.add(acc)
        if len(baits) >= N_HOST_CONTROLS:
            break
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
            bait_id=f"NTMBAIT_SHUF_CONTROL_{len(baits) + 1:02d}",
            bait_class="control_shuffled_negative", clade_id=b.clade_id,
            genome_id="", member_id="", contig=b.contig, start=b.start,
            end=b.end, seq=cand, module="markov_negative", dominated=False,
            expected_behavior="no_hits",
            notes=(f"source_bait={b.bait_id};seed={seed};"
                   f"shared_panel_kmers={shared};host_like_kmers={hl};"
                   f"mono_l1={dev['mono_l1']:.4f};"
                   f"dinuc_l1={dev['dinuc_l1']:.4f};redraws={tries}"),
            priority=(5, len(baits) + 1, 0))
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
                  proph_rows: Dict[str, dict]) -> List[dict]:
    tier = {r["genome_id"]: r["tier"] for r in selected}
    out = []
    for b in accepted:
        st = b.stats
        vals = st["_canon"][st["_valid_pos"]]
        out.append({
            "bait_id": b.bait_id, "bait_class": b.bait_class,
            "clade_id": b.clade_id, "genome_id": b.genome_id or "",
            "tier": tier.get(b.genome_id or "", ""),
            "cohort": "ntm2_ml_reconstructed" if b.genome_id else "",
            "member_prophage_id": b.member_id or "",
            "member_source_accession": proph_rows.get(
                b.member_id or "", {}).get("accession", ""),
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
               ml_seqs: Dict[str, str], proph_seqs: Dict[str, str]) -> None:
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
        for col in ("n_kmers_total", "n_usable_kmers", "gc_frac"):
            assert row[col] != "", f"{hid} missing {col}"
        assert 0 <= float(row["gc_frac"]) <= 1
    log(f"SELF-AUDIT OK: {len(records)} baits verified")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--v2-root", required=True)
    ap.add_argument("--external-root", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--verify-determinism", action="store_true")
    args = ap.parse_args(argv)

    if args.verify_determinism:
        return verify_determinism(args)

    inp = Inputs(args.v2_root)
    work = os.path.join(args.external_root, "work")
    os.makedirs(work, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)

    log("[select] pilot selection")
    selected, member_reps, ledger, medoids = select_pilot(inp, work)
    log(f"  {len(selected)} genomes, {len(member_reps)} member reps, "
        f"{sum(1 for r in ledger if r['eligible'] == 'yes')} eligible")

    proph_rows = load_master_prophage_rows(inp.master_manifest)
    host_intervals = load_host_intervals(proph_rows)
    panel = choose_host_panel(inp, selected, member_reps)
    log(f"[host-panel] {len(panel)} host assemblies (prophages masked)")

    log("[design] bait candidates")
    wanted_ml = {r["genome_id"] for r in selected}
    ml_seqs = load_selected_ml_sequences(inp.ml_fasta, wanted_ml)
    proph_seqs = {}
    for hid, seq in iter_fasta(inp.prophages_fa):
        if hid in {m["member_prophage_id"] for m in member_reps}:
            proph_seqs[hid] = seq.upper()
    cds_rows = load_cds_rows(inp.pharokka_cds, wanted_ml)
    baits = generate_baits(inp, selected, member_reps, ml_seqs, cds_rows,
                           proph_seqs)
    log(f"  {len(baits)} real baits; adding controls")

    for hb in make_host_control_baits(inp, selected, member_reps,
                                      proph_rows, host_intervals):
        compute_masks(hb)
        baits.append(hb)
    baits.sort(key=lambda b: b.priority)

    log("[mask] 31-mer statistics + host scan")
    for b in baits:
        compute_masks(b)
    host_like_union = host_scan(inp, baits, panel, host_intervals,
                                os.path.join(work, "host_scan.log"))
    accepted, rej_ledger = apply_gates(baits)
    log(f"[gates] {len(accepted)} accepted, "
        f"{len(rej_ledger) - len(accepted)} rejected")

    for sb in make_shuffled_controls(accepted, host_like_union):
        accepted.append(sb)
        rej_ledger.append(_ledger_row(sb, None))

    log("[report] outputs")
    mrows = manifest_rows(accepted, selected, proph_rows)
    records = [(b.bait_id, b.seq) for b in accepted]
    bl.write_fasta(records, os.path.join(args.out_dir, "baits.fa"))

    hcl = load_host_clades(inp.host_clades_tsv)
    ps_rows = []
    for row in selected:
        clade = row["clade_id"]
        reps = [m for m in member_reps if m["clade_id"] == clade]
        accs = sorted({proph_rows.get(m["member_prophage_id"], {})
                       .get("accession", "") for m in reps} - {""})
        ps_rows.append({**row,
                        "member_reps": ";".join(
                            m["member_prophage_id"] for m in reps),
                        "source_accessions": ";".join(accs),
                        "host_clade_ids_provenance_not_used": ";".join(
                            sorted({hcl[a][0] for a in accs if a in hcl})),
                        "species_provenance_not_used": ";".join(
                            sorted({hcl[a][1] for a in accs if a in hcl}))})
    write_tsv(os.path.join(args.out_dir, "pilot_selection.tsv"), ps_rows,
              ["selection_rank", "clade_id", "genome_id", "tier", "n_members",
               "length_bp", "release_flags", "medoid_member", "member_reps",
               "source_accessions", "host_clade_ids_provenance_not_used",
               "species_provenance_not_used", "selection_basis"])
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
    script = os.path.abspath(__file__)
    checksums = inp.checksums({
        "functional_qc": inp.functional_qc,
        "release_manifest": inp.release_manifest,
        "ml_fasta": inp.ml_fasta,
        "distances": inp.distances,
        "members_json": inp.members_json,
        "tight_clades": inp.tight_clades,
        "prophages_fa": inp.prophages_fa,
        "master_manifest": inp.master_manifest,
        "pharokka_cds": inp.pharokka_cds,
        "host_clades": inp.host_clades_tsv,
    })
    for row in selected:
        for f in ("ml.traversal.json", "partitions.bed"):
            p = os.path.join(inp.clade_dir(row["clade_id"]), f)
            checksums[f"clades/{row['clade_id']}/{f}"] = sha256_file(p)
    host_panel_prov = [{
        "accession": str(e["accession"]),
        "fasta": inp.host_fasta(str(e["accession"])),
        "fasta_sha256": sha256_file(inp.host_fasta(str(e["accession"]))),
        "prophage_intervals_masked": sum(
            len(v) for v in host_intervals.get(str(e["accession"]), {}).values()),
        "clade_id": e["clade_id"],
    } for e in panel]
    params = {
        "K": bl.K, "SELECT_N_A": SELECT_N_A, "SELECT_N_B": SELECT_N_B,
        "MIN_GENOME_BP": MIN_GENOME_BP, "MEMBER_MIN_BP": MEMBER_MIN_BP,
        "MODULE_W": MODULE_W, "JUNCTION_W": JUNCTION_W, "MEMBER_W": MEMBER_W,
        "HOSTCTL_W": HOSTCTL_W, "MODULE_CAP": MODULE_CAP,
        "JUNCTION_CAP": JUNCTION_CAP,
        "N_SHUFFLED_CONTROLS": N_SHUFFLED_CONTROLS,
        "N_HOST_CONTROLS": N_HOST_CONTROLS, "SHUFFLED_SEED": SHUFFLED_SEED,
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
         "script_sha256": sha256_file(script)},
        sort_keys=True).encode()).hexdigest()[:16]
    provenance = {
        "run_id": run_id,
        "design": "ntm/v2/bait/DESIGN.md",
        "task": "design-ntm-prophage",
        "inputs_sha256": checksums,
        "host_panel": host_panel_prov,
        "params": params,
        "script_sha256": sha256_file(script),
        "counts": {
            "selected_genomes": len(selected),
            "member_reps": len(member_reps),
            "baits_accepted": len(accepted),
            "baits_rejected": len(rej_ledger) - len(accepted),
            "eligibility_ledger_rows": len(ledger),
            "host_panel_genomes": len(panel),
        },
        "environment": {"python": sys.version.split()[0],
                        "numpy": np.__version__},
        "reserved_control_class": {
            "control_public_reference": (
                "deferred to task curate-public-mycobacteriophage; the panel "
                "schema accepts these rows without redesign")},
    }
    with open(os.path.join(args.out_dir, "provenance.json"), "w") as fh:
        json.dump(provenance, fh, indent=2, sort_keys=True)
        fh.write("\n")

    tier_counts: Dict[str, int] = {}
    for r in selected:
        tier_counts[r["tier"]] = tier_counts.get(r["tier"], 0) + 1
    cls_counts: Dict[str, int] = {}
    for r in mrows:
        cls_counts[r["bait_class"]] = cls_counts.get(r["bait_class"], 0) + 1
    n_recon = sum(1 for r in mrows if r["bait_class"] == "junction"
                  and r["adjacency_observed"] == "no")
    with open(os.path.join(args.out_dir, "README.md"), "w") as fh:
        fh.write(f"""# NTM prophage bait & junction panel (pilot)

Generated by `../scripts/design_ntm_baits.py` per the preregistered
`../DESIGN.md`. Deterministic run id: `{run_id}`.

- Pilot: {len(selected)} reconstructed ML genomes (tiers: {tier_counts});
  host clade/species labels are recorded as provenance only and never used
  as selection features.
- Baits accepted: {len(mrows)} by class: {cls_counts}.
- Junction baits carry both contributing partitions, the exact breakpoint,
  co-occurrence and observed-adjacency counts; {n_recon} accepted junction
  baits span joins with no observed member adjacency (purely reconstructed).
- Negative controls: {cls_counts.get('control_shuffled_negative', 0)} Markov
  composition-matched negatives (expected: no 31-mer hits) and
  {cls_counts.get('control_host_negative', 0)} host-derived controls
  (expected host-like; validates the host mask).
- `control_public_reference` rows are reserved for the curated public
  mycobacteriophage panel (sibling task `curate-public-mycobacteriophage`).

Files: `baits.fa`, `pilot_selection.tsv`, `bait_manifest.tsv`,
`rejection_ledger.tsv`, `filter_stats.tsv`, `provenance.json`.
Reruns are byte-identical (`--verify-determinism`).
""")
    with open(os.path.join(work, "host_panel.json"), "w") as fh:
        json.dump(host_panel_prov, fh, indent=2, sort_keys=True)

    self_audit(records, mrows, ml_seqs, proph_seqs)
    log(f"[done] panel in {args.out_dir} (run_id {run_id})")
    return 0


def verify_determinism(args: argparse.Namespace) -> int:
    """Run the whole pipeline twice into temp dirs; compare committed bytes."""
    hashes: List[Dict[str, str]] = []
    for _ in range(2):
        with tempfile.TemporaryDirectory() as td:
            out = os.path.join(td, "panel")
            rc = main(["--v2-root", args.v2_root, "--external-root",
                       os.path.join(td, "work"), "--out-dir", out])
            assert rc == 0
            hashes.append({f: sha256_file(os.path.join(out, f))
                           for f in sorted(os.listdir(out))})
    if hashes[0] != hashes[1]:
        for k in sorted(set(hashes[0]) | set(hashes[1])):
            if hashes[0].get(k) != hashes[1].get(k):
                log(f"DETERMINISM FAIL: {k}")
        return 1
    log("DETERMINISM OK: both runs byte-identical")
    return 0


if __name__ == "__main__":
    sys.exit(main())
