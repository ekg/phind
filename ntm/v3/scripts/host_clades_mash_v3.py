#!/usr/bin/env python3
"""NTM v3 — host clades on the all-inclusive cohort (MASH clades, nohup/resume-protected).

Cohort: every canonical genome object from `acquire-ntm-v3`
(`ntm/v3/genomes/canonical_objects/<acc>/<acc>.pansn.fa.gz`, 16,243 objects as
of 2026-09-12). The 18,055 blocked run assemblies have no FASTA on this host
(collaborator upload pending — see `ntm/v3/download_report.md`) and cannot be
sketched; when they land, re-running this driver over the extended filelist
produces fresh clade labels unchanged.

Method parity with v2 (`ntm/v2/scripts/host_clades_mash_v2.py`):
  * mash sketch params `-k 21 -s 10000` (identical to v1/v2)
  * clustering: union-find connected components (single-linkage) at
    `dist <= 0.05` (~95% ANI, species level), components sorted size-desc,
    fresh clade ids `host_clade_0001..`

Differences from v2 mechanics (deliberate, per task ntm-v3-host):
  * **Resume-protected individual sketching.** v2's single
    `mash sketch -l filelist -o host` is one all-or-nothing job; if the
    process (or the LLM stream that owns it) dies, all compute is lost.
    v3 sketches every genome individually into `sketches/<acc>.msh` and
    `mash paste`s them into `host.msh`. Re-invoking the driver skips every
    sketch already on disk, so an interrupted run resumes instead of
    restarting. Disk artifacts are the source of truth (state.json is
    bookkeeping only).
  * **`mash triangle` instead of `mash dist msh msh`.** v2 emitted the full
    N^2 `mash dist` matrix (10.08 GB); v3 emits the complete lower triangle
    (`n*(n-1)/2` pairs) via `mash triangle` — same pairs, half the bytes,
    explicit pair-count completeness check.
  * **nohup-protected heavy phase.** Lesson from operations.jsonl 2026-08-16:
    a provider stream timeout killed the v2 host attempt AFTER compute
    completed. Run the heavy phase (sketch+paste+triangle) detached:
        nohup setsid python3 host_clades_mash_v3.py --phase heavy \
            > heavy.log 2>&1 &
    The agent session may die; compute survives and any later attempt only
    parses/clusters/emits from artifacts already on disk (`--phase light`).

Phases (each idempotent, resumable from disk state):
  filelist  regenerate filelist.txt (sorted canonical genome paths)
  sketch    individual `mash sketch -k 21 -s 10000` per genome (parallel,
            resume: skips non-empty sketches/<acc>.msh)
  paste     `mash paste` sketches -> host.msh (atomic rename; skipped when
            host.msh exists and no sketch is newer)
  triangle  `mash triangle -k 21 -s 10000 -p T host.msh` -> host.dist
            (atomic rename; existence implies complete)
  cluster   parse host.dist (strict completeness check), cluster at <= 0.05,
            annotate species (v3 manifest first, v2 QC list fallback),
            write host_clades.tsv + host_clade_summary.tsv + cluster.log
  spotcheck 20 sampled pairs re-measured with `mash dist` vs triangle values
  receipts  sha256 of all artifacts -> SHA256SUMS
  heavy     filelist + sketch + paste + triangle   (run under nohup)
  light     cluster + spotcheck + receipts          (parse/emit only)
  all       heavy + light

Annotations: species from the v3 acquisition manifest (BV-BRC/NCBI export
label) for objects whose manifest row carries one; v2 QC-passed accession
list (`ntm/v2/inputs/NTM_QC_passed_accession_list.tsv`) for V2_NUMERIC-only
objects (8,212 rows have empty manifest species — measured 0 conflicts).
genus = first token of species; organism = species (v2 schema).

Usage:
    python3 host_clades_mash_v3.py --phase all
    python3 host_clades_mash_v3.py --phase light --base <BASE>
(mash must be on PATH.)
"""
import argparse
import csv
import datetime as dt
import json
import os
import random
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

THRESHOLD = 0.05          # mash distance ~ 95% ANI (v1/v2 convention)
K = 21
S = 10000
EXT = ".pansn.fa.gz"

DEFAULT_BASE = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3"
DEFAULT_MANIFEST = os.path.join(DEFAULT_BASE, "genomes", "v3_acquisition_manifest.tsv")
DEFAULT_V2_QC = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/inputs/NTM_QC_passed_accession_list.tsv"


def utcnow():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log(msg):
    print(f"[{utcnow()}] {msg}", flush=True)


def atomic_write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


class State:
    """Phase bookkeeping in state.json (disk artifacts stay authoritative)."""

    def __init__(self, path):
        self.path = path
        self.data = {}
        if os.path.exists(path):
            try:
                with open(path) as f:
                    self.data = json.load(f)
            except Exception as e:  # corrupt/stale state never blocks resume
                log(f"WARN state.json unreadable ({e}); continuing with empty state")

    def record(self, phase, **kw):
        rec = self.data.get(phase, {})
        rec.update(kw, updated=utcnow())
        self.data[phase] = rec
        atomic_write_json(self.path, self.data)

    def record_start(self, phase):
        self.record(phase, status="running", started=utcnow())
        self.data.setdefault("_order", []).append(phase)


def run(cmd, **kw):
    """subprocess.run with logged command; returns CompletedProcess."""
    log("$ " + " ".join(cmd) if len(cmd) < 20 else "$ " + cmd[0] + " ... (" + str(len(cmd)) + " args)")
    return subprocess.run(cmd, **kw)


# ---------------------------------------------------------------- phases

def phase_filelist(cfg, state):
    outdir = cfg.outdir
    os.makedirs(outdir, exist_ok=True)
    canon = os.path.join(cfg.base, "genomes", "canonical_objects")
    paths = []
    for dp, _, fns in os.walk(canon):
        for fn in fns:
            if fn.endswith(EXT):
                paths.append(os.path.join(dp, fn))
    paths.sort()
    filelist = os.path.join(outdir, "filelist.txt")
    old = None
    if os.path.exists(filelist):
        with open(filelist) as f:
            old = f.read().splitlines()
    if old == paths:
        log(f"filelist: {len(paths)} genomes (unchanged; resume)")
        state.record("filelist", status="ok", genomes=len(paths), changed=False)
        return paths
    with open(filelist + ".tmp", "w") as f:
        f.write("\n".join(paths) + "\n")
    os.replace(filelist + ".tmp", filelist)
    log(f"filelist: wrote {len(paths)} genomes -> {filelist}")
    state.record("filelist", status="ok", genomes=len(paths), changed=True)
    return paths


def _sketch_one(cfg, acc, path):
    """Sketch one genome to sketches/<acc>.msh; returns (acc, rc, seconds, skipped)."""
    target = os.path.join(cfg.outdir, "sketches", acc + ".msh")
    if os.path.exists(target) and os.path.getsize(target) > 0:
        return acc, 0, 0.0, True
    t0 = time.time()
    cmd = ["mash", "sketch", "-k", str(K), "-s", str(S), "-p", "1", "-o", target[:-4], path]
    for attempt in (1, 2):  # one retry per genome
        p = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if p.returncode == 0 and os.path.exists(target) and os.path.getsize(target) > 0:
            return acc, 0, time.time() - t0, False
        if attempt == 1:
            time.sleep(2)
    sys.stderr.write(f"sketch FAILED {acc}: {p.stderr.strip()[:300]}\n")
    return acc, p.returncode, time.time() - t0, False


def phase_sketch(cfg, state, paths):
    sdir = os.path.join(cfg.outdir, "sketches")
    os.makedirs(sdir, exist_ok=True)
    todo = []
    for p in paths:
        acc = os.path.basename(p)[: -len(EXT)]
        target = os.path.join(sdir, acc + ".msh")
        if os.path.exists(target) and os.path.getsize(target) > 0:
            continue
        todo.append((acc, p))
    log(f"sketch: {len(paths)} genomes total, {len(paths) - len(todo)} already on disk (resume), {len(todo)} to sketch")
    state.record("sketch", status="running", total=len(paths), resumed=len(paths) - len(todo), todo=len(todo))
    if not todo:
        state.record("sketch", status="ok", total=len(paths), resumed=len(paths), failed=[])
        return True
    failures = []
    done = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=cfg.workers) as ex:
        futs = {ex.submit(_sketch_one, cfg, acc, p): acc for acc, p in todo}
        for fut in as_completed(futs):
            acc, rc, secs, skipped = fut.result()
            done += 1
            if rc != 0:
                failures.append(acc)
            if done % 250 == 0 or done == len(todo):
                rate = done / max(time.time() - t0, 1e-9)
                log(f"sketch: {done}/{len(todo)} done ({rate:.1f}/s), failures so far: {len(failures)}")
    if failures:
        log(f"sketch: FAILED for {len(failures)} genomes: {sorted(failures)[:10]}{' ...' if len(failures) > 10 else ''}")
        state.record("sketch", status="failed", failed=sorted(failures))
        return False
    state.record("sketch", status="ok", total=len(paths), resumed=len(paths) - len(todo), failed=[])
    log(f"sketch: all {len(paths)} sketches present on disk")
    return True


def phase_paste(cfg, state, paths):
    host = os.path.join(cfg.outdir, "host.msh")
    sdir = os.path.join(cfg.outdir, "sketches")
    sketches = [os.path.join(sdir, os.path.basename(p)[: -len(EXT)] + ".msh") for p in paths]
    missing = [s for s in sketches if not (os.path.exists(s) and os.path.getsize(s) > 0)]
    if missing:
        log(f"paste: ABORT — {len(missing)} sketches missing (run --phase sketch first): {missing[:3]}")
        state.record("paste", status="failed", missing=len(missing))
        return False
    if os.path.exists(host) and os.path.getsize(host) > 0:
        host_mtime = os.path.getmtime(host)
        stale = [s for s in sketches if os.path.getmtime(s) > host_mtime]
        if not stale:
            log(f"paste: host.msh exists and no sketch is newer (resume; {len(sketches)} sketches)")
            state.record("paste", status="ok", sketches=len(sketches), resumed=True)
            return True
        log(f"paste: host.msh stale ({len(stale)} newer sketches) — re-pasting")
    tmp = os.path.join(cfg.outdir, "host.tmp.msh")  # mash appends .msh only if absent
    for junk in (tmp, tmp + ".msh"):
        if os.path.exists(junk):
            os.remove(junk)
    p = run(["mash", "paste", "host.tmp.msh"] + [os.path.relpath(s, cfg.outdir) for s in sketches], cwd=cfg.outdir)
    written = tmp if os.path.exists(tmp) else (tmp + ".msh" if os.path.exists(tmp + ".msh") else None)
    if p.returncode != 0 or written is None:
        log("paste: mash paste failed")
        state.record("paste", status="failed")
        return False
    os.replace(written, host)
    n = count_sketches(host)
    log(f"paste: host.msh = {len(sketches)} sketches ({n} counted by mash info), {os.path.getsize(host) / 1e6:.0f} MB")
    state.record("paste", status="ok", sketches=len(sketches), mash_info_count=n)
    return n == len(sketches)


def count_sketches(msh):
    """Sketch count via `mash info -t` (one row per sketch after a # header)."""
    p = subprocess.run(["mash", "info", "-t", msh], capture_output=True, text=True)
    if p.returncode != 0:
        return -1
    return sum(1 for ln in p.stdout.splitlines() if ln and not ln.startswith("#"))


def phase_triangle(cfg, state, paths):
    host = os.path.join(cfg.outdir, "host.msh")
    dist = os.path.join(cfg.outdir, "host.dist")
    n = len(paths)
    expected = n * (n - 1) // 2
    if os.path.exists(dist) and os.path.getsize(dist) > 0 and os.path.getmtime(dist) >= os.path.getmtime(host):
        # host.dist is only ever published via atomic rename after mash exit 0
        log(f"triangle: host.dist exists and >= host.msh mtime (resume; completeness re-verified in cluster phase)")
        state.record("triangle", status="ok", resumed=True, expected_pairs=expected)
        return True
    tmp = dist + ".tmp"
    if os.path.exists(tmp):
        os.remove(tmp)
    log(f"triangle: mash triangle -k {K} -s {S} -p {cfg.threads} host.msh -> host.dist ({n} genomes, {expected} pairs expected)")
    t0 = time.time()
    with open(tmp, "w") as out:
        p = subprocess.run(["mash", "triangle", "-k", str(K), "-s", str(S), "-p", str(cfg.threads), host],
                           stdout=out, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        log(f"triangle: mash triangle FAILED: {p.stderr.strip()[:300]}")
        os.remove(tmp)
        state.record("triangle", status="failed")
        return False
    os.replace(tmp, dist)
    secs = time.time() - t0
    log(f"triangle: wrote host.dist ({os.path.getsize(dist) / 1e9:.2f} GB) in {secs / 60:.1f} min; completeness checked in cluster phase")
    state.record("triangle", status="ok", seconds=round(secs, 1), expected_pairs=expected,
                 bytes=os.path.getsize(dist))
    return True


# ------------------------------------------------- clustering (v2 logic)

def load_species_annotations(manifest_path, v2qc_path):
    """acc -> species. Manifest (BV-BRC/NCBI export label) first; v2 QC list
    fallback for V2_NUMERIC-only objects whose manifest species is empty."""
    ann = {}
    with open(manifest_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            ca = (row.get("canonical_acc") or "").strip()
            sp = (row.get("species") or "").strip()
            if ca and sp and ca not in ann:
                ann[ca] = sp
    n_manifest = len(ann)
    fallback = 0
    with open(v2qc_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            acc = (row.get("accession") or "").strip()
            sp = (row.get("species") or "").strip()
            if acc and sp and acc not in ann:
                ann[acc] = sp
                fallback += 1
    return ann, n_manifest, fallback


def parse_triangle_iter(dist_path, n_expected=None):
    """Yield (index, line) data rows from a `mash triangle` output file.

    File format (verified empirically, mash 2.3): first line "\t<n genomes>",
    then n rows — row i is name_i followed by i tab-separated distances to
    genomes 0..i-1 in order. ("Max p-value:" is written to stderr, not here.)
    """
    with open(dist_path) as f:
        line = f.readline().rstrip("\n")
        if line.startswith("#") or line.startswith("Max p-value"):
            line = f.readline().rstrip("\n")
        if line.startswith("\t") and line[1:].isdigit():
            if n_expected is not None and int(line[1:]) != n_expected:
                raise ValueError(f"triangle header genome count {line[1:]} != expected {n_expected}")
            line = f.readline().rstrip("\n")
        idx = 0
        while line:
            yield idx, line
            idx += 1
            line = f.readline().rstrip("\n")


def phase_cluster(cfg, state, paths):
    outdir = cfg.outdir
    dist_path = os.path.join(outdir, "host.dist")
    n = len(paths)
    expected_pairs = n * (n - 1) // 2
    filelist_accs = [os.path.basename(p)[: -len(EXT)] for p in paths]

    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    # --- strict parse + union-find at <= THRESHOLD (single streaming pass)
    names = []
    pairs = 0
    unions = 0
    malformed = 0
    t0 = time.time()
    spot = pick_spot_pairs(cfg, n)  # {row_idx: [col idx,...]} for spotcheck
    spot_values = {}
    for i, line in parse_triangle_iter(dist_path, n_expected=n):
        toks = line.split("\t")
        name = toks[0]
        dists = toks[1:]
        if len(toks) != i + 1:
            malformed += 1
            if malformed <= 5:
                log(f"cluster: MALFORMED row {i}: {len(toks)} fields, expected {i + 1}")
            continue
        acc = os.path.basename(name)[: -len(EXT)] if name.endswith(EXT) else os.path.basename(name)
        names.append(acc)
        parent.setdefault(acc, acc)
        if i in spot:
            keep = {}
        for j, ds in enumerate(dists):
            pairs += 1
            if i in spot:
                keep[j] = ds
            if ds[0] != "0":        # any value not starting '0' is >= 0.1
                continue
            if float(ds) <= THRESHOLD:
                other = names[j]
                parent.setdefault(other, other)
                union(acc, other)
                unions += 1
        if i in spot:
            spot_values[i] = keep
    log(f"cluster: parsed {len(names)} rows, {pairs} pairs "
        f"(expected {expected_pairs}) in {time.time() - t0:.0f}s; unions={unions}, malformed={malformed}")

    # --- completeness gate
    problems = []
    if len(names) != n:
        problems.append(f"row count {len(names)} != {n}")
    if pairs != expected_pairs:
        problems.append(f"pair count {pairs} != {expected_pairs} (n*(n-1)/2)")
    if malformed:
        problems.append(f"{malformed} malformed rows")
    if set(names) != set(filelist_accs):
        missing = set(filelist_accs) - set(names)
        extra = set(names) - set(filelist_accs)
        problems.append(f"name mismatch: missing={len(missing)} extra={len(extra)}")
    if problems:
        for p_ in problems:
            log(f"cluster: FAIL {p_}")
        state.record("cluster", status="failed", problems=problems)
        return False

    # --- connected components -> clades, size-desc
    comps = defaultdict(list)
    for acc in filelist_accs:
        comps[find(acc) if acc in parent else acc].append(acc)
    clades = sorted(comps.values(), key=lambda ms: (-len(ms), ms))
    clade_of = {}
    clade_meta = []
    for i, members in enumerate(clades, 1):
        cid = f"host_clade_{i:04d}"
        clade_meta.append((cid, members))
        for m in members:
            clade_of[m] = cid

    # --- annotations
    ann, n_manifest, n_fallback = load_species_annotations(cfg.manifest, cfg.v2qc)
    missing_ann = [a for a in filelist_accs if a not in ann]

    # --- host_clades.tsv (v2 schema)
    rows = sorted((a, clade_of[a], ann.get(a, ""), ann.get(a, "").split(" ", 1)[0] if ann.get(a) else "", ann.get(a, ""))
                  for a in filelist_accs)
    with open(os.path.join(outdir, "host_clades.tsv"), "w") as f:
        f.write("accession\thost_clade_id\tspecies\tgenus\torganism\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    # --- host_clade_summary.tsv
    clade_counter = {}
    for cid, members in clade_meta:
        clade_counter[cid] = Counter(ann.get(m, "unknown") or "unknown" for m in members)
    with open(os.path.join(outdir, "host_clade_summary.tsv"), "w") as f:
        f.write("host_clade_id\tcount\tdominant_species\tspecies_distribution\n")
        for cid, members in clade_meta:
            cc = clade_counter[cid]
            dist_str = ";".join(f"{s_}:{n_}" for s_, n_ in cc.most_common())
            f.write(f"{cid}\t{len(members)}\t{cc.most_common(1)[0][0]}\t{dist_str}\n")

    # --- cluster.log (v2 style summary)
    sizes = Counter(len(m) for _, m in clade_meta)
    singletons = sizes.get(1, 0)
    cohort_species = Counter(ann.get(a, "unknown") or "unknown" for a in filelist_accs)
    tb = sum(v for k, v in cohort_species.items() if "tuberculosis" in k.lower() and "paratuberculosis" not in k.lower())
    absc = sum(v for k, v in cohort_species.items() if "abscessus" in k.lower())
    lines = []
    lines.append(f"filelist genomes: {n}")
    lines.append(f"host.dist pairs: {pairs}; union ops: {unions}; malformed: {malformed}")
    lines.append(f"clades: {len(clade_meta)}; assigned: {len(clade_of)} accessions; singletons: {singletons}")
    lines.append(f"annotation sources: {n_manifest} from v3 manifest, {n_fallback} v2-QC fallback")
    lines.append("")
    lines.append("=== summary ===")
    lines.append("clade size distribution (size:count, top):")
    for size, cnt in sorted(sizes.items(), reverse=True)[:15]:
        lines.append(f"  n={size}: {cnt} clades")
    lines.append(f"  (sizes <= threshold omitted: {sum(c for s, c in sizes.items() if s <= 6)} clades of size 1-6)")
    lines.append("")
    lines.append("=== composition ===")
    lines.append(f"total genomes: {n}")
    lines.append(f"TB-complex ('tuberculosis', excl. paratuberculosis): {tb} ({tb / n:.1%})")
    lines.append(f"abscessus complex ('abscessus'): {absc} ({absc / n:.1%})")
    lines.append("top 15 species:")
    for s_, n_ in cohort_species.most_common(15):
        lines.append(f"  {n_:6d}  {s_}")
    lines.append("")
    lines.append("top clades (dominant species):")
    for cid, members in clade_meta[:8]:
        lines.append(f"  {cid}: n={len(members)} dominant={clade_counter[cid].most_common(1)[0][0]}")
    lines.append("")
    lines.append("=== validation ===")
    ok = True
    # every genome exactly one clade, one row
    if len(rows) != n or len(set(r[0] for r in rows)) != n or set(r[0] for r in rows) != set(filelist_accs):
        ok = False
        lines.append("FAIL row completeness (rows != filelist or dup/missing accessions)")
    else:
        lines.append(f"row completeness: {n} rows == {n} filelist accessions; 0 duplicates")
    if missing_ann:
        ok = False
        lines.append(f"FAIL {len(missing_ann)} accessions missing species annotation: {missing_ann[:5]}")
    else:
        lines.append("annotation join: all accessions resolve to a species (manifest or v2 QC list)")
    # no genome in >1 clade guaranteed by construction; re-assert from rows
    multi = [a for a, c in Counter(r[0] for r in rows).items() if c > 1]
    if multi:
        ok = False
        lines.append(f"FAIL genomes in >1 clade: {multi[:5]}")
    else:
        lines.append("no genome in >1 clade")
    top_cid, top_members = clade_meta[0]
    frac = len(top_members) / n
    dom = clade_counter[top_cid].most_common(1)[0][0].lower()
    lines.append(f"top clade {top_cid} n={len(top_members)} ({frac:.1%}) dominant={clade_counter[top_cid].most_common(1)[0][0]}")
    if frac > 0.40 and not ("abscessus" in dom or "tuberculosis" in dom):
        ok = False
        lines.append(f"FAIL top clade {frac:.1%} > 40% and not abscessus/TB-complex-dominated")
    if ok:
        lines.append("validation: PASS")
    else:
        lines.append("validation: FAIL")
    with open(os.path.join(outdir, "cluster.log"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)

    state.record("cluster", status="ok" if ok else "failed", genomes=n, clades=len(clade_meta),
                 pairs=pairs, unions=unions, singletons=singletons,
                 annotation_manifest=n_manifest, annotation_v2qc=n_fallback)
    # stash spot values for the spotcheck phase
    spot_out = {str(i): spot_values[i] for i in spot_values} if spot_values else {}
    atomic_write_json(os.path.join(outdir, "spot_pairs.json"),
                      {"row_names": names, "spot_values": spot_out} if spot_values else {})
    return ok and pairs == expected_pairs


def pick_spot_pairs(cfg, n):
    """Deterministic sample: 20 distinct row indices (>=1), one column each."""
    if n < 3:
        return {}
    r = random.Random(cfg.sample_seed)
    rows = r.sample(range(1, n), min(cfg.sample_pairs, n - 1))
    return {i: True for i in rows}


def phase_spotcheck(cfg, state, paths):
    outdir = cfg.outdir
    spot_path = os.path.join(outdir, "spot_pairs.json")
    if not os.path.exists(spot_path):
        log("spotcheck: no spot_pairs.json (run --phase cluster first)")
        state.record("spotcheck", status="failed")
        return False
    with open(spot_path) as f:
        data = json.load(f)
    names = data["row_names"]
    spot = data["spot_values"]
    r = random.Random(cfg.sample_seed + 1)
    results = []
    ok = True
    for i_str, cols in sorted(spot.items(), key=lambda kv: int(kv[0])):
        i = int(i_str)
        col = r.randrange(i)
        ds = cols[str(col)]
        tri_val = float(ds)
        a, b = names[i], names[col]
        sa = os.path.join(outdir, "sketches", a + ".msh")
        sb = os.path.join(outdir, "sketches", b + ".msh")
        p = subprocess.run(["mash", "dist", "-p", "1", sa, sb],
                           capture_output=True, text=True)
        if p.returncode != 0:
            results.append((a, b, tri_val, None, "mash dist failed"))
            ok = False
            continue
        line = p.stdout.splitlines()[0].split("\t")
        dist_val = float(line[2])
        delta = abs(dist_val - tri_val)
        passed = delta <= 1e-4
        ok = ok and passed
        results.append((a, b, tri_val, dist_val, "OK" if passed else f"MISMATCH delta={delta}"))
    with open(os.path.join(outdir, "spotcheck.tsv"), "w") as f:
        f.write("genome_a\tgenome_b\ttriangle_dist\tmash_dist\tdelta\tresult\n")
        for a, b, t_, d_, res in results:
            delta = "" if d_ is None else f"{abs(d_ - t_):.2e}"
            f.write(f"{a}\t{b}\t{t_}\t{d_ if d_ is not None else 'NA'}\t{delta}\t{res}\n")
    n_ok = sum(1 for *_, res in results if res == "OK")
    log(f"spotcheck: {n_ok}/{len(results)} pairs match direct mash dist (<=1e-4)")
    state.record("spotcheck", status="ok" if ok else "failed", pairs=len(results), matched=n_ok)
    return ok


def phase_receipts(cfg, state, paths):
    import hashlib
    outdir = cfg.outdir
    files = ["filelist.txt", "host.msh", "host.dist", "host_clades.tsv",
             "host_clade_summary.tsv", "cluster.log", "spotcheck.tsv"]
    lines = []
    for fn in files:
        p = os.path.join(outdir, fn)
        if not os.path.exists(p):
            continue
        h = hashlib.sha256()
        with open(p, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 22), b""):
                h.update(chunk)
        lines.append(f"{h.hexdigest()}  {os.path.getsize(p):>12}  {fn}")
    with open(os.path.join(outdir, "SHA256SUMS"), "w") as f:
        f.write("# sha256  size(bytes)  file — NTM v3 host_clades artifacts (task ntm-v3-host)\n")
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)
    state.record("receipts", status="ok", files=len(lines))
    return True


PHASES = ["filelist", "sketch", "paste", "triangle", "cluster", "spotcheck", "receipts"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phase", default="all",
                    choices=["all", "heavy", "light"] + PHASES)
    ap.add_argument("--base", default=DEFAULT_BASE)
    ap.add_argument("--manifest", default=DEFAULT_MANIFEST)
    ap.add_argument("--v2qc", default=DEFAULT_V2_QC)
    ap.add_argument("--workers", type=int, default=32, help="parallel sketch processes")
    ap.add_argument("--threads", type=int, default=64, help="mash triangle threads")
    ap.add_argument("--sample-pairs", type=int, default=20)
    ap.add_argument("--sample-seed", type=int, default=42)
    cfg = ap.parse_args()
    cfg.outdir = os.path.join(cfg.base, "host_clades")
    os.makedirs(cfg.outdir, exist_ok=True)
    state = State(os.path.join(cfg.outdir, "state.json"))
    log(f"host_clades_mash_v3 phase={cfg.phase} base={cfg.base} outdir={cfg.outdir}")

    if cfg.phase in ("heavy", "all"):
        seq = PHASES if cfg.phase == "all" else ["sketch", "paste", "triangle"]
    elif cfg.phase == "light":
        seq = ["cluster", "spotcheck", "receipts"]
    else:
        seq = [cfg.phase]
    # filelist is cheap and idempotent; always run it first so downstream
    # phases have the genome list (and resume detection works)
    if seq[0] != "filelist":
        seq = ["filelist"] + seq

    if not shutil.which("mash"):
        log("FATAL: mash not on PATH")
        sys.exit(3)

    paths = None
    rc = 0
    for phase in seq:
        state.record_start(phase)
        t0 = time.time()
        if phase == "filelist":
            paths = phase_filelist(cfg, state)
            ok = True
        elif phase == "sketch":
            ok = phase_sketch(cfg, state, paths)
        elif phase == "paste":
            ok = phase_paste(cfg, state, paths)
        elif phase == "triangle":
            ok = phase_triangle(cfg, state, paths)
        elif phase == "cluster":
            ok = phase_cluster(cfg, state, paths)
        elif phase == "spotcheck":
            ok = phase_spotcheck(cfg, state, paths)
        elif phase == "receipts":
            ok = phase_receipts(cfg, state, paths)
        else:
            ok = False
        log(f"phase {phase}: {'OK' if ok else 'FAILED'} ({time.time() - t0:.0f}s)")
        rc = 0 if ok else 2
        if not ok:
            log(f"stopping at failed phase: {phase} (rerun this driver to resume from disk state)")
            break
    sys.exit(rc)


if __name__ == "__main__":
    main()
