#!/usr/bin/env python3
"""Stratified Tier E1 + E2 design from the frozen WGS universe (RUN_PLAN §4-§6).

Inputs
  --frozen        metadata/ecoli_wgs_frozen.tsv        (frozen WGS universe)
  --titles        metadata/bioproject_titles.tsv       (projects >= 100 runs)
Outputs
  metadata/bioproject_table.tsv       per-BioProject WGS run counts
  manifests/tierE1_manifest.json + tierE1_frozen.tsv
  manifests/tierE2_manifest.json + tierE2_frozen.tsv

Frozen rules (RUN_PLAN §5/§6; seed 20260822 everywhere):
  E1: BioProjects with >= 100 frozen WGS runs matched to archetypes by
      case-insensitive title substring:
        ST131: 'st131'
        global-diversity: 'global'|'worldwide'|'international'|'diversity'
        longitudinal: 'longitudinal'|'timecourse'|'time course'|'time-series'|
                      'persistence'|'follow-up'
      Greedy whole-project take in archetype order, cumulative union <= 12000;
      accept 8000-12000 (or take-all + document if less).
  E2: frame = frozen WGS runs minus all runs of E1 BioProjects;
      stratify by BioProject, cap C=60; per-project RNG seeded
      f"{seed}|{bioproject}"; project order = random.Random(seed).shuffle of
      projects sorted by (run count desc, accession asc); accumulate to
      exactly 12000.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SEED = 20260822
E1_MIN, E1_MAX = 8000, 12000
E1_HARD_CAP = 12000
E2_TARGET = 12000
E2_PER_PROJECT_CAP = 60
MIN_PROJECT_FOR_TITLE = 100

ARCHETYPES = [
    ("ST131", [re.compile(p, re.I) for p in (r"st131",)]),
    ("global-diversity", [re.compile(p, re.I) for p in
                          (r"global", r"worldwide", r"international", r"diversity")]),
    ("longitudinal", [re.compile(p, re.I) for p in
                      (r"longitudinal", r"timecourse", r"time course", r"time-series",
                       r"persistence", r"follow-up")]),
]


def http_get(url, retries=3, backoff=(5, 20, 60)):
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "local-ecoli-logan/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            if attempt < retries:
                time.sleep(backoff[min(attempt, len(backoff) - 1)])
    raise RuntimeError(f"GET failed after {retries} retries: {url}: {last}")


def load_frozen(path):
    by_proj = {}
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    for r in rows:
        by_proj.setdefault(r["bioproject"] or "(none)", []).append(r["run"])
    for p in by_proj:
        by_proj[p] = sorted(set(by_proj[p]))
    return rows, by_proj


def fetch_titles(accs, log):
    """db=bioproject esummary titles, <=1 req/s, sha256-logged."""
    out = {}
    todo = [a for a in accs if a != "(none)"]
    for i in range(0, len(todo), 200):
        batch = todo[i:i + 200]
        q = urllib.parse.urlencode({"db": "bioproject", "id": ",".join(batch),
                                    "retmode": "json"})
        body = http_get(f"{EUTILS}/esummary.fcgi?{q}")
        try:
            data = json.loads(body).get("result", {})
        except Exception:  # noqa: BLE001
            data = {}
        for a in batch:
            e = data.get(a)
            if e:
                out[a] = e.get("title", "")
        log.write(json.dumps({"event": "bioproject_esummary", "n": len(batch),
                              "n_titles": len(out),
                              "sha256": hashlib.sha256(body).hexdigest(),
                              "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}) + "\n")
        time.sleep(1.1)
    return out


def match_archetype(title):
    for name, pats in ARCHETYPES:
        for p in pats:
            if p.search(title or ""):
                return name
    return None


def write_frozen_tsv(path, runs, row_by_run, tier, extra_proj=None):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "scientific_name", "queries", "library_strategy",
                    "platform", "bioproject", "uids", "tier_e1_archetype"])
        for run in runs:
            r = row_by_run[run]
            w.writerow([run, r["scientific_name"], tier, r["library_strategy"],
                        r["platform"], r["bioproject"], r["uids"],
                        (extra_proj or {}).get(r["bioproject"], "")])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frozen", required=True)
    ap.add_argument("--out-dir", required=True, help="metadata dir (NVMe)")
    ap.add_argument("--manifest-dir", required=True, help="repo manifests dir")
    ap.add_argument("--log", required=True)
    ap.add_argument("--skip-titles", action="store_true",
                    help="reuse existing metadata/bioproject_titles.tsv")
    args = ap.parse_args()

    rows, by_proj = load_frozen(args.frozen)
    row_by_run = {r["run"]: r for r in rows}
    n_runs = len(row_by_run)
    print(f"frozen WGS universe: {n_runs} runs, {len(by_proj)} BioProjects")

    # per-project table (frozen order: count desc, accession asc)
    proj_counts = sorted(((p, len(rs)) for p, rs in by_proj.items()),
                         key=lambda kv: (-kv[1], kv[0]))
    with open(os.path.join(args.out_dir, "bioproject_table.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["bioproject", "n_wgs_runs"])
        for p, n in proj_counts:
            w.writerow([p, n])

    titles_path = os.path.join(args.out_dir, "bioproject_titles.tsv")
    big = [p for p, n in proj_counts if n >= MIN_PROJECT_FOR_TITLE]
    if args.skip_titles and os.path.exists(titles_path):
        titles = {r["bioproject"]: r["title"] for r in
                  csv.DictReader(open(titles_path), delimiter="\t")}
    else:
        with open(args.log, "a") as log:
            titles = fetch_titles([p for p, _ in proj_counts if p != "(none)"
                                   and dict(proj_counts)[p] >= MIN_PROJECT_FOR_TITLE], log)
        with open(titles_path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["bioproject", "n_wgs_runs", "title"])
            for p, n in proj_counts:
                if p in titles:
                    w.writerow([p, n, titles[p]])
    print(f"titles fetched for {len(titles)} projects >= {MIN_PROJECT_FOR_TITLE} runs")

    # ---- Tier E1 (RUN_PLAN §5) ----
    e1_projects, e1_taken, e1_runs, cumulative = [], {a: [] for a, _ in ARCHETYPES}, [], 0
    e1_candidates = []  # (archetype, project, n, title)
    for arch, _ in ARCHETYPES:
        for p, n in proj_counts:
            if p in titles and match_archetype(titles[p]) == arch:
                e1_candidates.append((arch, p, n, titles[p]))
    for arch, p, n, title in e1_candidates:
        if cumulative + n <= E1_HARD_CAP:
            e1_projects.append(p)
            e1_taken[arch].append(p)
            e1_runs.extend(by_proj[p])
            cumulative += n
    e1_runs = sorted(set(e1_runs))
    if not (E1_MIN <= len(e1_runs) <= E1_MAX):
        print(f"NOTE: E1 union = {len(e1_runs)} runs (outside {E1_MIN}-{E1_MAX}); "
              f"taking whole archetype union and documenting (RUN_PLAN §5)")

    e1_proj_arch = {p: arch for arch in e1_taken for p in e1_taken[arch]}
    e1_frozen = os.path.join(args.manifest_dir, "tierE1_frozen.tsv")
    write_frozen_tsv(e1_frozen, e1_runs, row_by_run, "E1", e1_proj_arch)
    e1_sha = hashlib.sha256(open(e1_frozen, "rb").read()).hexdigest()
    e1_manifest = {
        "tier": "E1", "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "rule": "RUN_PLAN §5: BioProjects >= 100 frozen WGS runs, title-archetype match, "
                "greedy whole-project take in archetype order, cumulative <= 12000",
        "archetypes": {a: {"match_strings": [p.pattern for p in pats]}
                       for a, pats in ARCHETYPES},
        "projects": [{"bioproject": p, "archetype": e1_proj_arch[p],
                      "n_wgs_runs": len(by_proj[p]), "title": titles.get(p, ""),
                      "justification": f"{e1_proj_arch[p]} archetype by frozen title match"}
                     for p in e1_projects],
        "n_runs": len(e1_runs), "frozen_tsv_sha256": e1_sha, "seed": SEED,
    }
    with open(os.path.join(args.manifest_dir, "tierE1_manifest.json"), "w") as fh:
        json.dump(e1_manifest, fh, indent=1, sort_keys=True)
    print(f"Tier E1: {len(e1_projects)} projects, {len(e1_runs)} runs, sha256={e1_sha[:16]}…")
    for arch, _ in ARCHETYPES:
        print(f"  {arch}: {len(e1_taken[arch])} projects "
              f"({', '.join(e1_taken[arch][:6])}{'…' if len(e1_taken[arch]) > 6 else ''})")

    # ---- Tier E2 (RUN_PLAN §6) ----
    e1_proj_set = set(e1_projects)
    frame = {p: rs for p, rs in by_proj.items() if p not in e1_proj_set}
    project_order = [p for p, _ in sorted(frame.items(), key=lambda kv: (-len(kv[1]), kv[0]))]
    random.Random(SEED).shuffle(project_order)
    e2_runs, strata = [], {}
    for p in project_order:
        if len(e2_runs) >= E2_TARGET:
            break
        avail = frame[p]
        rng = random.Random(f"{SEED}|{p}")
        take = rng.sample(avail, min(len(avail), E2_PER_PROJECT_CAP))
        quota = E2_TARGET - len(e2_runs)
        take = sorted(take[:quota])
        e2_runs.extend(take)
        strata[p] = {"n_avail": len(avail), "cap": E2_PER_PROJECT_CAP, "n_taken": len(take)}
    e2_runs = sorted(set(e2_runs))
    # keep exact target if frame supports it
    if len(e2_runs) > E2_TARGET:
        e2_runs = e2_runs[:E2_TARGET]

    e2_frozen = os.path.join(args.manifest_dir, "tierE2_frozen.tsv")
    write_frozen_tsv(e2_frozen, e2_runs, row_by_run, "E2")
    e2_sha = hashlib.sha256(open(e2_frozen, "rb").read()).hexdigest()
    e2_manifest = {
        "tier": "E2", "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "rule": "RUN_PLAN §6: stratified random sample by BioProject, cap 60/project, "
                "frame = frozen WGS universe minus E1 BioProjects, target 12000",
        "seed": SEED, "per_project_rng": "random.Random(f'{seed}|{bioproject}')",
        "project_order_rule": "random.Random(seed).shuffle of projects sorted by "
                              "(n_wgs_runs desc, accession asc)",
        "n_projects": len(strata), "n_runs": len(e2_runs),
        "frozen_tsv_sha256": e2_sha,
        "exclusion_rules": ["library_strategy != WGS (universe filter)",
                            "scientific_name not startswith 'Escherichia coli'",
                            "runs of Tier E1 BioProjects excluded entirely",
                            "duplicate run accessions"],
        "strata_sizes": {p: s["n_taken"] for p, s in sorted(strata.items())},
        "strata_meta": strata,
    }
    with open(os.path.join(args.manifest_dir, "tierE2_manifest.json"), "w") as fh:
        json.dump(e2_manifest, fh, indent=1, sort_keys=True)
    print(f"Tier E2: {len(strata)} projects, {len(e2_runs)} runs, sha256={e2_sha[:16]}…")
    overlaps = set(e1_runs) & set(e2_runs)
    print(f"E1∩E2 run overlap: {len(overlaps)} (must be 0)")
    return 0 if not overlaps else 1


if __name__ == "__main__":
    sys.exit(main())
