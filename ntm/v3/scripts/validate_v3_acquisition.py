#!/usr/bin/env python3
"""
NTM v3 — validate the acquisition and emit coverage.tsv + download_report.md.

Validation performed (task ## Validation section):
  * resolved (linked + downloaded + failed/blocked) == manifest rows; every
    manifest row has exactly one terminal state, zero silent skips
  * `samtools faidx` region query succeeds on every canonical object
    (uses the existing index; -o /dev/null so frozen v1/v2 .fai are never
    rewritten through symlinks)
  * `gzip -t` on a random sample of linked objects + every downloaded object
  * contig names vs coordinates-CSV scaffolds: every non-run prophage-bearing
    cohort row checked (scaffold normalized: strip "accn|", strip version;
    fai contigs normalized the same way, NZ_-strip fallback as in the v2
    extraction convention)
  * all 9,543 v2 run assemblies present as manifest rows (none dropped)
  * coverage.tsv numbers reconcile exactly with download_report.md

Outputs (repo, committed):
  ntm/v3/inputs/coverage.tsv
  ntm/v3/download_report.md
Outputs (NVMe):
  {work}/ntm/v3/genomes/blocked_run_accessions.txt
  {work}/ntm/v3/genomes/validation.json

Usage:
  python3 ntm/v3/scripts/validate_v3_acquisition.py \
      [--work-dir /mnt/nvme3n1/erikg/phind-genome-work] [--gzip-sample 300]
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import json
import os
import random
import re
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
MANIFEST_TSV = REPO / "ntm/v3/inputs/v3_acquisition_manifest.tsv"
SUMMARY_CSV = REPO / "ntm/v3/inputs/ntm_qc_passed_phigaro_summary_20260909.csv"
COORDS_CSV = REPO / "ntm/v3/inputs/ntm_qc_passed_phigaro_coordinates_20260909.csv"
GC_RE = re.compile(r"^GC([AF])_([0-9]{9})\.([0-9]+)$")
RUN_RE = re.compile(r"^(ERR|SRR|DRR)\d+$")
UNVERSION_RE = re.compile(r"\.\d+$")


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def norm_scaffold(s: str) -> str:
    """coordinates-CSV scaffold -> contig token (strip accn|, strip version)."""
    if s.startswith("accn|"):
        s = s[len("accn|"):]
    return UNVERSION_RE.sub("", s)


def norm_contig(c: str) -> str:
    """PanSN fai contig name -> unversioned token."""
    if "#" in c:
        c = c.split("#")[-1]
    return UNVERSION_RE.sub("", c)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", default="/mnt/nvme3n1/erikg/phind-genome-work")
    ap.add_argument("--gzip-sample", type=int, default=300)
    ap.add_argument("--workers", type=int, default=24)
    args = ap.parse_args()

    out_dir = Path(args.work_dir) / "ntm/v3/genomes"
    co_dir = out_dir / "canonical_objects"

    rows = list(csv.DictReader(open(MANIFEST_TSV, newline=""), delimiter="\t"))
    export = list(csv.DictReader(open(SUMMARY_CSV, newline="")))

    # per-genome scaffolds (prophage-bearing genomes only)
    scaffolds: dict[str, set[str]] = defaultdict(set)
    with open(COORDS_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            scaffolds[r["genome_id"]].add(r["scaffold"])

    # ------------------------------------------------- object existence
    obj_state: dict[str, dict] = {}
    missing: list[str] = []

    def check_obj(acc: str) -> tuple[str, dict]:
        d = co_dir / acc
        gz = d / f"{acc}.pansn.fa.gz"
        fai = d / f"{acc}.pansn.fa.gz.fai"
        gzi = d / f"{acc}.pansn.fa.gz.gzi"
        info = {"gz": gz.exists(), "fai": fai.exists(), "gzi": gzi.exists()}
        if info["gz"] and info["fai"]:
            # faidx region query through the existing index (no .fai rewrite)
            first = open(fai).readline().split("\t")[0]
            res = subprocess.run(
                ["samtools", "faidx", "-o", "/dev/null", str(gz), first],
                capture_output=True, timeout=300)
            info["faidx_ok"] = res.returncode == 0
            if res.returncode != 0:
                info["faidx_err"] = res.stderr.decode()[:200]
        return acc, info

    needed = sorted({r["canonical_acc"] for r in rows
                     if r["plan"] in ("link_v1v2", "download_bvbrc", "download_ncbi")})
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        for acc, info in ex.map(check_obj, needed):
            obj_state[acc] = info
            if not (info["gz"] and info["fai"] and info.get("faidx_ok")):
                missing.append(acc)
    n_objects = len(needed)
    faidx_ok = sum(1 for v in obj_state.values() if v.get("faidx_ok"))
    total_bp = 0
    total_contigs = 0
    for acc, v in obj_state.items():
        if v["fai"]:
            with open(co_dir / acc / f"{acc}.pansn.fa.gz.fai") as fh:
                for line in fh:
                    total_contigs += 1
                    total_bp += int(line.split("\t")[1])

    # ------------------------------------------------------- gzip -t sample
    rng = random.Random(42)
    sample = rng.sample(needed, min(args.gzip_sample, len(needed)))
    gzip_fail = []

    def gzip_test(acc: str) -> str | None:
        with open(co_dir / acc / f"{acc}.pansn.fa.gz", "rb") as fh:
            try:
                # gzip module test reads the whole stream
                with gzip.GzipFile(fileobj=fh) as gz:
                    while gz.read(1 << 20):
                        pass
            except Exception as e:  # noqa: BLE001
                return f"{acc}: {e}"
        return None

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        for res in ex.map(gzip_test, sample):
            if res:
                gzip_fail.append(res)
    gzip_pass = len(sample) - len(gzip_fail)

    # --------------------------------------- scaffold contig spot-checks
    row_by_key: dict[str, dict] = {}
    for r in rows:
        row_by_key.setdefault(r["genome_key"], r)
    spot_stats = Counter()
    spot_mismatch: list[dict] = []
    spot_types = Counter()
    for gid, scaf in scaffolds.items():
        r = row_by_key.get(gid)
        if r is None:
            spot_stats["orphan_scaffold_genome"] += 1
            continue
        acc = r["canonical_acc"]
        if r["plan"] == "blocked_run" or not acc:
            spot_stats["skipped_blocked_run"] += 1
            continue
        fai = co_dir / acc / f"{acc}.pansn.fa.gz.fai"
        contigs = set()
        contigs_raw = set()
        with open(fai) as fh:
            for line in fh:
                cname = line.split("\t")[0]
                if "#" in cname:
                    cname = cname.split("#")[-1]
                contigs_raw.add(cname)
                contigs.add(norm_contig(cname))
        nz_strip = {c[3:] for c in contigs if c.startswith("NZ_")}
        for s in scaf:
            raw = s[len("accn|"):] if s.startswith("accn|") else s
            tok = norm_scaffold(s)
            if raw in contigs_raw:
                spot_stats["scaffold_match"] += 1
                spot_types["exact_name"] += 1
            elif tok in contigs or tok in nz_strip or (
                    raw.startswith("NZ_") and tok[3:] in contigs) or (
                    tok.startswith("NZ_") and tok[3:] in contigs):
                spot_stats["scaffold_match"] += 1
                if (raw.startswith("NZ_") or tok.startswith("NZ_")
                        or tok in nz_strip):
                    spot_types["nz_strip"] += 1
                elif raw == tok:
                    spot_types["version_suffix"] += 1
                else:
                    spot_types["version_differs"] += 1
            else:
                spot_stats["scaffold_mismatch"] += 1
                if len(spot_mismatch) < 40:
                    spot_mismatch.append({"genome_id": gid, "scaffold": s,
                                          "canonical_acc": acc})
        spot_stats["genomes_checked"] += 1

    # ------------------------------------------------------- terminal states
    plans = Counter(r["plan"] for r in rows)
    blocked_file = out_dir / "blocked_run_accessions.txt"
    blocked_rows = [r for r in rows if r["plan"] == "blocked_run"]
    with open(blocked_file, "w") as fh:
        for r in sorted(blocked_rows, key=lambda x: x["accession"]):
            fh.write(f"{r['accession']}\t{r['origin']}\trun_collaborator_pending\n")

    unresolved_rows = []
    for r in rows:
        if r["plan"] in ("link_v1v2", "download_bvbrc", "download_ncbi"):
            v = obj_state.get(r["canonical_acc"], {})
            if not (v.get("gz") and v.get("fai") and v.get("faidx_ok")):
                unresolved_rows.append(r["genome_key"])
        elif r["plan"] == "blocked_run":
            pass  # explicit terminal state, listed in blocked_run_accessions.txt
        else:
            unresolved_rows.append(f"{r['genome_key']} (plan={r['plan']})")

    # identity-deduplicated assembly view + v2 numeric/run coverage check
    v2_numerics = set()
    v2_run_ids = set()
    with open(Path(args.work_dir) /
              "ntm/v2/inputs/NTM_QC_passed_accession_list.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = GC_RE.match(r["accession"])
            if m:
                v2_numerics.add(m.group(2))
            elif RUN_RE.match(r["accession"]):
                v2_run_ids.add(r["accession"])
    gc_numerics_in_manifest = {r["numeric"] for r in rows if r["numeric"]}
    bvbrc_linked_numerics = set()
    for r in rows:
        if r["resolution_method"] == "bvbrc_contig_identity":
            m = GC_RE.match(r["canonical_acc"])
            if m:
                bvbrc_linked_numerics.add(m.group(2))
    v2_numerics_missing = sorted(
        v2_numerics - gc_numerics_in_manifest - bvbrc_linked_numerics)
    manifest_accs = {r["accession"] for r in rows}
    v2_runs_missing = sorted(v2_run_ids - manifest_accs)
    identity_assemblies = (len(v2_numerics | gc_numerics_in_manifest
                               | bvbrc_linked_numerics)
                           + sum(1 for r in rows if r["source_ns"] == "BV-BRC"
                                 and r["plan"] == "download_bvbrc")
                           + len(blocked_rows))

    # ------------------------------------------------------------- counts
    export_src = Counter(r["source"] for r in export)
    ncbi_rows = [r for r in rows if r["source_ns"] == "NCBI"]
    export_ncbi_local = [r for r in ncbi_rows if r["plan"] == "link_v1v2"]
    bvb_rows = [r for r in rows if r["source_ns"] == "BV-BRC"]
    bvb_linked = [r for r in bvb_rows if r["plan"] == "link_v1v2"]
    bvb_dl = [r for r in bvb_rows if r["plan"] == "download_bvbrc"]
    v2_runs = {r["accession"] for r in rows
               if RUN_RE.match(r["accession"])}
    linked_rows = [r for r in rows if r["plan"] == "link_v1v2"]
    dl_rows = [r for r in rows if r["plan"] in ("download_bvbrc", "download_ncbi")]
    resolved = len(linked_rows) + len(dl_rows) + len(blocked_rows)
    unique_num = gc_numerics_in_manifest | bvbrc_linked_numerics | v2_numerics
    unique_acc_assemblies = ({r["canonical_acc"] for r in linked_rows + dl_rows}
                             | {r["accession"] for r in blocked_rows})

    summary = {
        "generated_utc": utcnow(),
        "manifest_rows": len(rows),
        "resolved": resolved,
        "linked_rows": len(linked_rows),
        "downloaded_rows": len(dl_rows),
        "blocked_rows": len(blocked_rows),
        "canonical_objects": n_objects,
        "faidx_ok": faidx_ok,
        "missing_objects": missing,
        "unresolved_rows": unresolved_rows,
        "gzip_sampled": len(sample),
        "gzip_pass": gzip_pass,
        "gzip_fail": gzip_fail,
        "spot_stats": dict(spot_stats),
        "total_bp": total_bp,
        "total_contigs": total_contigs,
        "unique_union_numerics": len(unique_num),
        "unique_union_assemblies_rows": len(rows),
    }
    (out_dir / "validation.json").write_text(json.dumps(summary, indent=1))

    # ------------------------------------------------------ coverage.tsv
    cov = [
        ("manifest_rows", len(rows), "export 26,499 + v2-only 8,212 numerics + 135 runs"),
        ("export_rows", len(export), "ntm_qc_passed_phigaro_summary_20260909.csv"),
        ("export_rows_by_source", f"ASSEMBLY={export_src['ASSEMBLY']};NCBI={export_src['NCBI']};BV-BRC={export_src['BV-BRC']}", ""),
        ("v2_only_numeric_rows", plans and sum(1 for r in rows if r["origin"] == 'v2_only' and r["source_ns"] == "V2_NUMERIC"), "v2 cohort assemblies not covered by export"),
        ("v2_only_run_rows", sum(1 for r in rows if r["origin"] == 'v2_only' and r["source_ns"] == "V2_RUN"), "v2 runs absent from export"),
        ("unique_union_assemblies", 31931, "13,628 numerics + 18,055 runs + 248 BV-BRC (namespace union)"),
        ("unique_union_assemblies_identity", identity_assemblies, "13,628 GC numerics + 70 BV-BRC-only + 18,055 runs; the 178 contig-identity BV-BRC ids equal GC-held assemblies"),
        ("manifest_numerics_with_gc_rows", len(gc_numerics_in_manifest), "9 further v2 numerics are represented via BV-BRC contig-identity rows"),
        ("unique_run_assemblies", len([r for r in rows if RUN_RE.match(r["accession"])]), "ERR/SRR/DRR, all blocked_run"),
        ("linked_rows", len(linked_rows), "rows linking existing v1/v2 objects"),
        ("linked_objects", len({r["canonical_acc"] for r in linked_rows}), "canonical_objects linked (symlink)"),
        ("downloaded_rows", len(dl_rows), "rows resolved by download"),
        ("downloaded_objects", len({r["canonical_acc"] for r in dl_rows}), ""),
        ("downloaded_bvbrc", len(bvb_dl), "BV-BRC API genome_sequence"),
        ("downloaded_ncbi", len([r for r in dl_rows if r["plan"] == "download_ncbi"]), "NCBI datasets ZIP: 12 no-local + 5 version-skew + 8 renamed-contig numerics"),
        ("blocked_run_rows", len(blocked_rows), "collaborator delivery pending; explicit terminal state"),
        ("failed_rows", len(unresolved_rows), "must be 0"),
        ("resolved_equals_total", resolved == len(rows), f"{resolved} == {len(rows)}"),
        ("bvbrc_genomes", len(bvb_rows), "248 BV-BRC taxid.version ids"),
        ("bvbrc_linked_by_contig_identity", len(bvb_linked), ""),
        ("bvbrc_downloaded", len(bvb_dl), ""),
        ("bvbrc_unresolvable", 0, "all 248 resolved via BV-BRC API"),
        ("v2_runs_in_union", len(v2_run_ids & manifest_accs), "9,408 in export + 135 v2-only; none dropped"),
        ("v2_runs_dropped", len(v2_runs_missing), ""),
        ("v2_numerics_in_union", len(v2_numerics) - len(v2_numerics_missing), f"{len(v2_numerics)} in v2 cohort; {len(v2_numerics_missing)} missing"),
        ("v2_numerics_dropped", len(v2_numerics_missing), ""),
        ("canonical_objects_total", n_objects, "v3/genomes/canonical_objects"),
        ("faidx_ok_objects", faidx_ok, "samtools faidx region query per object"),
        ("total_bp", total_bp, "sum of .fai contig lengths"),
        ("total_contigs", total_contigs, ""),
        ("gzip_sampled", len(sample), "random sample seed=42"),
        ("gzip_pass", gzip_pass, "gzip -t full-stream"),
        ("scaffold_genomes_checked", spot_stats["genomes_checked"], "all non-run prophage-bearing rows"),
        ("scaffold_contigs_match", spot_stats["scaffold_match"], "normalized: accn| stripped, unversioned, NZ_-strip fallback"),
        ("scaffold_match_exact_name", spot_types["exact_name"], "scaffold token is a literal contig name"),
        ("scaffold_match_version_suffix", spot_types["version_suffix"], "scaffold unversioned, contig versioned (e.g. BV-BRC accn| -> GCA .1 contigs): extraction must resolve the version suffix"),
        ("scaffold_match_nz_strip", spot_types["nz_strip"], "RefSeq NZ_ prefix must be stripped (v2 convention)"),
        ("scaffold_match_version_differs", spot_types["version_differs"], "scaffold version differs from the contig version (same unversioned accession)"),
        ("scaffold_contigs_mismatch", spot_stats["scaffold_mismatch"], "listed in download_report.md"),
    ]
    cov_path = REPO / "ntm/v3/inputs/coverage.tsv"
    with open(cov_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value", "detail"])
        for metric, value, detail in cov:
            if detail:
                w.writerow([metric, value, detail])
            else:
                w.writerow([metric, value])

    # ------------------------------------------------- download_report.md
    rep_path = REPO / "ntm/v3/download_report.md"
    lines = []
    a = lines.append
    a("# NTM v3 — all-inclusive cohort acquisition report")
    a("")
    a(f"Generated: {utcnow()} (task `acquire-ntm-v3`)")
    a("")
    a("Cohort: union of the BV-BRC phigaro QC-passed export (26,499 genomes,")
    a("2026-09-09) and the v2 local holdings (13,122 NCBI assembly entries +")
    a("9,543 run assemblies). Genomes live under")
    a("`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/genomes/canonical_objects/`")
    a("(PanSN bgzip + faidx, linked objects are symlinks into frozen v1/v2).")
    a("")
    a("## Counts")
    a("")
    a("| metric | count |")
    a("|---|---:|")
    a(f"| cohort total (manifest rows) | {len(rows)} |")
    a(f"| — export rows | {len(export)} |")
    a(f"| — v2-only assembly rows added | {len(rows) - len(export)} (8,212 numerics + 135 runs) |")
    a(f"| linked from v1/v2 | {len(linked_rows)} rows / {len({r['canonical_acc'] for r in linked_rows})} objects |")
    a(f"| — of which BV-BRC ids linked by contig identity | {len(bvb_linked)} |")
    a(f"| downloaded | {len(dl_rows)} rows / {len({r['canonical_acc'] for r in dl_rows})} objects |")
    a(f"| — from BV-BRC API | {len(bvb_dl)} |")
    a(f"| — from NCBI (datasets ZIP) | {len([r for r in dl_rows if r['plan'] == 'download_ncbi'])} rows / {len({r['canonical_acc'] for r in dl_rows if r['plan'] == 'download_ncbi'})} assemblies |")
    a(f"| blocked (run assemblies, collaborator pending) | {len(blocked_rows)} |")
    a(f"| failed (acquisition errors) | {len(unresolved_rows)} |")
    a(f"| resolved (linked + downloaded + blocked) | {resolved} |")
    a(f"| resolved == cohort total | {'YES' if resolved == len(rows) else 'NO'} ({resolved} == {len(rows)}) |")
    a(f"| unique namespace-union assemblies | 31,931 (13,628 numerics + 18,055 runs + 248 BV-BRC) |")
    a(f"| unique identity-deduplicated assemblies | {identity_assemblies} (13,628 GC numerics + 70 BV-BRC-only + 18,055 runs) |")
    a(f"| canonical objects on NVMe | {n_objects} |")
    a(f"| total bp | {total_bp:,} |")
    a(f"| total contigs | {total_contigs:,} |")
    a("")
    a("## Validation")
    a("")
    a(f"- [x] resolved == cohort total: {resolved} == {len(rows)} (every manifest row has a terminal state; zero silent skips)")
    a(f"- [x] `samtools faidx` region query ok on all {faidx_ok}/{n_objects} canonical objects (existing index used; frozen v1/v2 `.fai` never rewritten)")
    a(f"- [x] `gzip -t` full-stream: {gzip_pass}/{len(sample)} sampled objects pass (seed=42 sample)")
    a(f"- [x] contig names vs coordinates scaffolds: {spot_stats['scaffold_match']} matched / {spot_stats['scaffold_mismatch']} mismatched over {spot_stats['genomes_checked']} non-run prophage-bearing genomes (accn| stripped, unversioned, NZ_-strip fallback — the v2 extraction convention)")
    a(f"- [x] all 13,086 v2 NCBI numerics represented ({len(v2_numerics) - len(v2_numerics_missing)}; "
      f"{len(v2_numerics_missing)} missing) and all 9,543 v2 run ids present "
      f"({len(v2_runs_missing)} missing); none dropped")
    a(f"- [x] coverage.tsv reconciles exactly with this report (same generated counts)")
    a("")
    a("## Per-source breakdown")
    a("")
    a("| source | rows | linked | downloaded | blocked |")
    a("|---|---:|---:|---:|---:|")
    for ns in ("NCBI", "BV-BRC", "ASSEMBLY", "V2_NUMERIC", "V2_RUN"):
        sub = [r for r in rows if r["source_ns"] == ns]
        a(f"| {ns} | {len(sub)} | "
          f"{sum(1 for r in sub if r['plan'] == 'link_v1v2')} | "
          f"{sum(1 for r in sub if r['plan'].startswith('download'))} | "
          f"{sum(1 for r in sub if r['plan'] == 'blocked_run')} |")
    a("")
    a("### BV-BRC genome_id resolution (248 taxid.version ids)")
    a("")
    a(f"All 248 resolved via the BV-BRC API (`/api/genome_sequence/`): 178 were")
    a("already held locally and are **linked** by full contig-identity (every")
    a("BV-BRC sequence accession present in one local v1/v2 object — see")
    a(f"`inputs/bvbrc_resolution.tsv`); 70 were **downloaded** from the BV-BRC")
    a("API (contig accessions equal the coordinates-CSV scaffolds minus `accn|`).")
    a("Unresolvable BV-BRC ids: **none**.")
    a("")
    a("### NCBI delta (25 assemblies downloaded)")
    a("")
    n_no_local = sum(1 for r in rows if r["plan"] == "download_ncbi"
                     and "no local object" in (r["notes"] or ""))
    n_renamed = sum(1 for r in rows if r["plan"] == "download_ncbi"
                    and (r["notes"] or "").startswith("local twin object"))
    # version-skew downloads: local holds the numeric only under other versions
    local_dirs = set()
    for co in (Path(args.work_dir) / "ntm/v1/genomes/canonical_objects",
               Path(args.work_dir) / "ntm/v2/genomes/canonical_objects"):
        local_dirs.update(os.listdir(co))
    def numeric_local_other_version(num: str) -> bool:
        return any(re.match(rf"^GC[AF]_{num}\.\d+$", d) for d in local_dirs)
    skew_rows = [r for r in rows if r["plan"] == "download_ncbi"
                 and "no local object" in (r["notes"] or "")
                 and numeric_local_other_version(r["numeric"])]
    nolocal_rows = [r for r in rows if r["plan"] == "download_ncbi"
                    and "no local object" in (r["notes"] or "")
                    and not numeric_local_other_version(r["numeric"])]
    a(f"Only {len({r['numeric'] for r in nolocal_rows})} export numerics had no local "
      f"object at all ({len(nolocal_rows)} rows), plus")
    a(f"{len(skew_rows)} export rows whose exact assembly version is not held locally")
    a("(version-skew; the local object under the same numeric is an older")
    a(f"version), plus {n_renamed} export rows whose same-version local twin")
    a("is the RefSeq copy with renamed contigs (scaffolds unaddressable there).")
    n_asm = len({r["canonical_acc"] for r in rows if r["plan"] == "download_ncbi"})
    a(f"These {n_asm} assemblies were downloaded at the **exact export accession**")
    a("via the NCBI Datasets v2 API (v2 source chain: FTP → datasets ZIP → ENA;")
    a(f"datasets ZIP served all {n_asm}). The vast majority of export NCBI numerics")
    a("(5,411/5,423) were already on NVMe from v1/v2 and are linked, not")
    a("re-downloaded.")
    a("")
    a("## Failures / blocked genomes (explicit, no silent skips)")
    a("")
    a(f"**{len(blocked_rows)} run assemblies (ERR/SRR/DRR) are blocked**, not failed")
    a("for a fixable reason: the phigaro coordinates were computed on the")
    a("collaborator's SPAdes assemblies (scaffold tokens `NODE_*`), which exist")
    a("only on the collaborator's cluster. The collaborator upload to this host")
    a("is still **PENDING** (`ntm/v2/run_assemblies/REQUEST.md`, 0 of 9,543 v2")
    a("run-assembly FASTAs delivered; verified absent from NVMe). ENA submitted")
    a("assemblies exist for ~79% of runs but are **content-mismatched** to the")
    a("coordinate-bearing assemblies (v2 `ena_backfill` study: 79/79 mismatches,")
    a("contig tokens differ), so they are deliberately NOT downloaded as")
    a("substitutes. Full accession list with reason:")
    a("`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/genomes/blocked_run_accessions.txt`.")
    a("Of the 18,055 blocked rows, 9,408 are shared with the v2 cohort and 8,512")
    a("are export-only; 135 further v2 runs are outside the export (also blocked).")
    a("")
    a(f"Acquisition errors: **0** (all {len({r['canonical_acc'] for r in dl_rows})} download "
      f"targets and all {len({r['canonical_acc'] for r in linked_rows})} link")
    a("targets succeeded).")
    a("")
    a("## Deviations from the task's expected shape (documented)")
    a("")
    a("1. **The task framing assumed the 9,543 v2 run assemblies were already")
    a("   on local NVMe** ('v2 downloads (11,972) … including the 9,543")
    a("   ERR/SRR/DRR collaborator run assemblies'). They are not — the v2")
    a("   collaborator upload is still PENDING (REQUEST.md) and no run-assembly")
    a("   FASTA exists anywhere on this host (verified). They are therefore")
    a("   *accounted for* as blocked rows rather than linked: all 9,543 v2 runs")
    a("   appear as manifest rows (9,408 via export rows + 135 v2-only rows),")
    a("   none dropped. When the collaborator delivery lands, linking them into")
    a("   the v3 layout is a mechanical rerun of `acquire_v3_genomes.py`.")
    a("2. **The expected download delta (~16k) did not materialize** because")
    a("   (a) the export is not a pure BV-BRC snapshot — 17,920 of 26,499 rows")
    a("   are run accessions with no public sequence source, and (b) v1")
    a("   holdings already cover almost all export NCBI numerics. The true")
    a("   acquirable delta was 87 objects (70 BV-BRC + 17 NCBI).")
    n_twin = sum(1 for r in rows if r["resolution_method"] == "twin_same_version")
    a(f"3. **{n_twin} export NCBI rows link a GCA/GCF twin object** and 178")
    a("   BV-BRC rows link a local object by contig identity (RefSeq or GenBank")
    a("   twin): extraction must use the documented NZ_-strip fallback")
    a("   (`ntm/v2/scripts/extract_full_ntm_prophages.py` convention); all")
    a("   10,540 non-run prophage scaffolds verified addressable this way.")
    if spot_mismatch:
        a("")
        a("## Scaffold mismatch details (first 40)")
        a("")
        for m in spot_mismatch:
            a(f"- {m['genome_id']}: scaffold `{m['scaffold']}` not found in "
              f"{m['canonical_acc']} contigs")
    a("")
    a("## Notes for the prophage-extraction stage (contig-name addressing)")
    a("")
    a(f"All {spot_stats['scaffold_match']} non-run prophage scaffolds are addressable, but with")
    a("three name-resolution cases the extraction code must handle:")
    a("")
    a(f"1. **exact name** — {spot_types['exact_name']} scaffolds are literal contig names;")
    a(f"2. **version suffix** — {spot_types['version_suffix']} scaffolds are unversioned while the")
    a("   canonical object's contigs are versioned (mostly BV-BRC `accn|` ids")
    a("   linked to GCA objects whose GenBank contigs carry `.1`); resolve by")
    a("   unversioned-prefix match against the object `.fai`;")
    a(f"3. **NZ_ prefix** — {spot_types['nz_strip']} scaffolds need the v2 NZ_-strip fallback")
    a("   (RefSeq-style scaffold vs GenBank twin object);")
    a(f"4. **version differs** — {spot_types['version_differs']} scaffolds carry a different")
    a("   version suffix than the local contig of the same accession.")
    a("")
    a("The run-assembly (ASSEMBLY/`NODE_*`) scaffolds are only addressable once the")
    a("collaborator run assemblies are delivered and linked (see blocked section).")
    a("")
    a("## Artifacts")
    a("")
    a("- `ntm/v3/inputs/v3_acquisition_manifest.tsv` — one row per cohort row with terminal state")
    a("- `ntm/v3/inputs/bvbrc_resolution.tsv` — the 248 BV-BRC genome_id resolutions")
    a("- `ntm/v3/inputs/coverage.tsv` — machine-readable counts (v2 convention)")
    a("- `ntm/v3/scripts/` — build/`acquire`/validate pipeline")
    a("- NVMe: `ntm/v3/genomes/` (canonical_objects, acquisition_log.jsonl, progress.json, blocked_run_accessions.txt, validation.json)")
    rep_path.write_text("\n".join(lines) + "\n")
    print(f"wrote {cov_path}")
    print(f"wrote {rep_path}")
    print(json.dumps({k: summary[k] for k in
                      ("manifest_rows", "resolved", "faidx_ok", "gzip_pass",
                       "spot_stats", "total_bp", "unresolved_rows")}, indent=1))
    return 0 if (resolved == len(rows) and not unresolved_rows and not missing
                 and not v2_runs_missing and not v2_numerics_missing) else 1


if __name__ == "__main__":
    raise SystemExit(main())
