#!/usr/bin/env python3
"""
NTM v3 — build the all-inclusive acquisition manifest (export ∪ v2 holdings).

Cohort definition (v3 = union):
  * every row of the BV-BRC phigaro QC-passed export
    (ntm/v3/inputs/ntm_qc_passed_phigaro_summary_20260909.csv, 26,499 rows)
  * plus every v2 cohort assembly not already covered by an export row
    (v2 accession list: 13,122 NCBI assembly entries + 9,543 run assemblies)

Identifier namespaces in the export (see inputs/INPUTS.md §4):
  ASSEMBLY  ERR/SRR/DRR run accessions (collaborator SPAdes assemblies)
  NCBI      GCA_/GCF_ assembly accessions (GCA/GCF twins collapse by numeric)
  BV-BRC    BV-BRC genome_id (taxid.version), resolved via the BV-BRC API

Resolution per manifest row -> one of:
  plan=link_v1v2       object already on NVMe (v1/v2 canonical_objects)
  plan=download_ncbi   exact assembly version downloadable from NCBI
  plan=download_bvbrc  BV-BRC genome not held locally -> BV-BRC API sequences
  plan=blocked_run     run assembly: collaborator files not delivered locally
                       (v2 REQUEST.md still PENDING); ENA substitutions are
                       content-mismatched (v2 ena_backfill study) and are NOT
                       used. Explicit terminal state, never a silent skip.

Local-object matching rules for export NCBI rows (in order):
  1. exact accession dir (prefix+version)                      -> exact_acc
  2. same version, other GCA/GCF prefix (twin assembly)          -> twin_same_version
  3. local object under a different version                      -> download exact
  4. no local object for the numeric                             -> download

BV-BRC genomes: all sequence accessions (BV-BRC genome_sequence API) must be
found (unversioned) in one local object -> link (bvbrc_contig_identity);
otherwise download from the BV-BRC API (bvbrc_api_sequence). The BV-BRC
sequence accessions are exactly the coordinates-CSV scaffolds minus the
"accn|" prefix, so downloaded contig names address the prophage coordinates.

Outputs (repo, committed):
  ntm/v3/inputs/v3_acquisition_manifest.tsv.gz  one row per cohort row (gzip -n;
                                               deterministic; plain copy kept on
                                               NVMe — repo precedent: large
                                               tables are committed gzipped)
  ntm/v3/inputs/bvbrc_resolution.tsv           the 248 BV-BRC resolutions
Outputs (NVMe):
  {work}/ntm/v3/genomes/manifest_summary.json       machine-readable summary
  {work}/ntm/v3/genomes/v3_acquisition_manifest.tsv plain-text manifest copy
  {work}/ntm/v3/scratch/bvbrc_resolution.json        full API cache (no sequences)

Usage:
  python3 ntm/v3/scripts/build_v3_acquisition_manifest.py \
      [--work-dir /mnt/nvme3n1/erikg/phind-genome-work] [--reuse-bvbrc-cache]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]  # .../phind
SUMMARY_CSV = REPO / "ntm/v3/inputs/ntm_qc_passed_phigaro_summary_20260909.csv"
COORDS_CSV = REPO / "ntm/v3/inputs/ntm_qc_passed_phigaro_coordinates_20260909.csv"
GC_RE = re.compile(r"^GC([AF])_([0-9]{9})\.([0-9]+)$")
RUN_RE = re.compile(r"^(ERR|SRR|DRR)(\d+)$")
BVBRC_API = "https://www.bv-brc.org/api"
USER_AGENT = "phind-ntm-v3-acquire/1.0"
UNVERSION_RE = re.compile(r"\.\d+$")


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def unversion(token: str) -> str:
    return UNVERSION_RE.sub("", token)


def http_json(url: str, retries: int = 4, timeout: int = 90) -> list | dict:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read())
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"BV-BRC API failed after {retries} tries: {url}\n{last}")


def bvbrc_sequences(genome_id: str) -> list[dict]:
    """Sequence metadata (no sequence payload) for one BV-BRC genome."""
    q = urllib.parse.quote(genome_id)
    url = (f"{BVBRC_API}/genome_sequence/?eq(genome_id,{q})"
           "&select(sequence_id,accession,length,sequence_md5,sequence_type,topology)")
    return http_json(url)


def load_local_inventory(v1_co: Path, v2_co: Path) -> tuple[dict, dict]:
    """Scan v1/v2 canonical_objects.

    Returns (by_acc, by_num):
      by_acc: accession -> {"dir": name, "gz": path, "co": "v1"|"v2", "fai": bool}
      by_num: 9-digit numeric -> sorted list of accessions
    v2 objects win over v1 for the same accession (v2 is the newer freeze).
    """
    by_acc: dict[str, dict] = {}
    for label, co in (("v1", v1_co), ("v2", v2_co)):
        for d in sorted(os.listdir(co)):
            if not GC_RE.match(d):
                continue
            gz = co / d / f"{d}.pansn.fa.gz"
            fai = co / d / f"{d}.pansn.fa.gz.fai"
            if label == "v1" and d in by_acc:
                continue
            by_acc[d] = {"dir": d, "gz": str(gz), "co": label,
                         "has_gz": gz.exists(), "has_fai": fai.exists()}
    by_num: dict[str, list[str]] = defaultdict(list)
    for acc in by_acc:
        m = GC_RE.match(acc)
        by_num[m.group(2)].append(acc)
    for num in by_num:
        by_num[num].sort()
    return by_acc, by_num


def choose_local(rows_accessions: list[str], by_acc: dict, by_num: dict) -> tuple[str | None, str]:
    """Pick the local object for one export NCBI row set (per accession).

    Returns (canonical_acc | None, method).
    Preference: exact accession; same version other prefix (GCA first);
    then nothing (caller decides download — version-skew is downloaded, not
    silently substituted with a different version).
    """
    for acc in rows_accessions:
        if acc in by_acc and by_acc[acc]["has_gz"]:
            return acc, "exact_acc"
    for acc in rows_accessions:
        m = GC_RE.match(acc)
        cands = by_num.get(m.group(2), [])
        same_ver = [c for c in cands
                    if GC_RE.match(c).group(3) == m.group(3)
                    and by_acc[c]["has_gz"]]
        if same_ver:
            gca = [c for c in same_ver if c.startswith("GCA_")]
            pick = sorted(gca)[0] if gca else sorted(same_ver)[0]
            return pick, "twin_same_version"
    return None, ""


def load_local_contig_index(v1_co: Path, v2_co: Path) -> dict[str, set[str]]:
    """unversioned contig name -> set of local object accessions (from .fai)."""
    idx: dict[str, set[str]] = defaultdict(set)
    for co in (v1_co, v2_co):
        for d in sorted(os.listdir(co)):
            fai = co / d / f"{d}.pansn.fa.gz.fai"
            if not fai.exists():
                continue
            with open(fai) as fh:
                for line in fh:
                    name = line.split("\t", 1)[0]
                    if "#" in name:
                        name = name.split("#")[-1]
                    idx[unversion(name)].add(d)
    return idx


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", default="/mnt/nvme3n1/erikg/phind-genome-work")
    ap.add_argument("--reuse-bvbrc-cache", action="store_true")
    args = ap.parse_args()

    work = Path(args.work_dir)
    v1_co = work / "ntm/v1/genomes/canonical_objects"
    v2_co = work / "ntm/v2/genomes/canonical_objects"
    v2_acc_list = work / "ntm/v2/inputs/NTM_QC_passed_accession_list.tsv"
    out_dir = work / "ntm/v3/genomes"
    scratch = work / "ntm/v3/scratch"
    out_dir.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)

    manifest_path = REPO / "ntm/v3/inputs/v3_acquisition_manifest.tsv.gz"
    manifest_nvme = out_dir / "v3_acquisition_manifest.tsv"
    bvbrc_tsv = REPO / "ntm/v3/inputs/bvbrc_resolution.tsv"
    cache_path = scratch / "bvbrc_resolution.json"

    # ------------------------------------------------------------------ load
    export_rows = list(csv.DictReader(open(SUMMARY_CSV, newline="")))
    assert len(export_rows) == 26499, len(export_rows)

    by_acc, by_num = load_local_inventory(v1_co, v2_co)
    print(f"local objects: {len(by_acc)} accession dirs "
          f"({sum(1 for v in by_acc.values() if v['co'] == 'v2')} v2 / "
          f"{sum(1 for v in by_acc.values() if v['co'] == 'v1')} v1-only)")

    # per-genome scaffolds from the coordinates CSV (identity cross-checks)
    scaffolds: dict[str, set[str]] = defaultdict(set)
    with open(COORDS_CSV, newline="") as fh:
        for row in csv.DictReader(fh):
            scaffolds[row["genome_id"]].add(row["scaffold"])

    # ------------------------------------------------- export row resolution
    manifest: list[dict] = []
    covered_versions: set[tuple[str, str]] = set()  # (numeric, version) linked by export rows
    _contig_cache: dict[str, set[str]] = {}
    # group NCBI rows per numeric for twin handling of downloads
    ncbi_rows_by_num: dict[str, list[dict]] = defaultdict(list)

    for row in export_rows:
        gid, src, acc = row["genome_id"], row["source"], row["accession"]
        mrow = {
            "genome_key": gid,
            "origin": "export",
            "source_ns": src,
            "accession": acc,
            "numeric": "",
            "species": row["species"],
            "prophage_count": row["prophage_count"],
            "canonical_acc": "",
            "plan": "",
            "resolution_method": "",
            "link_gz": "",
            "bvbrc_contigs": "",
            "notes": "",
        }
        if src == "ASSEMBLY":
            assert RUN_RE.match(acc), acc
            mrow["plan"] = "blocked_run"
            mrow["resolution_method"] = "run_collaborator_pending"
            mrow["notes"] = ("collaborator run assembly not delivered locally "
                             "(ntm/v2/run_assemblies/REQUEST.md PENDING); "
                             "ENA substitute assemblies are content-mismatched")
        elif src == "NCBI":
            m = GC_RE.match(acc)
            assert m, acc
            mrow["numeric"] = m.group(2)
            ncbi_rows_by_num[m.group(2)].append(mrow)
        elif src == "BV-BRC":
            mrow["plan"] = "download_bvbrc"   # finalized after API resolution
        else:
            raise AssertionError(src)
        manifest.append(mrow)

    # NCBI rows: resolve per numeric (twins share the chosen object only when
    # downloading; linked rows keep their own exact/twin object per row)
    download_plan: dict[str, str] = {}  # numeric -> chosen accession to download

    def object_contigs(acc: str) -> set[str]:
        """Normalized contig tokens of a local object (cached)."""
        if acc not in _contig_cache:
            toks: set[str] = set()
            fai = by_acc[acc]["gz"] + ".fai"
            if os.path.exists(fai):
                with open(fai) as fh:
                    for line in fh:
                        name = line.split("\t", 1)[0]
                        if "#" in name:
                            name = name.split("#")[-1]
                        tok = unversion(name)
                        toks.add(tok)
                        if tok.startswith("NZ_"):
                            toks.add(tok[3:])
            _contig_cache[acc] = toks
        return _contig_cache[acc]

    def scaffolds_addressable(mrow: dict, acc: str) -> bool:
        """Every coordinates-CSV scaffold of this genome must be nameable in
        the object (normalized both sides, NZ_-strip fallback)."""
        scaf = scaffolds.get(mrow["genome_key"])
        if not scaf:
            return True
        toks = object_contigs(acc)
        for s in scaf:
            t = s[len("accn|"):] if s.startswith("accn|") else s
            t = unversion(t)
            if t in toks or (t.startswith("NZ_") and t[3:] in toks):
                continue
            return False
        return True

    for num, rows_ in ncbi_rows_by_num.items():
        # per-row link first
        for mrow in rows_:
            acc = mrow["accession"]
            pick, method = choose_local([acc], by_acc, by_num)
            if pick and not scaffolds_addressable(mrow, pick):
                # local twin object has renamed contigs (e.g. RefSeq NC_ names);
                # the prophage coordinates are unaddressable there -> download
                # the exact export accession instead
                mrow["notes"] = ("local twin object lacks the scaffold contig "
                                 "names (RefSeq renaming); downloading the "
                                 "exact export accession")
                pick = None
            if pick:
                mrow["plan"] = "link_v1v2"
                mrow["canonical_acc"] = pick
                mrow["resolution_method"] = method
                mrow["link_gz"] = by_acc[pick]["gz"]
                pm = GC_RE.match(pick)
                covered_versions.add((pm.group(2), pm.group(3)))
        pending = [r for r in rows_ if not r["canonical_acc"]]
        if pending:
            # one downloaded object per numeric, GCA-preferred (v2 convention)
            gca = sorted(r["accession"] for r in pending if r["accession"].startswith("GCA_"))
            gcf = sorted(r["accession"] for r in pending if r["accession"].startswith("GCF_"))
            chosen = (gca or gcf)[0]
            for mrow in pending:
                mrow["plan"] = "download_ncbi"
                mrow["canonical_acc"] = chosen
                mrow["resolution_method"] = "ncbi_exact_version"
                if not mrow["notes"]:
                    mrow["notes"] = ("no local object for this exact assembly "
                                     "version; downloading the export accession")
            download_plan[num] = chosen
            cm = GC_RE.match(chosen)
            covered_versions.add((cm.group(2), cm.group(3)))

    # --------------------------------------------------- BV-BRC resolution
    bvbrc_ids = [r["accession"] for r in manifest if r["source_ns"] == "BV-BRC"]
    cache = {}
    if args.reuse_bvbrc_cache and cache_path.exists():
        cache = json.loads(cache_path.read_text())
    missing = [g for g in bvbrc_ids if g not in cache]
    print(f"BV-BRC genomes: {len(bvbrc_ids)}; API calls needed: {len(missing)}")
    for i, gid in enumerate(missing):
        seqs = bvbrc_sequences(gid)
        cache[gid] = seqs
        if (i + 1) % 50 == 0:
            print(f"  {i + 1}/{len(missing)}", flush=True)
        time.sleep(0.3)
    cache_path.write_text(json.dumps(cache, indent=1))

    contig_idx = load_local_contig_index(v1_co, v2_co)
    bvbrc_report: list[dict] = []
    bvbrc_linked = 0
    for mrow in manifest:
        if mrow["source_ns"] != "BV-BRC":
            continue
        gid = mrow["accession"]
        seqs = cache[gid]
        accs = [s.get("accession") or s.get("sequence_id") for s in seqs]
        if not accs or any(a is None for a in accs):
            mrow["plan"] = "unresolvable"
            mrow["resolution_method"] = "bvbrc_api_error"
            mrow["notes"] = f"BV-BRC API returned {len(seqs)} sequences, unusable"
            continue
        matched = None
        if accs:
            sets = [contig_idx.get(unversion(a), set()) for a in accs]
            inter = set.intersection(*sets) if sets else set()
            if inter:
                matched = sorted(inter)[0] if len(inter) == 1 else sorted(
                    a for a in inter if by_acc.get(a, {}).get("has_gz"))[0]
        bp = sum(s.get("length", 0) for s in seqs)
        rep = {"genome_id": gid, "n_seqs": len(seqs), "seq_accessions": ",".join(accs),
               "total_bp": bp, "local_match": matched or ""}
        if matched and by_acc[matched]["has_gz"]:
            mrow["plan"] = "link_v1v2"
            mrow["canonical_acc"] = matched
            mrow["resolution_method"] = "bvbrc_contig_identity"
            mrow["link_gz"] = by_acc[matched]["gz"]
            mrow["bvbrc_contigs"] = ",".join(accs)
            mrow["notes"] = "all BV-BRC sequence accessions present in local object"
            pm = GC_RE.match(matched)
            if pm:
                covered_versions.add((pm.group(2), pm.group(3)))
            bvbrc_linked += 1
        else:
            mrow["plan"] = "download_bvbrc"
            mrow["canonical_acc"] = gid
            mrow["resolution_method"] = "bvbrc_api_sequence"
            mrow["bvbrc_contigs"] = ",".join(accs)
        bvbrc_report.append(rep)

    # ------------------------------------------------------- v2-only rows
    v2_runs: set[str] = set()
    v2_entries: dict[str, list[str]] = defaultdict(list)  # numeric.version -> accessions
    with open(v2_acc_list, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = row["accession"]
            if RUN_RE.match(acc):
                v2_runs.add(acc)
                continue
            m = GC_RE.match(acc)
            assert m, acc
            v2_entries[m.group(2) + "." + m.group(3)].append(acc)

    export_runs = {r["accession"] for r in manifest if r["source_ns"] == "ASSEMBLY"}
    n_v2_runs_added = 0
    for run in sorted(v2_runs - export_runs):
        manifest.append({
            "genome_key": run, "origin": "v2_only", "source_ns": "V2_RUN",
            "accession": run, "numeric": "", "species": "", "prophage_count": "",
            "canonical_acc": "", "plan": "blocked_run",
            "resolution_method": "run_collaborator_pending",
            "link_gz": "", "bvbrc_contigs": "", "notes":
            "v2-only run assembly; collaborator files not delivered locally",
        })
        n_v2_runs_added += 1

    n_v2_num_added = 0
    for key in sorted(v2_entries):
        accs = v2_entries[key]
        gca = sorted(a for a in accs if a.startswith("GCA_"))
        chosen = (gca or sorted(accs))[0]
        mch = GC_RE.match(chosen)
        if (mch.group(2), mch.group(3)) in covered_versions:
            continue  # this assembly (numeric.version) is already covered by an
                      # export row (possibly under the other GCA/GCF prefix)
        info = by_acc.get(chosen)
        if info is None or not info["has_gz"]:
            raise RuntimeError(f"v2 cohort entry {chosen} has no local FASTA")
        manifest.append({
            "genome_key": chosen, "origin": "v2_only", "source_ns": "V2_NUMERIC",
            "accession": chosen, "numeric": GC_RE.match(chosen).group(2),
            "species": "", "prophage_count": "",
            "canonical_acc": chosen, "plan": "link_v1v2",
            "resolution_method": "v2_cohort_entry",
            "link_gz": info["gz"], "bvbrc_contigs": "",
            "notes": (f"v2 cohort assembly not in export (twins in v2 list: "
                      f"{','.join(sorted(set(accs) - {chosen}))})"),
        })
        mch2 = GC_RE.match(chosen)
        covered_versions.add((mch2.group(2), mch2.group(3)))
        n_v2_num_added += 1

    # ------------------------------------------------------------ outputs
    fields = ["genome_key", "origin", "source_ns", "accession", "numeric",
              "species", "prophage_count", "canonical_acc", "plan",
              "resolution_method", "link_gz", "bvbrc_contigs", "notes"]
    import gzip as _gzip
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=fields, delimiter="\t",
                       lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
    w.writeheader()
    for mrow in manifest:
        w.writerow(mrow)
    with open(manifest_path, "wb") as _raw:
        with _gzip.GzipFile(fileobj=_raw, mode="wb", compresslevel=9,
                            mtime=0) as fh:
            fh.write(buf.getvalue().encode("utf-8"))
    with open(manifest_nvme, "w", newline="") as fh:
        fh.write(buf.getvalue())

    with open(bvbrc_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome_id", "n_seqs",
                                           "seq_accessions", "total_bp",
                                           "local_match", "plan"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        for rep in bvbrc_report:
            rep = dict(rep)
            rep["plan"] = ("link_v1v2" if rep["local_match"] else "download_bvbrc")
            w.writerow(rep)

    plans = Counter(m["plan"] for m in manifest)
    methods = Counter(m["resolution_method"] for m in manifest)
    linked_accs = {m["canonical_acc"] for m in manifest if m["plan"] == "link_v1v2"}
    summary = {
        "generated_utc": utcnow(),
        "manifest_rows": len(manifest),
        "export_rows": 26499,
        "v2_only_numeric_rows": n_v2_num_added,
        "v2_only_run_rows": n_v2_runs_added,
        "plans": dict(plans),
        "resolution_methods": dict(methods),
        "unique_canonical_objects_linked": len(linked_accs),
        "ncbi_download_numerics": len(download_plan),
        "ncbi_download_accessions": sorted(download_plan.values()),
        "bvbrc_linked": bvbrc_linked,
        "bvbrc_download": plans.get("download_bvbrc", 0),
        "blocked_runs": plans.get("blocked_run", 0),
    }
    (out_dir / "manifest_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    print(f"wrote {manifest_path} ({len(manifest)} rows) + plain copy {manifest_nvme}")
    print(f"wrote {bvbrc_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
