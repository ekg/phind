#!/usr/bin/env python3
"""
fetch_sources.py — cache all raw source responses for the public
mycobacteriophage reference panel and record provenance.

Stages (all idempotent; cached artifacts are never re-downloaded):
  bulk        PhagesDB metadata TSV + bulk actinobacteriophage FASTA,
              INPHARED release table + genomes fasta.
  ncbi        E-utilities esearch (mycobacteriophage nucleotide records),
              batched efetch FASTA (versioned accessions), batched esummary.
  taxonomy    resolve host-genus taxids + Mycobacteriaceae lineage via
              NCBI taxonomy efetch (tiny, cached per genus).
  gb          GenBank flat files for NCBI-only records that have no host
              label from PhagesDB/INPHARED (bounded remainder; the set is
              computed with the same merge logic as the builder).
  manifest    (re)write cache/fetch_manifest.json with sha256, sizes and
              retrieval timestamps for every artifact.

Usage:
  python3 fetch_sources.py [--data-dir DIR] [--stage all|bulk|ncbi|taxonomy|gb|manifest]
"""

from __future__ import annotations

import argparse
import datetime as dt
import glob
import gzip
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import panel_lib as pl  # noqa: E402

DEFAULT_DATA_DIR = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/external_validation"

PHAGESDB_METADATA_URL = "https://phagesdb.org/data/?set=seq&type=full"
PHAGESDB_FASTA_URL = "https://phagesdb.org/media/Actinobacteriophages-All.fasta"
INPHARED_TABLE_URL = "https://s3.climb.ac.uk/millardlab-inphared/2026/7Apr2026_millardlab_website_table.txt.gz"
INPHARED_FASTA_URL = "https://s3.climb.ac.uk/millardlab-inphared/2026/7Apr2026_genomes.fa.gz"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_TERM = 'mycobacteriophage[All Fields] OR "Mycobacterium phage"[Title]'


def http_get(url: str, out: str, max_time: int = 1800, retries: int = 3) -> bool:
    for attempt in range(retries):
        r = subprocess.run(["curl", "-sS", "--max-time", str(max_time), "--retry", "2",
                            "-o", out + ".part", url])
        if r.returncode == 0 and os.path.getsize(out + ".part") > 0:
            os.replace(out + ".part", out)
            return True
        time.sleep(3 * (attempt + 1))
    sys.stderr.write(f"[fetch] FAILED: {url} -> {out}\n")
    return False


def http_post(url: str, params: dict, out: str, max_time: int = 600,
              retries: int = 3) -> bool:
    """POST form params (required for long efetch id lists)."""
    for attempt in range(retries):
        cmd = ["curl", "-sS", "--max-time", str(max_time), "-o", out + ".part", "-X", "POST"]
        for k, v in params.items():
            cmd += ["--data-urlencode", f"{k}={v}"]
        cmd.append(url)
        r = subprocess.run(cmd)
        if (r.returncode == 0 and os.path.getsize(out + ".part") > 0
                and not open(out + ".part", errors="ignore").read(200).lstrip().startswith("<?xml")):
            os.replace(out + ".part", out)
            return True
        time.sleep(5 * (attempt + 1))
    sys.stderr.write(f"[fetch] FAILED POST: {url} -> {out}\n")
    return False


def head_last_modified(url: str) -> str:
    r = subprocess.run(["curl", "-sI", "--max-time", "30", url],
                       capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.lower().startswith("last-modified:"):
            return line.split(":", 1)[1].strip()
    return ""


def now_utc() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def manifest_path(cache_dir: str) -> str:
    return os.path.join(cache_dir, "fetch_manifest.json")


def load_manifest(cache_dir: str) -> dict:
    p = manifest_path(cache_dir)
    return json.load(open(p)) if os.path.exists(p) else {}


def record(man: dict, key: str, url: str, path: str, retrieved_at: str = "",
           note: str = "", last_modified: str = "") -> None:
    man[key] = {
        "url": url,
        "file": os.path.abspath(path),
        "sha256": pl.sha256_file(path),
        "bytes": os.path.getsize(path),
        "retrieved_at": retrieved_at or now_utc(),
        "http_last_modified": last_modified,
        "note": note,
    }


def stage_bulk(cache_dir: str) -> dict:
    man = load_manifest(cache_dir)
    jobs = [
        ("phagesdb_metadata", PHAGESDB_METADATA_URL,
         os.path.join(cache_dir, "phagesdb", "phagesdb_sequenced_phages_full.tsv")),
        ("phagesdb_bulk_fasta", PHAGESDB_FASTA_URL,
         os.path.join(cache_dir, "phagesdb", "Actinobacteriophages-All.fasta")),
        ("inphared_table", INPHARED_TABLE_URL,
         os.path.join(cache_dir, "inphared", "7Apr2026_millardlab_website_table.txt.gz")),
        ("inphared_genomes_fa", INPHARED_FASTA_URL,
         os.path.join(cache_dir, "inphared", "7Apr2026_genomes.fa.gz")),
    ]
    for key, url, path in jobs:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if not os.path.exists(path):
            print(f"[bulk] downloading {key}: {url}", flush=True)
            ok = http_get(url, path)
            if not ok:
                raise RuntimeError(f"download failed: {url}")
            retrieved = now_utc()
            lm = head_last_modified(url)
        else:
            retrieved = man.get(key, {}).get("retrieved_at", "") or dt.datetime.fromtimestamp(
                os.path.getmtime(path), dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            lm = man.get(key, {}).get("http_last_modified", "")
        record(man, key, url, path, retrieved_at=retrieved,
               last_modified=lm,
               note="cached bulk artifact" if os.path.exists(path) else "")
        print(f"[bulk] {key}: {man[key]['bytes']} bytes sha256={man[key]['sha256'][:12]}…")
    return man


def stage_ncbi(cache_dir: str) -> dict:
    man = load_manifest(cache_dir)
    ncbi = os.path.join(cache_dir, "ncbi")
    os.makedirs(ncbi, exist_ok=True)

    esearch = os.path.join(ncbi, "esearch_mycobacteriophage.xml")
    ids_xml = os.path.join(ncbi, "esearch_ids.xml")
    if not (os.path.exists(esearch) and os.path.exists(ids_xml)):
        print("[ncbi] esearch", flush=True)
        http_get(f"{EUTILS}/esearch.fcgi?db=nucleotide&usehistory=y&retmax=0"
                 f"&term={urllib_quote(NCBI_TERM)}", esearch)
        time.sleep(0.6)
        http_get(f"{EUTILS}/esearch.fcgi?db=nucleotide&retmax=10000"
                 f"&term={urllib_quote(NCBI_TERM)}", ids_xml)
    webenv_qkey = os.path.join(ncbi, "webenv.txt")
    if not os.path.exists(webenv_qkey):
        t = open(esearch).read()
        m = re.search(r"<WebEnv>(\w+)</WebEnv>", t)
        q = re.search(r"<QueryKey>(\d+)</QueryKey>", t)
        open(webenv_qkey, "w").write(f"{m.group(1)}\n{q.group(1)}\n")

    uids = re.findall(r"<Id>(\d+)</Id>", open(ids_xml).read())
    print(f"[ncbi] {len(uids)} uids", flush=True)

    webenv, qkey = open(webenv_qkey).read().split()
    got_fa = len(glob.glob(os.path.join(ncbi, "efetch_fasta_*.txt")))
    for start in range(0, len(uids), 500):
        out = os.path.join(ncbi, f"efetch_fasta_{start:05d}.txt")
        if os.path.exists(out):
            continue
        print(f"[ncbi] efetch fasta {start}…", flush=True)
        http_get(f"{EUTILS}/efetch.fcgi?db=nucleotide&WebEnv={webenv}"
                 f"&query_key={qkey}&rettype=fasta&retmode=text"
                 f"&retstart={start}&retmax=500", out, max_time=600)
        time.sleep(0.6)
    for start in range(0, len(uids), 500):
        out = os.path.join(ncbi, f"esummary_{start:05d}.json")
        if os.path.exists(out):
            continue
        print(f"[ncbi] esummary {start}…", flush=True)
        http_get(f"{EUTILS}/esummary.fcgi?db=nuccore&query_key={qkey}"
                 f"&WebEnv={webenv}&retmode=json&retstart={start}&retmax=500",
                 out, max_time=600)
        time.sleep(0.6)

    retrieved = man.get("ncbi_esearch", {}).get("retrieved_at", "") or now_utc()
    record(man, "ncbi_esearch",
           f"{EUTILS}/esearch.fcgi?db=nucleotide&term=…&usehistory=y",
           esearch, retrieved_at=retrieved,
           note=f"term={NCBI_TERM!r}; uids={len(uids)}")
    fa_files = sorted(glob.glob(os.path.join(ncbi, "efetch_fasta_*.txt")))
    ss_files = sorted(glob.glob(os.path.join(ncbi, "esummary_*.json")))
    record(man, "ncbi_efetch", f"{EUTILS}/efetch.fcgi?db=nucleotide&rettype=fasta",
           os.path.join(ncbi, "efetch_fasta_00000.txt"),
           retrieved_at=retrieved,
           note=f"batched efetch: {len(fa_files)} batches of <=500 uids; "
                "versioned accessions in FASTA headers")
    record(man, "ncbi_esummary", f"{EUTILS}/esummary.fcgi?db=nuccore",
           os.path.join(ncbi, "esummary_00000.json"),
           retrieved_at=retrieved,
           note=f"{len(ss_files)} JSON batches")
    return man


def urllib_quote(s: str) -> str:
    import urllib.parse
    return urllib.parse.quote(s)


def distinct_host_genera(cache_dir: str) -> set[str]:
    """Host genus labels mentioned by any in-scope source (from caches only)."""
    genera = set()
    pdb_tsv = os.path.join(cache_dir, "phagesdb", "phagesdb_sequenced_phages_full.tsv")
    if os.path.exists(pdb_tsv):
        for r in pl.parse_phagesdb_tsv(pdb_tsv):
            g = pl.host_genus_of(r.get("Host", ""))
            if g:
                genera.add(g)
    pdb_fa = os.path.join(cache_dir, "phagesdb", "Actinobacteriophages-All.fasta")
    if os.path.exists(pdb_fa):
        with open(pdb_fa, errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    g = line[1:].split()[0]
                    if g in pl.ACTINOPHAGE_HOST_GENERA:
                        genera.add(g)
    inp_tsv = os.path.join(cache_dir, "inphared",
                           "7Apr2026_millardlab_website_table.txt.gz")
    if os.path.exists(inp_tsv):
        for row in pl.parse_inphared_table(inp_tsv):
            if pl.inphared_scope(row):
                g = pl.host_genus_of(row["host"])
                if g:
                    genera.add(g)
    genera |= pl.MYCOBACTERIACEAE_GENERA
    return genera


def stage_taxonomy(cache_dir: str) -> dict:
    man = load_manifest(cache_dir)
    out_path = os.path.join(cache_dir, "taxonomy_hosts.json")
    tax = json.load(open(out_path)) if os.path.exists(out_path) else {}
    genera = sorted(distinct_host_genera(cache_dir))
    print(f"[taxonomy] {len(genera)} distinct host genera", flush=True)
    changed = False
    for g in genera:
        key = g.lower()
        if key in tax:
            continue
        # esearch by Scientific Name, then efetch lineage
        r = subprocess.run(
            ["curl", "-s", "--max-time", "60", "--get",
             f"{EUTILS}/esearch.fcgi",
             "--data-urlencode", "db=taxonomy",
             "--data-urlencode", f"term={g}[Scientific Name]"],
            capture_output=True, text=True)
        m = re.search(r"<Id>(\d+)</Id>", r.stdout)
        time.sleep(0.4)
        if not m:
            tax[key] = {"taxid": "", "scientific_name": g, "rank": "no_match",
                        "family_taxid": "", "lineage": [],
                        "retrieved_at": now_utc()}
            changed = True
            continue
        taxid = m.group(1)
        r = subprocess.run(
            ["curl", "-s", "--max-time", "60", "--get",
             f"{EUTILS}/efetch.fcgi",
             "--data-urlencode", "db=taxonomy",
             "--data-urlencode", f"id={taxid}",
             "--data-urlencode", "retmode=xml"],
            capture_output=True, text=True)
        names = re.findall(r"<ScientificName>([^<]+)</ScientificName>", r.stdout)
        ranks = re.findall(r"<Rank>([^<]+)</Rank>", r.stdout)
        taxids = re.findall(r"<TaxId>(\d+)</TaxId>", r.stdout)
        # first Taxon block is the node itself; lineage follows
        lineage = names[1:] if names else []
        lineage_t = taxids[1:] if taxids else []
        family_taxid = ""
        for n, t in zip(lineage, lineage_t):
            if n == "Mycobacteriaceae":
                family_taxid = t
        tax[key] = {
            "taxid": taxid,
            "scientific_name": names[0] if names else g,
            "rank": ranks[0] if ranks else "",
            "family_taxid": family_taxid or (str(pl.MYCOBACTERIACEAE_FAMILY_TAXID)
                                             if g in pl.MYCOBACTERIACEAE_GENERA else ""),
            "lineage": lineage,
            "retrieved_at": now_utc(),
        }
        changed = True
        print(f"[taxonomy] {g} -> {taxid} ({tax[key]['rank']})")
        time.sleep(0.4)
    if changed or not os.path.exists(out_path):
        json.dump(tax, open(out_path, "w"), indent=1, sort_keys=True)
    record(man, "taxonomy_hosts",
           f"{EUTILS}/esearch.fcgi?db=taxonomy + efetch.fcgi?db=taxonomy",
           out_path,
           retrieved_at=min(v["retrieved_at"] for v in tax.values()) if tax else now_utc(),
           note="host-genus taxid + lineage lookups (cached per genus)")
    return man


def compute_merge(cache_dir: str) -> dict:
    pdb_tsv = os.path.join(cache_dir, "phagesdb", "phagesdb_sequenced_phages_full.tsv")
    pdb_fa = os.path.join(cache_dir, "phagesdb", "Actinobacteriophages-All.fasta")
    inp_tsv = os.path.join(cache_dir, "inphared",
                           "7Apr2026_millardlab_website_table.txt.gz")
    inp_fa = os.path.join(cache_dir, "inphared", "7Apr2026_genomes.fa.gz")
    ncbi_dir = os.path.join(cache_dir, "ncbi")
    man = load_manifest(cache_dir)

    print("[merge-prep] parsing phagesdb…", flush=True)
    pdb_meta = pl.parse_phagesdb_tsv(pdb_tsv)
    pdb_fasta, pdb_extras = pl.parse_phagesdb_fasta(pdb_fa)
    print("[merge-prep] parsing inphared table…", flush=True)
    inp_rows = pl.parse_inphared_table(inp_tsv)
    print("[merge-prep] parsing ncbi summaries…", flush=True)
    summ = pl.parse_ncbi_esummary(
        sorted(glob.glob(os.path.join(ncbi_dir, "esummary_*.json")))
        + sorted(glob.glob(os.path.join(ncbi_dir, "backfill_esummary_*.json"))))
    ncbi_fa = pl.parse_ncbi_fasta(
        sorted(glob.glob(os.path.join(ncbi_dir, "efetch_fasta_*.txt")))
        + sorted(glob.glob(os.path.join(ncbi_dir, "backfill_efetch_*.txt"))))
    taxonomy = {}
    tpath = os.path.join(cache_dir, "taxonomy_hosts.json")
    if os.path.exists(tpath):
        taxonomy = json.load(open(tpath))
    gb = pl.parse_gb_remainder(os.path.join(ncbi_dir, "gb_remainder.txt"))
    # inphared sequences only needed for in-scope rows
    wanted = {r["accession"] for r in inp_rows if pl.inphared_scope(r)}
    print(f"[merge-prep] streaming inphared fasta ({len(wanted)} in-scope)…", flush=True)
    inp_seqs = pl.load_inphared_seqs(inp_fa, wanted)
    print(f"[merge-prep] got {len(inp_seqs)} sequences", flush=True)
    return pl.merge_records(pdb_meta, pdb_fasta, pdb_extras, inp_rows, inp_seqs, summ,
                            ncbi_fa, taxonomy, gb, man)


def _valid_json_result(path: str) -> bool:
    try:
        json.load(open(path))["result"]
        return True
    except Exception:
        return False


def _external_accessions(cache_dir: str) -> set[str]:
    """Accessions claimed by PhagesDB/INPHARED rows (base-normalized)."""
    accs = set()
    pdb_tsv = os.path.join(cache_dir, "phagesdb", "phagesdb_sequenced_phages_full.tsv")
    if os.path.exists(pdb_tsv):
        for r in pl.parse_phagesdb_tsv(pdb_tsv):
            a = (r.get("Accession #") or "").strip()
            if a and a.lower() not in ("none", "n/a", "-") and re.match(r"^[A-Za-z]{1,6}_?\d+$", a):
                accs.add(pl.accession_base(a))
    inp_tsv = os.path.join(cache_dir, "inphared",
                           "7Apr2026_millardlab_website_table.txt.gz")
    if os.path.exists(inp_tsv):
        for row in pl.parse_inphared_table(inp_tsv):
            if pl.inphared_scope(row) and re.match(r"^[A-Za-z]{1,6}_?\d+(\.\d+)?$",
                                                   row["accession"]):
                accs.add(pl.accession_base(row["accession"]))
    return accs


def _ncbi_known_bases(cache_dir: str) -> set[str]:
    summ = pl.parse_ncbi_esummary(
        sorted(glob.glob(os.path.join(cache_dir, "ncbi", "esummary_*.json"))))
    return {pl.accession_base(r["accessionversion"]) for r in summ.values()}


def stage_backfill(cache_dir: str) -> dict:
    """Fetch NCBI records for PhagesDB/INPHARED accessions the esearch
    term missed (very recent deposits are not yet in the search index,
    but efetch/esummary-by-id resolve them)."""
    man = load_manifest(cache_dir)
    ncbi = os.path.join(cache_dir, "ncbi")
    missing = sorted(_external_accessions(cache_dir) - _ncbi_known_bases(cache_dir))
    print(f"[backfill] {len(missing)} external accessions not in esearch set",
          flush=True)
    missing_path = os.path.join(ncbi, "backfill_missing.txt")
    open(missing_path, "w").write("\n".join(missing) + "\n")

    found: list[str] = []
    for i in range(0, len(missing), 150):
        chunk = missing[i:i + 150]
        out = os.path.join(ncbi, f"backfill_esummary_{i:06d}.json")
        for attempt in range(4):
            if not _valid_json_result(out):
                if os.path.exists(out):
                    os.remove(out)
                http_post(f"{EUTILS}/esummary.fcgi",
                          {"db": "nuccore", "retmode": "json",
                           "id": ",".join(chunk)}, out)
                time.sleep(1.0)
            else:
                break
        try:
            d = json.load(open(out))["result"]
            for u in d.get("uids", []):
                found.append(d[u]["accessionversion"])
        except Exception:
            sys.stderr.write(f"[backfill] WARNING: unresolved esummary batch {i}\n")
    found = sorted(found)
    print(f"[backfill] {len(found)} resolvable by id", flush=True)
    json.dump(found, open(os.path.join(ncbi, "backfill_found.json"), "w"))

    for i in range(0, len(found), 150):
        chunk = found[i:i + 150]
        out = os.path.join(ncbi, f"backfill_efetch_{i:06d}.txt")
        if os.path.exists(out):
            continue
        if not http_post(f"{EUTILS}/efetch.fcgi",
                         {"db": "nucleotide", "rettype": "fasta",
                          "retmode": "text", "id": ",".join(chunk)},
                         out, max_time=900):
            raise RuntimeError(f"backfill efetch batch {i} failed")
        time.sleep(0.7)
    if found:
        record(man, "ncbi_backfill",
               f"{EUTILS}/esummary+efetch by accession id",
               os.path.join(ncbi, "backfill_efetch_000000.txt"),
               note=f"esearch-missed accessions resolved by id: {len(found)}; "
                    f"requested {len(missing)}")
    return man


def stage_gb(cache_dir: str) -> dict:
    """Fetch GenBank flat files for NCBI-only records lacking host labels."""
    man = load_manifest(cache_dir)
    merged = compute_merge(cache_dir)
    needs = merged["needs_gb_remainder"]
    gb_path = os.path.join(cache_dir, "ncbi", "gb_remainder.txt")
    print(f"[gb] {len(needs)} records need GenBank host metadata", flush=True)
    existing = pl.parse_gb_remainder(gb_path) if os.path.exists(gb_path) else {}
    delta = [n for n in needs if n not in existing]
    print(f"[gb] {len(delta)} not yet cached", flush=True)
    if delta:
        chunks = [",".join(delta[i:i + 150]) for i in range(0, len(delta), 150)]
        parts = []
        for c in chunks:
            out = gb_path + f".part{len(parts)}"
            if not http_post(f"{EUTILS}/efetch.fcgi",
                             {"db": "nucleotide", "id": c,
                              "rettype": "gb", "retmode": "text"}, out):
                raise RuntimeError(f"gb efetch failed for chunk of {len(c.split(','))}")
            parts.append(out)
            time.sleep(0.8)
        with open(gb_path, "a") as out:
            for p in parts:
                out.write("\n//\n")
                out.write(open(p, errors="replace").read())
                os.remove(p)
    if os.path.exists(gb_path):
        record(man, "ncbi_gb_remainder",
               f"{EUTILS}/efetch.fcgi?db=nucleotide&rettype=gb&retmode=text",
               gb_path,
               note=f"remainder GenBank records for host/topology qualifiers: "
                    f"{len(needs)} requested")
    return man


def write_manifest(cache_dir: str, man: dict) -> None:
    json.dump(man, open(manifest_path(cache_dir), "w"), indent=1, sort_keys=True)
    print(f"[manifest] {len(man)} artifacts -> {manifest_path(cache_dir)}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=DEFAULT_DATA_DIR)
    ap.add_argument("--stage", default="all",
                    choices=["all", "bulk", "ncbi", "taxonomy", "backfill",
                             "gb", "manifest"])
    args = ap.parse_args(argv)
    cache_dir = os.path.join(args.data_dir, "cache")
    os.makedirs(cache_dir, exist_ok=True)

    man = load_manifest(cache_dir)
    stages = (["bulk", "ncbi", "taxonomy", "backfill", "gb", "manifest"]
              if args.stage == "all" else [args.stage])
    for st in stages:
        if st == "bulk":
            man = stage_bulk(cache_dir) or man
        elif st == "ncbi":
            man = stage_ncbi(cache_dir) or man
        elif st == "taxonomy":
            man = stage_taxonomy(cache_dir) or man
        elif st == "backfill":
            man = stage_backfill(cache_dir) or man
        elif st == "gb":
            man = stage_gb(cache_dir) or man
        elif st == "manifest":
            pass
        if man:
            write_manifest(cache_dir, man)
    return 0


if __name__ == "__main__":
    sys.exit(main())
