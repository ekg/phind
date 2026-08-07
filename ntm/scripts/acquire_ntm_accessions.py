#!/usr/bin/env python3
"""
acquire_ntm_accessions.py — scope & list NTM (non-tuberculous Mycobacterium)
assemblies from NCBI Assembly via E-utilities.

NTM definition used here: every assembly under family *Mycobacteriaceae* —
which under the post-2018 Gupta reclassification spans Mycobacterium (s.s.),
Mycobacteroides (abscessus/chelonae), Mycolicibacterium (fortuitum/smegmatis),
Mycolicibacter and Mycolicibacillus — EXCEPT the *M. tuberculosis* complex
(MTC) and *M. leprae* / *M. lepromatosis*. Querying only "Mycobacterium" misses
the reclassified genera, so the family-level scope is required for completeness.

Outputs (in --outdir):
  ntm_accessions.txt         one assembly accession per line (GCF_/GCA_ vN)
  ntm_accession_manifest.tsv full per-assembly metadata
  species_summary.tsv        species -> count, median genome size, assembly levels
  scoping_report.json        run params + summary stats

Usage:
  python3 acquire_ntm_accessions.py --outdir /mnt/nvme3n1/erikg/phind-genome-work/ntm/v1/accessions
  # optional: restrict to assembly_level Complete/Chromosome/Scaffold
  python3 acquire_ntm_accessions.py --outdir ... --min-level Scaffold
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "phind-ntm-scoping/1.0 (research)"

# Assembly levels in increasing contiguity; --min-level keeps this and above.
LEVEL_ORDER = ["Contig", "Scaffold", "Chromosome", "Complete Genome"]


def http_get(url: str, retries: int = 4) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"GET failed after {retries} tries: {url}\n{last}")


def esearch_uids(query: str) -> list[str]:
    uids: list[str] = []
    retmax = 500
    start = 0
    while True:
        q = urllib.parse.urlencode({
            "db": "assembly", "term": query, "retmode": "json",
            "retmax": retmax, "retstart": start, "usehistory": "n",
        })
        data = json.loads(http_get(f"{EUTILS}/esearch.fcgi?{q}"))
        res = data.get("esearchresult", {})
        batch = res.get("idlist", [])
        total = int(res.get("count", "0"))
        uids.extend(batch)
        start += retmax
        if start >= total or not batch:
            break
        time.sleep(0.34)  # respect 3 req/s
    return uids


def parse_meta(meta_text: str) -> dict:
    """Parse the CDATA <Meta> blob for total_length / contig_count / assembly-status."""
    m = {}
    try:
        meta = ET.fromstring("<root>" + meta_text.strip() + "</root>")
    except ET.ParseError:
        return m
    for st in meta.iter("Stat"):
        cat = st.get("category")
        if cat in ("total_length", "contig_count", "scaffold_count"):
            m[cat] = int((st.text or "0").strip() or 0)
    a = meta.find("assembly-status")
    if a is not None and a.text:
        m["assembly-status"] = a.text.strip()
    return m


def parse_docsum(xml_bytes: bytes) -> list[dict]:
    """Parse esummary assembly DocumentSummarySet -> records."""
    out = []
    root = ET.fromstring(xml_bytes)
    for ds in root.iter("DocumentSummary"):
        organism = ds.findtext("Organism") or ""
        species = (ds.findtext("SpeciesName") or "").strip()
        rec = {
            "uid": ds.findtext("Uid") or ds.attrib.get("uid", ""),
            "assembly_acc": ds.findtext("AssemblyAccession") or "",
            "biosample_acc": ds.findtext("BioSampleAccn") or "",
            "organism": organism,
            "species": species.split(" subspecies ")[0].strip() or species,
            "genus": species.split()[0] if species else "",
            "assembly_name": ds.findtext("AssemblyName") or "",
            "refseq_category": ds.findtext("RefSeq_category") or "",
            "ftp_refseq": ds.findtext("FtpPath_RefSeq") or "",
            "ftp_genbank": ds.findtext("FtpPath_GenBank") or "",
        }
        # assembly level: prefer <AssemblyStatus>, fall back to Meta assembly-status
        rec["assembly_level"] = ds.findtext("AssemblyStatus") or ""
        # genome stats live in the CDATA <Meta> blob
        meta = ds.find("Meta")
        mv = parse_meta(meta.text) if meta is not None and meta.text else {}
        rec["genome_size"] = mv.get("total_length", 0)
        rec["contig_count"] = mv.get("contig_count", 0)
        if not rec["assembly_level"] and mv.get("assembly-status"):
            rec["assembly_level"] = mv["assembly-status"]
        rec["gc_percent"] = 0.0  # computed downstream if needed
        # primary accession: RefSeq (GCF) preferred when an FTP exists, else GCA
        rec["primary_acc"] = rec["assembly_acc"]
        out.append(rec)
    return out


def fetch_docsums(uids: list[str]) -> list[dict]:
    recs: list[dict] = []
    batch = 250
    for i in range(0, len(uids), batch):
        chunk = uids[i:i + batch]
        q = urllib.parse.urlencode({"db": "assembly", "id": ",".join(chunk)})
        xml = http_get(f"{EUTILS}/esummary.fcgi?{q}")
        recs.extend(parse_docsum(xml))
        sys.stderr.write(f"\r  fetched {min(i + batch, len(uids))}/{len(uids)} docsums")
        sys.stderr.flush()
        time.sleep(0.34)
    sys.stderr.write("\n")
    return recs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-level", default=None,
                    choices=LEVEL_ORDER,
                    help="drop assemblies below this contiguity level")
    ap.add_argument("--exclude-substrings", default="",
                    help="comma extra organism-name substrings to exclude")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # NTM query: family Mycobacteriaceae (covers reclassified genera too),
    # exclude TB complex + leprae + lepromatosis.
    query = (
        '"Mycobacteriaceae"[Organism] '
        'NOT "Mycobacterium tuberculosis complex"[Organism] '
        'NOT "Mycobacterium leprae"[Organism] '
        'NOT "Mycobacterium lepromatosis"[Organism]'
    )

    t0 = time.time()
    print(f"esearch: {query}", flush=True)
    uids = esearch_uids(query)
    print(f"  {len(uids)} assembly UIDs", flush=True)

    recs = fetch_docsums(uids)

    # hard exclusion on organism name (safety net for taxonomy drift)
    excl = ["tuberculosis", "bovis", "africanum", "microti", "canettii",
            "caprae", "pinnipedii", "orygis", "suricattae", "mungi",
            "dassie", "leprae", "lepromatosis"]
    if args.exclude_substrings:
        excl += [s.strip() for s in args.exclude_substrings.split(",") if s.strip()]
    before = len(recs)
    recs = [r for r in recs if not any(x in r["organism"].lower() for x in excl)]
    print(f"  organism-name exclusion: {before} -> {len(recs)}", flush=True)

    # dedupe on primary accession (keep RefSeq over GenBank duplicates)
    seen = {}
    for r in recs:
        a = r["primary_acc"]
        if a not in seen:
            seen[a] = r
        else:
            # keep the one with a RefSeq FTP if mixed
            if r["ftp_refseq"] and not seen[a]["ftp_refseq"]:
                seen[a] = r
    recs = list(seen.values())

    # assembly-level filter
    if args.min_level:
        floor = LEVEL_ORDER.index(args.min_level)
        recs = [r for r in recs if r["assembly_level"] in LEVEL_ORDER
                and LEVEL_ORDER.index(r["assembly_level"]) >= floor]
        print(f"  after --min-level {args.min_level}: {len(recs)}", flush=True)

    recs.sort(key=lambda r: (r["species"], -r["genome_size"], r["primary_acc"]))

    # ---- write accession list (same format as 26k_ecoli_accession.txt) ----
    acc_path = outdir / "ntm_accessions.txt"
    with acc_path.open("w") as f:
        for r in recs:
            f.write(r["primary_acc"] + "\n")

    # ---- full manifest ----
    man_path = outdir / "ntm_accession_manifest.tsv"
    cols = ["primary_acc", "assembly_acc", "biosample_acc", "organism", "species",
            "genus", "assembly_level", "assembly_name", "genome_size", "contig_count",
            "gc_percent", "refseq_category", "uid", "ftp_refseq", "ftp_genbank"]
    with man_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in recs:
            w.writerow({c: r.get(c, "") for c in cols})

    # ---- species summary ----
    by_sp = defaultdict(list)
    for r in recs:
        by_sp[r["species"]].append(r)
    sp_path = outdir / "species_summary.tsv"
    sizes_all = [r["genome_size"] for r in recs if r["genome_size"]]
    with sp_path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["species", "genus", "n_assemblies", "median_genome_size",
                    "min_genome_size", "max_genome_size", "median_contigs",
                    "median_gc", "assembly_levels"])
        for sp, rs in sorted(by_sp.items(), key=lambda kv: -len(kv[1])):
            sizes = sorted(r["genome_size"] for r in rs if r["genome_size"])
            contigs = sorted(r["contig_count"] for r in rs if r["contig_count"])
            gcs = sorted(r["gc_percent"] for r in rs if r["gc_percent"])
            levels = Counter(r["assembly_level"] for r in rs)
            def med(x): return x[len(x) // 2] if x else ""
            genus = rs[0]["genus"] if rs else ""
            w.writerow([sp, genus, len(rs), med(sizes), sizes[0] if sizes else "",
                        sizes[-1] if sizes else "", med(contigs),
                        f"{med(gcs):.2f}" if gcs else "",
                        ",".join(f"{k}:{v}" for k, v in levels.most_common())])

    # ---- report ----
    level_counts = Counter(r["assembly_level"] for r in recs)
    report = {
        "query": query,
        "n_uids": len(uids),
        "n_after_exclusion_and_dedup": len(recs),
        "n_species": len(by_sp),
        "top_species": {sp: len(rs) for sp, rs in
                        sorted(by_sp.items(), key=lambda kv: -len(kv[1]))[:25]},
        "assembly_levels": dict(level_counts),
        "median_genome_size_mb": round(sizes_all[len(sizes_all) // 2] / 1e6, 2) if sizes_all else None,
        "total_bases_gb": round(sum(sizes_all) / 1e9, 2) if sizes_all else None,
        "min_level_filter": args.min_level,
        "elapsed_s": round(time.time() - t0, 1),
    }
    (outdir / "scoping_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)
    print(f"\nwrote: {acc_path}\n       {man_path}\n       {sp_path}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
