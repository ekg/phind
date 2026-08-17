#!/usr/bin/env python3
"""
panel_lib.py — shared parsing / merge / emit logic for the public
mycobacteriophage reference panel (PHIND NTM v2 external validation).

Everything here is a PURE function of the on-disk cache under
  <data_dir>/cache/{phagesdb,inphared,ncbi}/
plus <data_dir>/cache/fetch_manifest.json and taxonomy_hosts.json.
No network access, no wall-clock reads: reruns from the same cache are
byte-identical (verified by verify_panel.py).

Scope rules (documented in REPORT.md):
  * PhagesDB   — ALL sequenced Actinobacteriophages (broad host taxonomy,
                 per task spec: Mycobacterium, Gordonia, Arthrobacter,
                 Microbacterium, Streptomyces, Rhodococcus, ...).
  * INPHARED   — mycobacteriophages only: host genus in Mycobacteriaceae
                 OR description matches mycobacteriophage naming.
  * NCBI       — E-utilities nucleotide search
                 `mycobacteriophage[All Fields] OR "Mycobacterium phage"[Title]`
                 minus structural/WGS-placeholder records (explicit reasons).

Dedup:
  * alias key    = accession base (NC_001900.1 -> NC_001900): joins the same
                   submission across sources (PhagesDB Accession #, INPHARED
                   accession, NCBI accession.version).
  * sequence key = sha256 of the uppercased sequence (exact duplicates only;
                   NOT reverse-complement aware — documented caveat).
  * When the same accession yields different bytes across sources (PhagesDB
    terminal-overhang convention vs GenBank), the NCBI versioned record is
    the canonical carrier; the PhagesDB variant is kept as a crosswalk alias
    with relation=sequence_variant and its own seq_sha256.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re

# ---------------------------------------------------------------------------
# constants

PANEL_ID_PREFIX = "PHPUB"

# Genera of the family Mycobacteriaceae (2018 Gupta split kept in NCBI).
MYCOBACTERIACEAE_GENERA = {
    "Mycobacterium",
    "Mycobacteroides",
    "Mycolicibacterium",
    "Mycolicibacter",
    "Mycolicibacillus",
}
MYCOBACTERIACEAE_FAMILY_TAXID = 1762

# Host genera of sequenced Actinobacteriophages (PhagesDB scope).
ACTINOPHAGE_HOST_GENERA = MYCOBACTERIACEAE_GENERA | {
    "Gordonia",
    "Arthrobacter",
    "Microbacterium",
    "Streptomyces",
    "Rhodococcus",
    "Propionibacterium",
    "Cutibacterium",
    "Curtobacterium",
    "Corynebacterium",
    "Brevibacterium",
    "Tsukamurella",
    "Tetrasphaera",
    "Cellulomonas",
    "Agromyces",
    "Clavibacter",
    "Leifsonia",
    "Nocardia",
    "Dietzia",
    "Williamsia",
    "Segniliparus",
}

MYCO_DESC_RE = re.compile(
    r"mycobacteriophage|mycobacterium\s+phage|mycobacterium\s+virus", re.I)

# Evidence classes (cultured-vs-MAG/provirus axis).
EVIDENCE_ISOLATED = "isolated_sequenced"           # PhagesDB plaque-purified isolate
EVIDENCE_DEPOSITED = "deposited_isolate_sequence"  # GenBank/INPHARED deposit
EVIDENCE_PROPHAGE = "predicted_prophage"
EVIDENCE_METAGENOME = "metagenome_assembled"
EVIDENCE_SYNTHETIC = "synthetic_construct"

PROPHAGE_HINT_RE = re.compile(
    r"prophage|provirus|integrated\s+element|lysogen|temperate\s+element", re.I)
METAGENOME_HINT_RE = re.compile(
    r"metagenome|metavirome|uncultured|enriched\s+culture|env_sample", re.I)
SYNTHETIC_HINT_RE = re.compile(
    r"synthetic\s+construct|engineered|mutant|knockout|insertion\s+derivative", re.I)

# Host label evidence is ALWAYS a reported label — never experimental host range.
HOST_EVIDENCE_NOTE = "reported_host_label_only"
HOST_EVIDENCE_GB = "gb_host_qualifier_reported"

# Classic literature reference phages (cluster representatives from Hatfull
# et al. reviews); resolved by name against the sources, not from memory.
LITERATURE_REFERENCE_PHAGES = {
    "L5": "Hatfull & Sarkis 1993 (L5, cluster A1 reference)",
    "D29": "Jacobs-Sera et al. 2010 (D29, cluster A2 reference)",
    "TM4": "Ford et al. 1998 (TM4, cluster F reference)",
    "Bxb1": "Mediavilla et al. 2000 (Bxb1, cluster B1 reference)",
    "Bxz1": "Pedulla et al. 2003 (Bxz1, cluster C1 reference)",
    "Che9d": "Pedulla et al. 2003 (Che9d, cluster F reference)",
    "Corndog": "Pedulla et al. 2003 (Corndog, cluster E reference)",
    "Omega": "Pedulla et al. 2003 (Omega, cluster K reference)",
    "PMC": "Morris et al. 2008 (PMC, cluster G reference)",
    "Boomer": "Pope et al. 2011 (Boomer, cluster J reference)",
    "Starlet": "Hatfull et al. 2010 (Starlet, cluster A9 reference)",
    "RedRock": "Hatfull et al. 2010 (RedRock, cluster B4 reference)",
    "BP": "Sampson et al. 2009 (BP/BPs, cluster A3)",
    "BPs": "Sampson et al. 2009 (BP/BPs, cluster A3)",
    "Lily": "Hendrix 2002 review (Lily, cluster E reference)",
    "Athena": "Pope et al. 2011 (Athena, cluster L reference)",
}


# ---------------------------------------------------------------------------
# small helpers

def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def seq_sha256(seq: str) -> str:
    """sha256 over the uppercased, whitespace-stripped sequence."""
    return sha256_bytes("".join(seq.split()).upper().encode("ascii"))


def accession_base(acc: str) -> str:
    """NC_001900.1 -> NC_001900 ; tolerant of stray whitespace."""
    return re.sub(r"\.\d+$", "", acc.strip())


def is_versioned_accession(acc: str) -> bool:
    return re.fullmatch(r"[A-Za-z]{1,6}_?\d+\.\d+", acc.strip()) is not None


def host_genus_of(host_label: str) -> str:
    """First token of a host label, normalized; '' when empty."""
    t = (host_label or "").strip().split()
    return t[0] if t else ""


def read_fetch_manifest(cache_dir: str) -> dict:
    return json.load(open(os.path.join(cache_dir, "fetch_manifest.json")))


# ---------------------------------------------------------------------------
# parsers

def parse_phagesdb_tsv(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_phagesdb_fasta(path: str) -> tuple[dict[str, dict], list[dict]]:
    """Parse the bulk Actinobacteriophages-All.fasta.

    Returns (records, extras):
      records — {phage_name: {header, seq, header_host_genus}} first occurrence
      extras  — subsequent/malformed headers as dicts with name, header, seq,
                genus, issue ('duplicate_name' | 'malformed_header')
    Header conventions: '<HostGenus> phage <Name> ...', '<Name> ...',
    'Mycobacteriophage <Name> ...', '<HostGenus> <Name> ...'.
    """
    records: dict[str, dict] = {}
    extras: list[dict] = []
    genus = ""
    name = header = None
    chunks: list[str] = []

    def _extract(hdr: str):
        toks = hdr.lstrip(">").strip().split()
        g, idx = "", 0
        if toks and toks[0] in ACTINOPHAGE_HOST_GENERA:
            g = toks[0]
            idx = 1
        elif toks and toks[0].lower() in ("mycobacteriophage", "mycobacteriophages"):
            g = "Mycobacterium"
            idx = 1
        if idx < len(toks) and toks[idx].lower() in ("phage", "phages"):
            idx += 1
        if idx >= len(toks) or toks[idx].lower() in (
                "complete", "nearly", "complete,", "whole"):
            return g, None, " ".join(toks)
        return g, toks[idx], " ".join(toks)

    pending = False

    def _flush(nm, hd, gn, sq):
        nonlocal pending
        if not pending:
            return
        pending = False
        if nm is None:
            extras.append({"name": "", "header": hd, "seq": sq,
                           "genus": gn, "issue": "malformed_header"})
        elif nm in records:
            extras.append({"name": nm, "header": hd, "seq": sq,
                           "genus": gn, "issue": "duplicate_name"})
        else:
            records[nm] = {"header": hd, "seq": sq, "header_host_genus": gn}

    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                _flush(name, header, genus, "".join(chunks))
                genus, name, header = _extract(line)
                pending = True
                chunks = []
            else:
                chunks.append(line.strip())
    _flush(name, header, genus, "".join(chunks))
    return records, extras


def parse_inphared_table(path: str) -> list[dict]:
    """INPHARED millardlab_website_table (gzipped, HTML-linked accessions)."""
    rows = []
    tag_re = re.compile(r"<[^>]+>")
    with gzip.open(path, "rt", newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append({
                "accession": tag_re.sub("", r.get("Accession", "")).strip(),
                "description": (r.get("Description") or "").strip(),
                "length_kb": (r.get("Genome Length (KB)") or "").strip(),
                "gc": (r.get("molGC (%)") or "").strip(),
                "genus": (r.get("Genus") or "").strip(),
                "subfamily": (r.get("Sub-family") or "").strip(),
                "family": (r.get("Family") or "").strip(),
                "host": (r.get("Host") or "").strip(),
            })
    return rows


def inphared_scope(row: dict) -> bool:
    return (host_genus_of(row["host"]) in MYCOBACTERIACEAE_GENERA
            or bool(MYCO_DESC_RE.search(row["description"])))


def load_inphared_seqs(path: str, wanted: set[str]) -> dict[str, dict]:
    """Stream the INPHARED fasta once; keep only wanted accessions.

    Returns {accession_base: {desc, seq}} — keyed by base accession so
    versioned FASTA ids (PV876933.1) match unversioned table accessions
    (PV876933) and vice versa.
    """
    wanted_bases = {accession_base(a) for a in wanted}
    out: dict[str, dict] = {}
    for fid, desc, seq in iter_fasta_gz(path):
        base = accession_base(fid)
        if fid in wanted or base in wanted_bases:
            out[base] = {"desc": desc, "seq": seq, "fasta_id": fid}
    return out


def iter_fasta_gz(path: str):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        fid = None
        desc = ""
        chunks: list[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if fid is not None:
                    yield fid, desc, "".join(chunks)
                parts = line[1:].split(None, 1)
                fid = parts[0] if parts else ""
                desc = parts[1].strip() if len(parts) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip())
        if fid is not None:
            yield fid, desc, "".join(chunks)


def parse_ncbi_esummary(paths: list[str]) -> dict[int, dict]:
    out: dict[int, dict] = {}
    for p in paths:
        d = json.load(open(p))["result"]
        for u in d["uids"]:
            r = d[u]
            out[int(u)] = {
                "uid": int(u),
                "accessionversion": r.get("accessionversion", ""),
                "caption": r.get("caption", ""),
                "title": r.get("title", ""),
                "slen": int(r.get("slen", 0) or 0),
                "organism": r.get("organism", ""),
                "taxid": r.get("taxid", ""),
                "sourcedb": r.get("sourcedb", ""),
            }
    return out


def parse_ncbi_fasta(paths: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for p in paths:
        fid = header = None
        chunks: list[str] = []
        with open(p, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if fid is not None:
                        out[fid] = {"header": header, "seq": "".join(chunks)}
                    header = line[1:].strip()
                    fid = header.split()[0] if header.split() else ""
                    chunks = []
                else:
                    chunks.append(line.strip())
        if fid is not None:
            out[fid] = {"header": header, "seq": "".join(chunks)}
    return out


PDB_ID_RE = re.compile(r"^pdb\|(\w+)\|(\w+)$", re.I)


def ncbi_uid_to_fasta_key(uid_rec: dict) -> str:
    av = uid_rec["accessionversion"]
    m = PDB_ID_RE.match(av)
    if m:
        return f"pdb|{m.group(1)}|{m.group(2)}"
    return av


def parse_gb_remainder(path: str) -> dict[str, dict]:
    """Parse GenBank flat files fetched for NCBI-only records lacking host.

    Returns accession_base -> {host_qualifier, topology, mol_type,
    completeness_note}. Only the FIRST record for each accession is kept.
    """
    out: dict[str, dict] = {}
    if not os.path.exists(path):
        return out
    loc_re = re.compile(r"^LOCUS\s+(\S+)\s+\d+\s+bp\s+(\S+)\s+(\S+)\s+(\S+)", re.M)
    acc_re = re.compile(r"^ACCESSION\s+(\S+)", re.M)
    ver_re = re.compile(r"^VERSION\s+(\S+)", re.M)
    host_re = re.compile(r'/(?:lab_)?host="([^"]+)"')
    comp_re = re.compile(r'COMPLETENESS[:\s]+([^"\n]+)', re.I)
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    for rec in text.split("//\n"):
        rec = rec.strip()
        if not rec or "LOCUS" not in rec[:200]:
            continue
        mv = ver_re.search(rec)
        ma = acc_re.search(rec)
        acc = mv.group(1) if mv else (ma.group(1) if ma else "")
        if not acc:
            continue
        ml = loc_re.search(rec)
        mh = host_re.search(rec)
        mc = comp_re.search(rec)
        out.setdefault(accession_base(acc), {
            "accession": acc,
            "host_qualifier": mh.group(1).strip() if mh else "",
            "topology": ml.group(3) if ml else "",
            "mol_type": ml.group(2) if ml else "",
            "completeness_note": mc.group(1).strip() if mc else "",
        })
    return out


# ---------------------------------------------------------------------------
# merge

def merge_records(pdb_meta, pdb_fasta, pdb_extras, inp_rows, inp_seqs,
                  ncbi_summ, ncbi_fa, taxonomy, gb_remainder, fetch_manifest) -> dict:
    """Pure merge of parsed sources -> panel structures.

    Returns {"records", "crosswalk", "exclusions", "source_counts",
             "needs_gb_remainder", "inphared_out_of_scope"}.
    """
    fm = fetch_manifest
    ts = {k: v.get("retrieved_at", "") for k, v in fm.items()}
    urls = {k: v.get("url", "") for k, v in fm.items()}

    exclusions: list[dict] = []

    def excl(source, ident, reason, detail=""):
        exclusions.append({"source": source, "identifier": ident,
                           "reason": reason, "detail": str(detail)[:200]})

    # ---------------- PhagesDB leg --------------------------------------
    pdb_by_name = {r["Phage Name"]: r for r in pdb_meta}
    phagesdb: dict[str, dict] = {}
    for x in pdb_extras:
        if x["issue"] == "malformed_header":
            excl("phagesdb_fasta", x["header"][:60],
                 "malformed_header_missing_name",
                 f"host_genus={x['genus']!r}; len={len(x['seq'])}")
            continue
        assert x["issue"] == "duplicate_name"
        first = pdb_fasta.get(x["name"])
        if first is not None and seq_sha256(first["seq"]) == seq_sha256(x["seq"]):
            excl("phagesdb_fasta", x["name"], "exact_duplicate_in_bulk_fasta",
                 "second identical record in bulk FASTA")
        else:
            nm = x["name"] + "#2"  # distinct variant of a duplicated name
            phagesdb[nm] = {"name": nm, "meta": pdb_by_name.get(x["name"]),
                            "seq": x["seq"], "header": x["header"],
                            "host_label": "", "host_genus": x["genus"],
                            "accession": (pdb_by_name.get(x["name"]) or {}).get(
                                "Accession #", "").strip(),
                            "metadata_source": "phagesdb_bulk_fasta_variant"}
            excl("phagesdb_fasta", x["name"],
                 "retained_as_suffixed_variant",
                 f"second record with different sequence kept as {nm}")
    for name, fa in sorted(pdb_fasta.items()):
        meta = pdb_by_name.get(name)
        host_label = (meta or {}).get("Host", "") or ""
        hgenus = host_genus_of(host_label) or fa.get("header_host_genus", "")
        if meta is None and hgenus not in ACTINOPHAGE_HOST_GENERA:
            excl("phagesdb_fasta", name, "no_metadata_row_and_host_genus_unresolvable",
                 f"header_host_genus={hgenus!r}")
            continue
        acc = ((meta or {}).get("Accession #") or "").strip()
        if acc.lower() in ("none", "n/a", "-", ""):
            acc = ""
        phagesdb[name] = {
            "name": name, "meta": meta, "seq": fa["seq"], "header": fa["header"],
            "host_label": host_label, "host_genus": hgenus, "accession": acc,
            "metadata_source": "phagesdb_metadata" if meta else "fasta_header_only",
        }
    for name, meta in sorted(pdb_by_name.items()):
        if name not in phagesdb:
            excl("phagesdb_metadata", name, "no_sequence_in_bulk_fasta",
                 f"in_genbank={meta.get('In GenBank?')}; "
                 f"accession={meta.get('Accession #')}")

    # ---------------- INPHARED leg ---------------------------------------
    inphared: dict[str, dict] = {}
    out_of_scope = 0
    for row in inp_rows:
        acc = row["accession"]
        if not acc:
            excl("inphared_table", row["description"][:60], "empty_accession")
            continue
        if not inphared_scope(row):
            out_of_scope += 1
            continue
        seq = inp_seqs.get(accession_base(acc))
        if seq is None:
            excl("inphared", acc, "no_sequence_in_release_fasta",
                 f"host={row['host']}; desc={row['description'][:60]}")
            seq = {"desc": "", "seq": ""}
        inphared[accession_base(acc)] = {"row": row, "seq": seq["seq"]}

    # ---------------- NCBI leg -------------------------------------------
    ncbi: dict[str, dict] = {}
    for uid, r in ncbi_summ.items():
        fkey = ncbi_uid_to_fasta_key(r)
        base = accession_base(r["accessionversion"])
        title = r["title"]
        if fkey.startswith("pdb|"):
            excl("ncbi", r["accessionversion"], "pdb_structure_record_not_nucleotide",
                 title[:100])
            continue
        fa = ncbi_fa.get(fkey)
        if fa is None:
            if r["slen"] <= 10 or "shotgun sequencing project" in title:
                excl("ncbi", r["accessionversion"], "wgs_master_or_placeholder_record",
                     f"slen={r['slen']}; {title[:100]}")
            else:
                excl("ncbi", r["accessionversion"], "no_fasta_from_efetch",
                     f"uid={uid}; slen={r['slen']}; {title[:100]}")
            continue
        if len("".join(fa["seq"].split())) == 0:
            excl("ncbi", r["accessionversion"], "empty_sequence")
            continue
        ncbi[base] = {"summ": r, "seq": fa["seq"]}

    # ---------------- union by accession base -----------------------------
    entries: list[dict] = []
    base2i: dict[str, int] = {}

    def new_entry() -> int:
        entries.append({"ncbi": None, "inphared": None, "phagesdb": None})
        return len(entries) - 1

    def attach(base: str, which: str, payload, source: str, ident: str = ""):
        if base in base2i:
            i = base2i[base]
        else:
            i = new_entry()
            base2i[base] = i
        if entries[i][which] is None:
            entries[i][which] = payload
        else:  # duplicate identifier within one source
            excl(source, ident or base, "duplicate_identifier_within_source",
                 f"first occurrence kept; colliding key={base}")
            return None
        return i

    for base in sorted(ncbi):
        attach(base, "ncbi", ncbi[base], "ncbi")
    for base in sorted(inphared):
        attach(base, "inphared", inphared[base], "inphared")
    for name in sorted(phagesdb):
        e = phagesdb[name]
        base = accession_base(e["accession"]) if e["accession"] else f"phagesdb_name:{name}"
        if base.startswith("phagesdb_name:"):
            if base in base2i:
                excl("phagesdb", name, "duplicate_phage_name_alias",
                     f"collides with {base}")
                continue
            base2i[base] = new_entry()
            entries[base2i[base]]["phagesdb"] = e
        else:
            attach(base, "phagesdb", e, "phagesdb", ident=name)

    # ---------------- curate NCBI-only records (text-search artifacts) ----
    # Records found ONLY by the NCBI esearch term must still look like
    # actinobacteriophage genomes; linked PhagesDB/INPHARED records are
    # trusted without this screen.
    ncbi_only_screened = 0
    for i, e in enumerate(entries):
        if e["ncbi"] is None or (e["phagesdb"] or e["inphared"]):
            continue
        av = e["ncbi"]["summ"]["accessionversion"]
        title = e["ncbi"]["summ"]["title"]
        seq = "".join(e["ncbi"]["seq"].split())
        gbr = gb_remainder.get(accession_base(av), {})
        host_gen = host_genus_of(gbr.get("host_qualifier", ""))
        if not re.search(r"phage|virus", title, re.I):
            excl("ncbi", av, "non_phage_record_matched_search",
                 f"title={title[:100]}")
            e["ncbi"] = None
        elif len(seq) < 10_000:
            excl("ncbi", av, "short_fragment_not_genome_scale",
                 f"len={len(seq)}; title={title[:80]}")
            e["ncbi"] = None
        elif not (MYCO_DESC_RE.search(title)
                  or host_gen in ACTINOPHAGE_HOST_GENERA
                  or any(g in title for g in ACTINOPHAGE_HOST_GENERA)):
            excl("ncbi", av, "no_actinobacteriophage_evidence",
                 f"title={title[:100]}; gb_host={gbr.get('host_qualifier','')}")
            e["ncbi"] = None
        else:
            ncbi_only_screened += 1
    # drop emptied entries and rebuild index mapping
    survivors_pre = [i for i, e in enumerate(entries)
                     if e["ncbi"] or e["inphared"] or e["phagesdb"]]
    old2new = {old: new for new, old in enumerate(survivors_pre)}
    entries = [entries[i] for i in survivors_pre]
    base2i = {b: old2new[i] for b, i in base2i.items() if i in old2new}

    # ---------------- exact-sequence dedup --------------------------------
    def carrier(e):
        if e["ncbi"] and e["ncbi"]["seq"]:
            return e["ncbi"]["seq"]
        if e["inphared"] and e["inphared"]["seq"]:
            return e["inphared"]["seq"]
        if e["phagesdb"] and e["phagesdb"]["seq"]:
            return e["phagesdb"]["seq"]
        return ""

    by_hash: dict[str, list[int]] = {}
    for i, e in enumerate(entries):
        c = carrier(e)
        if c:
            by_hash.setdefault(seq_sha256(c), []).append(i)

    merged_into: dict[int, int] = {}
    for h, idxs in sorted(by_hash.items()):
        if len(idxs) > 1:
            # prefer keeping the entry with an NCBI versioned record
            keep = min(idxs, key=lambda i: (entries[i]["ncbi"] is None, i))
            for j in idxs:
                if j != keep:
                    merged_into[j] = keep

    def absorb(i, j):
        """merge entry j's metadata into entry i (sequences identical)."""
        for which in ("ncbi", "inphared", "phagesdb"):
            if entries[i][which] is None:
                entries[i][which] = entries[j][which]

    for j in sorted(merged_into):
        tgt = merged_into[j]
        while tgt in merged_into:
            tgt = merged_into[tgt]
        merged_into[j] = tgt
        absorb(tgt, j)

    survivors = [i for i in range(len(entries)) if i not in merged_into]

    # ---------------- assemble + panel ids --------------------------------
    records = [_assemble(entries[i], ts, urls, taxonomy, gb_remainder)
               for i in survivors]
    assign_panel_ids(records)
    # map entry index -> panel id for crosswalk
    idx2pid = {i: rec["panel_id"] for i, rec in zip(survivors, records)}
    idx2rec = {i: rec for i, rec in zip(survivors, records)}

    # ---------------- crosswalk -------------------------------------------
    crosswalk: list[dict] = []

    def cw(pid, alias, source, relation, sha="", note=""):
        crosswalk.append({"panel_id": pid, "alias": alias, "source": source,
                          "relation": relation, "alias_seq_sha256": sha,
                          "note": note})

    for i in survivors:
        pid = idx2pid[i]
        rec = idx2rec[i]
        e = entries[i]
        if e["ncbi"]:
            pid_key = e["ncbi"]["summ"]["accessionversion"]
            base = accession_base(pid_key)
            if base != pid_key:
                cw(pid, base, "ncbi", "accession_base_alias")
            if e["inphared"] and e["inphared"]["row"]["accession"] != pid_key:
                cw(pid, e["inphared"]["row"]["accession"], "inphared",
                   "accession_alias")
        elif e["inphared"]:
            pid_key = e["inphared"]["row"]["accession"]
            base = accession_base(pid_key)
            if base != pid_key:
                cw(pid, base, "inphared", "accession_base_alias")
        if e["phagesdb"]:
            e2 = e["phagesdb"]
            cw(pid, f"phagesdb_name:{e2['name']}", "phagesdb", "name_alias")
            if e2["accession"] and not (
                    e["ncbi"] and accession_base(e2["accession"]) ==
                    accession_base(e["ncbi"]["summ"]["accessionversion"])):
                cw(pid, e2["accession"], "phagesdb", "accession_alias")
            former = ((e2["meta"] or {}).get("Former Name(s)") or "").strip()
            if former and former.lower() not in ("none", "n/a", "-"):
                cw(pid, f"phagesdb_former_name:{former}", "phagesdb",
                   "former_name_alias")
        # PhagesDB sequence variant (terminal overhang) vs NCBI carrier
        if e["ncbi"] and e["phagesdb"] and e["phagesdb"]["seq"]:
            sh = seq_sha256(e["phagesdb"]["seq"])
            if sh != rec["seq_sha256"]:
                cw(pid, f"phagesdb_seq:{e['phagesdb']['name']}", "phagesdb",
                   "sequence_variant", sh,
                   f"phagesdb len={len(''.join(e['phagesdb']['seq'].split()))} "
                   f"differs from carrier "
                   f"len={rec['length_bp']} "
                   "(terminal-overhang convention)")

    for j in sorted(merged_into):
        tgt = merged_into[j]
        pid = idx2pid[tgt]
        sha = seq_sha256(carrier(entries[j])) if carrier(entries[j]) else ""
        if entries[j]["ncbi"]:
            cw(pid, entries[j]["ncbi"]["summ"]["accessionversion"], "ncbi",
               "same_sequence", sha, "exact sequence duplicate merged")
        if entries[j]["inphared"]:
            cw(pid, entries[j]["inphared"]["row"]["accession"], "inphared",
               "same_sequence", sha, "exact sequence duplicate merged")
        if entries[j]["phagesdb"]:
            cw(pid, f"phagesdb_name:{entries[j]['phagesdb']['name']}",
               "phagesdb", "same_sequence", sha, "exact sequence duplicate merged")

    # ---------------- needs_gb_remainder -----------------------------------
    needs_gb = sorted(
        accession_base(e["ncbi"]["summ"]["accessionversion"])
        for i in survivors
        if (e := entries[i])["ncbi"] is not None
        and e["phagesdb"] is None and e["inphared"] is None
        and accession_base(e["ncbi"]["summ"]["accessionversion"]) not in gb_remainder
    )

    counts = {
        "phagesdb": {"metadata_rows": len(pdb_meta),
                     "fasta_records": len(pdb_fasta) + len(pdb_extras),
                     "retained": len(phagesdb)},
        "inphared": {"table_rows": len(inp_rows),
                     "in_scope_rows": len(inphared),
                     "out_of_scope": out_of_scope,
                     "with_sequence": sum(1 for v in inphared.values() if v["seq"])},
        "ncbi": {"uids": len(ncbi_summ), "with_fasta": len(ncbi),
                 "ncbi_only_retained": ncbi_only_screened},
        "entries": {"total": len(entries),
                    "merged_duplicates": len(merged_into),
                    "survivors": len(survivors)},
    }
    return {"records": records, "crosswalk": crosswalk, "exclusions": exclusions,
            "source_counts": counts, "needs_gb_remainder": needs_gb}


def _assemble(e, ts, urls, taxonomy, gb_remainder) -> dict:
    nc, inp, pdb = e["ncbi"], e["inphared"], e["phagesdb"]
    seq = _carrier_seq(e)
    seq_clean = "".join(seq.split())

    if nc:
        primary_acc = nc["summ"]["accessionversion"]
        primary_src, seq_artifact = "ncbi", "ncbi_efetch_fasta"
        retrieval, surl = ts.get("ncbi_efetch", ""), urls.get("ncbi_efetch", "")
    elif inp:
        primary_acc = inp["row"]["accession"]
        primary_src, seq_artifact = "inphared", "inphared_genomes_fasta"
        retrieval, surl = ts.get("inphared_genomes_fa", ""), urls.get("inphared_table", "")
    else:
        primary_acc = pdb["accession"] or pdb["name"]
        primary_src, seq_artifact = "phagesdb", "phagesdb_bulk_fasta"
        retrieval, surl = ts.get("phagesdb_bulk_fasta", ""), urls.get("phagesdb_bulk_fasta", "")

    # ---- host (reported labels only; never equated with verified range) --
    host_label, host_source = "", ""
    host_evidence = HOST_EVIDENCE_NOTE
    if pdb and pdb["host_label"]:
        host_label, host_source = pdb["host_label"], "phagesdb_metadata"
    elif pdb and not pdb["meta"] and pdb["host_genus"]:
        host_label, host_source = pdb["host_genus"], "phagesdb_fasta_header"
    elif inp and inp["row"]["host"]:
        host_label, host_source = inp["row"]["host"], "inphared_table"
    else:
        base = accession_base(primary_acc)
        gbr = gb_remainder.get(base)
        if gbr and gbr.get("host_qualifier"):
            host_label = gbr["host_qualifier"]
            host_source = "ncbi_genbank_host_qualifier"
            host_evidence = HOST_EVIDENCE_GB

    hgenus = host_genus_of(host_label)
    tax = taxonomy.get(hgenus.lower(), {}) if hgenus else {}
    host_taxid = str(tax.get("taxid", "")) if tax else ""
    host_taxid_source = "ncbi_taxonomy" if host_taxid else ""
    in_myco = hgenus in MYCOBACTERIACEAE_GENERA or (
        tax.get("family_taxid") == MYCOBACTERIACEAE_FAMILY_TAXID)
    title_desc = " ".join(filter(None, [
        inp["row"]["description"] if inp else "",
        nc["summ"]["title"] if nc else "",
    ]))
    if in_myco:
        scope = "mycobacteriophage"
    elif MYCO_DESC_RE.search(title_desc):
        scope = "mycobacteriophage"
    elif hgenus in ACTINOPHAGE_HOST_GENERA:
        scope = "other_actinobacteriophage"
    elif not host_label:
        scope = "phage_host_unresolved"
    else:
        scope = "other_host_phage"

    # ---- evidence class ---------------------------------------------------
    desc = title_desc
    if pdb:
        ev = EVIDENCE_ISOLATED
    elif SYNTHETIC_HINT_RE.search(desc):
        ev = EVIDENCE_SYNTHETIC
    elif PROPHAGE_HINT_RE.search(desc):
        ev = EVIDENCE_PROPHAGE
    elif METAGENOME_HINT_RE.search(desc):
        ev = EVIDENCE_METAGENOME
    else:
        ev = EVIDENCE_DEPOSITED

    # ---- topology / completeness ------------------------------------------
    end_type = (pdb["meta"].get("End Type", "") if pdb and pdb["meta"] else "") or ""
    topology, topology_source = "", ""
    if pdb and end_type and end_type.lower() not in ("unknown", "none", "n/a", "-"):
        topology, topology_source = "linear", f"phagesdb_end_type:{end_type}"
    base = accession_base(primary_acc)
    gbr = gb_remainder.get(base, {})
    if not topology and gbr.get("topology"):
        topology, topology_source = gbr["topology"], "genbank_locus"

    completeness, completeness_source = "", ""
    if pdb and pdb["meta"] and pdb["meta"].get("Finished Sequence?") == "True":
        completeness, completeness_source = "finished_sequence", "phagesdb_metadata"
    if not completeness and nc:
        t = nc["summ"]["title"].lower()
        if "complete genome" in t or "complete sequence" in t:
            completeness, completeness_source = "complete_genome_record", "ncbi_title"
    if not completeness and inp:
        completeness, completeness_source = "length_reported_only", "inphared_table"
    if gbr.get("completeness_note"):
        completeness = gbr["completeness_note"]
        completeness_source = "genbank_completeness_note"

    # ---- literature reference ----------------------------------------------
    lit = ""
    if pdb:
        lit = LITERATURE_REFERENCE_PHAGES.get(pdb["name"], "")
    if not lit and nc:
        m = re.match(r"Mycobacterium phage (\S+)", nc["summ"]["title"])
        if m:
            lit = LITERATURE_REFERENCE_PHAGES.get(m.group(1), "")

    sources = ";".join(s for s, v in (("ncbi", nc), ("inphared", inp),
                                      ("phagesdb", pdb)) if v)

    m = pdb["meta"] if (pdb and pdb["meta"]) else {}
    return {
        "panel_id": "",
        "primary_accession": primary_acc,
        "primary_source": primary_src,
        "phage_name": pdb["name"] if pdb else "",
        "scope": scope,
        "in_mycobacteriaceae_host": str(bool(in_myco)).lower(),
        "length_bp": str(len(seq_clean)),
        "seq_sha256": seq_sha256(seq),
        "host_reported": host_label,
        "host_reported_source": host_source,
        "host_evidence": host_evidence,
        "host_genus": hgenus,
        "host_taxid": host_taxid,
        "host_taxid_source": host_taxid_source,
        "evidence_class": ev,
        "topology": topology,
        "topology_source": topology_source,
        "phagesdb_end_type": end_type,
        "completeness": completeness,
        "completeness_source": completeness_source,
        "phagesdb_cluster": m.get("Cluster", ""),
        "phagesdb_subcluster": m.get("Subcluster", ""),
        "phagesdb_found_country": m.get("Found Country", ""),
        "phagesdb_found_year": m.get("Year Found", ""),
        "inphared_accession": inp["row"]["accession"] if inp else "",
        "inphared_host": inp["row"]["host"] if inp else "",
        "inphared_family": inp["row"]["family"] if inp else "",
        "ncbi_title": nc["summ"]["title"] if nc else "",
        "ncbi_taxid": str(nc["summ"]["taxid"]) if nc else "",
        "ncbi_sourcedb": nc["summ"]["sourcedb"] if nc else "",
        "gb_host_qualifier": gbr.get("host_qualifier", ""),
        "sources": sources,
        "seq_source": seq_artifact,
        "retrieval_timestamp": retrieval,
        "source_url": surl,
        "literature_reference": lit,
        "_seq": seq,
    }


def _carrier_seq(e: dict) -> str:
    if e.get("ncbi") and e["ncbi"].get("seq"):
        return e["ncbi"]["seq"]
    if e.get("inphared") and e["inphared"].get("seq"):
        return e["inphared"]["seq"]
    if e.get("phagesdb") and e["phagesdb"].get("seq"):
        return e["phagesdb"]["seq"]
    return ""


def assign_panel_ids(records: list[dict]) -> None:
    order = sorted(
        range(len(records)),
        key=lambda i: (
            0 if records[i]["scope"] == "mycobacteriophage" else 1,
            records[i]["host_genus"],
            records[i]["phage_name"] or accession_base(records[i]["primary_accession"]),
            records[i]["primary_accession"],
        ),
    )
    for n, idx in enumerate(order, start=1):
        records[idx]["panel_id"] = f"{PANEL_ID_PREFIX}-{n:06d}"


# ---------------------------------------------------------------------------
# deterministic emit

def ascii_safe(s: str) -> str:
    """NFKD-normalize + drop non-ASCII (host labels like mc²155 -> mc2155)."""
    import unicodedata
    return unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")


def fasta_bytes(rec: dict) -> bytes:
    label = rec["phage_name"] or accession_base(rec["primary_accession"])
    hdr = (f"{rec['panel_id']} {label}|{rec['primary_accession']}|"
           f"host={rec['host_reported']}")
    seq = "".join(rec["_seq"].split()).upper()
    lines = [">" + ascii_safe(hdr)]
    for k in range(0, len(seq), 60):
        lines.append(seq[k:k + 60])
    return ("\n".join(lines) + "\n").encode("ascii")


MANIFEST_COLUMNS = [
    "panel_id", "primary_accession", "primary_source", "phage_name", "scope",
    "in_mycobacteriaceae_host", "length_bp", "seq_sha256", "fasta_bytes_sha256",
    "host_reported", "host_reported_source", "host_evidence", "host_genus",
    "host_taxid", "host_taxid_source", "evidence_class", "topology",
    "topology_source", "phagesdb_end_type", "completeness",
    "completeness_source", "phagesdb_cluster", "phagesdb_subcluster",
    "phagesdb_found_country", "phagesdb_found_year", "inphared_accession",
    "inphared_host", "inphared_family", "ncbi_title", "ncbi_taxid",
    "ncbi_sourcedb", "gb_host_qualifier", "sources", "seq_source",
    "retrieval_timestamp", "source_url", "literature_reference",
]

CROSSWALK_COLUMNS = ["panel_id", "alias", "source", "relation",
                     "alias_seq_sha256", "note"]
EXCLUSION_COLUMNS = ["source", "identifier", "reason", "detail"]


def manifest_rows(records: list[dict]) -> list[dict]:
    rows = []
    for rec in sorted(records, key=lambda r: r["panel_id"]):
        row = dict(rec)
        row["fasta_bytes_sha256"] = sha256_bytes(fasta_bytes(rec))
        rows.append(row)
    return rows


def write_tsv(path: str, columns: list[str], rows: list[dict]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_panel_fasta(path: str, records: list[dict]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        for rec in sorted(records, key=lambda r: r["panel_id"]):
            fh.write(fasta_bytes(rec))
