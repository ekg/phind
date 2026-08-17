#!/usr/bin/env python3
"""
test_panel.py — unit tests for the public mycobacteriophage panel builder,
using tiny synthetic source fixtures (no network, no real data).

Run:  python3 -m pytest ntm/v2/external_validation/scripts/test_panel.py -v
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import panel_lib as pl  # noqa: E402
import build_panel as bp  # noqa: E402

FM = {
    "phagesdb_metadata": {"url": "https://phagesdb.org/data/?set=seq&type=full",
                          "retrieved_at": "2026-08-17T16:48:00Z"},
    "phagesdb_bulk_fasta": {"url": "https://phagesdb.org/media/Actinobacteriophages-All.fasta",
                            "retrieved_at": "2026-08-17T16:49:00Z"},
    "inphared_table": {"url": "https://s3.climb.ac.uk/...table.txt.gz",
                       "retrieved_at": "2026-08-17T16:48:10Z"},
    "inphared_genomes_fa": {"url": "https://s3.climb.ac.uk/...genomes.fa.gz",
                            "retrieved_at": "2026-08-17T16:49:20Z"},
    "ncbi_esearch": {"url": "https://eutils.ncbi.nlm.nih.gov/...",
                     "retrieved_at": "2026-08-17T16:44:00Z"},
    "ncbi_efetch": {"url": "https://eutils.ncbi.nlm.nih.gov/...efetch",
                    "retrieved_at": "2026-08-17T16:46:00Z"},
    "ncbi_esummary": {"url": "https://eutils.ncbi.nlm.nih.gov/...esummary",
                      "retrieved_at": "2026-08-17T16:45:00Z"},
    "ncbi_gb_remainder": {"url": "https://eutils...efetch gb",
                          "retrieved_at": "2026-08-17T17:08:00Z"},
    "ncbi_backfill": {"url": "https://eutils...by id",
                      "retrieved_at": "2026-08-17T17:15:00Z"},
}

TAXONOMY = {
    "mycobacterium": {"taxid": "1763", "rank": "genus",
                      "family_taxid": "1762", "lineage": ["Mycobacteriaceae"]},
    "gordonia": {"taxid": "79255", "rank": "genus", "family_taxid": "", "lineage": []},
}


# ---------------------------------------------------------------------------
# helpers

def _pdb_meta(rows):
    cols = ["Phage Name", "Host", "Cluster", "Subcluster", "Finished Sequence?",
            "Temperate?", "Archive Titer", "In GenBank?", "Accession #",
            "Published in a Paper?", "Genome Length(bp)", "GC%", "End Type",
            "Overhang Length (bp)", "Overhang Sequence", "Term Rpt Length",
            "# ORFs", "# tRNAs", "# tmRNAs", "Finding Institution", "Program",
            "Finder Name(s)", "Found City", "Found State", "Found Country",
            "Found GPS Lat", "Found GPS Long", "Year Found", "Morphotype",
            "Seq Facility", "Shotgun Seq Method", "From Enriched Sample?",
            "Former Name(s)", "Date Finished", "Phamerated?",
            "Annotation Status", "Submitted DNAM File?", "Submitted Minimal File?",
            "Author List?", "Cover Sheet?", "Annotation Submission Date"]
    out = ["\t".join(cols)]
    for r in rows:
        d = {c: "" for c in cols}
        d.update(r)
        out.append("\t".join(d[c] for c in cols))
    return "\n".join(out) + "\n"


def _fasta(recs, prefix=">"):
    return "".join(f"{prefix}{h}\n{s}\n" for h, s in recs)


def _inphared_table(rows):
    cols = ["Accession", "Description", "Genome Length (KB)", "molGC (%)",
            "Genus", "Sub-family", "Family", "Host"]
    out = ["\t".join(cols)]
    for r in rows:
        d = {c: "" for c in cols}
        d.update(r)
        d["Accession"] = f'<a href="https://www.ncbi.nlm.nih.gov/nuccore/{d["Accession"]}">{d["Accession"]}</a>'
        out.append("\t".join(d[c] for c in cols))
    return "\n".join(out) + "\n"


def _summ(uid, acc_ver, title, slen, organism="Mycobacterium phage T", taxid="28369"):
    return {"uid": uid, "accessionversion": acc_ver, "caption": acc_ver.split(".")[0],
            "title": title, "slen": slen, "organism": organism, "taxid": taxid,
            "sourcedb": "GenBank"}


S1 = "ACGTACGTAC" * 30        # 300 bp phage A
S2 = "TTTTGGGGCC" * 40        # 400 bp phage B
S3 = "ACGTACGTAC" * 30        # identical to S1 (RefSeq twin of phage A)


@pytest.fixture()
def cache(tmp_path):
    d = tmp_path / "cache"
    for sub in ("phagesdb", "inphared", "ncbi"):
        (d / sub).mkdir(parents=True)

    (d / "phagesdb" / "phagesdb_sequenced_phages_full.tsv").write_text(_pdb_meta([
        {"Phage Name": "Alpha", "Host": "Mycobacterium smegmatis mc²155",
         "Cluster": "A", "Subcluster": "A1", "Finished Sequence?": "True",
         "In GenBank?": "True", "Accession #": "JN000001",
         "Genome Length(bp)": "300", "End Type": "3' Sticky Overhang",
         "Found Country": "USA", "Year Found": "2010", "Former Name(s)": "OldAlpha"},
        {"Phage Name": "Beta", "Host": "Gordonia terrae 3612",
         "Cluster": "B", "Finished Sequence?": "True", "In GenBank?": "True",
         "Accession #": "JN000002", "Genome Length(bp)": "400",
         "End Type": "defined ends"},
        {"Phage Name": "Gamma", "Host": "Microbacterium foliorum NRRL B-24224",
         "Cluster": "EA", "Finished Sequence?": "True", "In GenBank?": "True",
         "Accession #": "JN000003", "Genome Length(bp)": "320",
         "End Type": "circularly permuted"},
        {"Phage Name": "NoSeq", "Host": "Mycobacterium smegmatis mc²155",
         "In GenBank?": "False", "Accession #": ""},
    ]))
    (d / "phagesdb" / "Actinobacteriophages-All.fasta").write_text(_fasta([
        ("Mycobacterium phage Alpha complete sequence, 300 bp, Cluster A1", S1),
        # repeated name, identical sequence -> exact duplicate exclusion
        ("Mycobacterium phage Alpha complete sequence, 300 bp", S1),
        ("Gordonia phage Beta complete sequence, 400 bp", S2),
        ("Microbacterium phage Gamma complete sequence, 320 bp", "ACGT" * 80),
        # Delta: same sequence as Alpha, no accession -> same_sequence merge
        ("Mycobacterium phage Delta complete sequence, 300 bp", S1),
        # no metadata row; header genus keeps it in scope
        ("Streptomyces phage Eps complete sequence, 300 bp", "GGCC" * 75),
        ("Arthrobacter phage complete sequence, 111 bp", "ACGTACGTAC"),  # malformed
        ("Mycobacterium phage Zeta complete sequence, 360 bp", "TTTT" * 90),
    ]))
    (d / "inphared" / "7Apr2026_millardlab_website_table.txt.gz").write_bytes(
        gzip.compress(_inphared_table([
            {"Accession": "JN000001", "Description": "Mycobacterium phage Alpha",
             "Host": "Mycobacterium"},
            {"Accession": "NC_999001.1", "Description": "Mycobacterium virus AlphaRef",
             "Host": "Mycobacterium"},  # RefSeq twin: identical seq, merges
            {"Accession": "JN000002", "Description": "Gordonia phage Beta",
             "Host": "Gordonia"},  # actino host but not myco -> out of scope
            {"Accession": "JN000003", "Description": "Microbacterium phage Gamma",
             "Host": "Microbacterium"},  # non-myco host + non-myco desc -> OUT of scope
            {"Accession": "PV777777", "Description": "Mycobacterium phage Future",
             "Host": "Mycobacterium"},  # no sequence in fasta -> exclusion
            {"Accession": "AB000001", "Description": "Escherichia phage Lambda",
             "Host": "Escherichia"},   # out of scope
        ]).encode()))
    (d / "inphared" / "7Apr2026_genomes.fa.gz").write_bytes(gzip.compress(_fasta([
        ("JN000001.1 Mycobacterium phage Alpha", S1),
        ("NC_999001.1 Mycobacterium virus AlphaRef", S1),
        ("JN000002.1 Gordonia phage Beta", S2),
        ("JN000003.1 Microbacterium phage Gamma", "ACGT" * 80),
        ("AB000001.1 Escherichia phage Lambda", "ACGT" * 100),
    ]).encode()))
    (d / "ncbi" / "esummary_00000.json").write_text(json.dumps({"result": {
        "uids": ["1", "2", "3", "4", "5", "6", "7"],
        "1": dict(_summ(1, "JN000001.1", "Mycobacterium phage Alpha, complete genome", 300)),
        "2": dict(_summ(2, "NC_999001.1", "Mycobacterium virus AlphaRef, complete genome", 300)),
        "3": dict(_summ(3, "JN000002.1", "Gordonia phage Beta, complete genome", 400)),
        "4": dict(_summ(4, "JN000003.1", "Microbacterium phage Gamma, complete genome", 320)),
        "5": dict(_summ(5, "AE999999.1", "Escherichia coli O157:H7 str. X, complete genome", 5500000)),
        "6": dict(_summ(6, "AB123456.1", "Mycobacteriophage L5 tail fiber gene, complete cds", 800)),
        "7": dict(_summ(7, "LR999999.1", "Pacmanvirus A23 genome assembly, complete genome", 395405)),
    }}))
    (d / "ncbi" / "efetch_fasta_00000.txt").write_text(_fasta([
        ("JN000001.1 Mycobacterium phage Alpha, complete genome", S1),
        ("NC_999001.1 Mycobacterium virus AlphaRef, complete genome", S1),
        ("JN000002.1 Gordonia phage Beta, complete genome", S2),
        ("JN000003.1 Microbacterium phage Gamma, complete genome", "ACGT" * 80),
        ("AE999999.1 Escherichia coli O157:H7 str. X, complete genome", "ACGT" * 1000),
        ("AB123456.1 Mycobacteriophage L5 tail fiber gene, complete cds", "ACGT" * 200),
        ("LR999999.1 Pacmanvirus A23 genome assembly, complete genome", "ACGT" * 5000),
    ]))
    (d / "ncbi" / "gb_remainder.txt").write_text("")
    (d / "fetch_manifest.json").write_text(json.dumps(FM, indent=1))
    return d


def load(cache):
    import glob as _g
    pdb_meta = pl.parse_phagesdb_tsv(str(cache / "phagesdb" / "phagesdb_sequenced_phages_full.tsv"))
    pdb_fasta, pdb_extras = pl.parse_phagesdb_fasta(str(cache / "phagesdb" / "Actinobacteriophages-All.fasta"))
    inp_rows = pl.parse_inphared_table(str(cache / "inphared" / "7Apr2026_millardlab_website_table.txt.gz"))
    wanted = {r["accession"] for r in inp_rows if pl.inphared_scope(r)}
    inp_seqs = pl.load_inphared_seqs(str(cache / "inphared" / "7Apr2026_genomes.fa.gz"), wanted)
    summ = pl.parse_ncbi_esummary(_g.glob(str(cache / "ncbi" / "esummary_*.json")))
    ncbi_fa = pl.parse_ncbi_fasta(_g.glob(str(cache / "ncbi" / "efetch_fasta_*.txt")))
    gb = pl.parse_gb_remainder(str(cache / "ncbi" / "gb_remainder.txt"))
    return pl.merge_records(pdb_meta, pdb_fasta, pdb_extras, inp_rows, inp_seqs,
                            summ, ncbi_fa, TAXONOMY, gb, FM)


# ---------------------------------------------------------------------------
# parser tests

def test_accession_helpers():
    assert pl.accession_base("NC_001900.1") == "NC_001900"
    assert pl.accession_base("KJ410132") == "KJ410132"
    assert pl.is_versioned_accession("NC_001900.1")
    assert not pl.is_versioned_accession("NC_001900")
    assert not pl.is_versioned_accession("phage_name")


def test_phagesdb_fasta_header_variants(cache):
    recs, extras = pl.parse_phagesdb_fasta(str(cache / "phagesdb" / "Actinobacteriophages-All.fasta"))
    assert "Alpha" in recs and "Beta" in recs and "Gamma" in recs and "Zeta" in recs
    assert recs["Gamma"]["header_host_genus"] == "Microbacterium"
    issues = sorted(x["issue"] for x in extras)
    assert issues == ["duplicate_name", "malformed_header"]
    dup = [x for x in extras if x["issue"] == "duplicate_name"][0]
    assert dup["name"] == "Alpha"
    malformed = [x for x in extras if x["issue"] == "malformed_header"][0]
    assert malformed["header"].startswith("Arthrobacter phage complete")


def test_inphared_scope_and_matching(cache):
    rows = pl.parse_inphared_table(str(cache / "inphared" / "7Apr2026_millardlab_website_table.txt.gz"))
    assert len(rows) == 6
    accs = {r["accession"] for r in rows}
    assert accs == {"JN000001", "NC_999001.1", "JN000002", "JN000003",
                    "PV777777", "AB000001"}
    in_scope = [r for r in rows if pl.inphared_scope(r)]
    assert {r["accession"] for r in in_scope} == {"JN000001", "NC_999001.1", "PV777777"}
    wanted = {r["accession"] for r in in_scope}
    seqs = pl.load_inphared_seqs(str(cache / "inphared" / "7Apr2026_genomes.fa.gz"), wanted)
    # base-keyed: versioned fasta id matches unversioned table accession
    assert set(seqs) == {"JN000001", "NC_999001"}
    assert seqs["JN000001"]["seq"] == S1


# ---------------------------------------------------------------------------
# merge tests

def test_merge_join_and_dedup(cache):
    m = load(cache)
    recs = {r["primary_accession"]: r for r in m["records"]}
    # Alpha: JN000001 joined across ncbi+inphared+phagesdb; RefSeq twin merged
    a = recs["JN000001.1"]
    assert a["sources"] == "ncbi;inphared;phagesdb"
    assert a["host_reported"] == "Mycobacterium smegmatis mc²155"
    assert a["host_taxid"] == "1763"
    assert a["host_evidence"] == "reported_host_label_only"
    assert a["evidence_class"] == "isolated_sequenced"
    assert a["completeness"] == "finished_sequence"
    assert a["topology"] == "linear"
    assert a["scope"] == "mycobacteriophage"
    assert a["in_mycobacteriaceae_host"] == "true"
    assert a["literature_reference"] == ""
    assert a["phagesdb_cluster"] == "A"
    # RefSeq twin merged as same_sequence duplicate with alias rows
    aliases = {(c["alias"], c["relation"]) for c in m["crosswalk"]}
    assert ("NC_999001.1", "same_sequence") in aliases
    # exact phagesdb duplicate merged
    assert ("phagesdb_name:Delta", "same_sequence") in aliases
    # Beta joins ncbi+phagesdb (Gordonia -> outside INPHARED myco scope)
    b = recs["JN000002.1"]
    assert b["sources"] == "ncbi;phagesdb"
    assert b["scope"] == "other_actinobacteriophage"
    assert b["host_genus"] == "Gordonia"
    assert b["host_taxid"] == "79255"
    # Gamma: phagesdb metadata host Microbacterium + ncbi record
    g = recs["JN000003.1"]
    assert g["sources"] == "ncbi;phagesdb"
    assert g["host_genus"] == "Microbacterium"
    assert g["host_taxid"] == ""  # no taxonomy entry in fixture
    # Zeta: phagesdb-only (no accession), host genus from FASTA header
    z = [r for r in m["records"] if r["phage_name"] == "Zeta"][0]
    assert z["primary_source"] == "phagesdb"
    assert z["host_genus"] == "Mycobacterium"
    assert z["scope"] == "mycobacteriophage"
    assert z["primary_source"] == "phagesdb"


def test_merge_exclusions_explicit_reasons(cache):
    m = load(cache)
    reasons = {(e["source"], e["reason"]) for e in m["exclusions"]}
    # NoSeq metadata row (no bulk fasta record)
    assert ("phagesdb_metadata", "no_sequence_in_bulk_fasta") in reasons
    # PV777777 listed in table but missing from release fasta
    assert ("inphared", "no_sequence_in_release_fasta") in reasons
    # Delta exact duplicate
    assert ("phagesdb_fasta", "exact_duplicate_in_bulk_fasta") in reasons
    # malformed header
    assert ("phagesdb_fasta", "malformed_header_missing_name") in reasons
    # Eps: metadata row missing + Streptomyces resolvable from header -> retained
    recs = {r["phage_name"]: r for r in m["records"]}
    assert "Eps" in recs and recs["Eps"]["scope"] == "other_actinobacteriophage"
    # NCBI screening rules
    assert ("ncbi", "non_phage_record_matched_search") in reasons      # E. coli genome
    assert ("ncbi", "short_fragment_not_genome_scale") in reasons      # gene fragment
    assert ("ncbi", "no_actinobacteriophage_evidence") in reasons      # Pacmanvirus
    assert "AE999999.1" not in {r["primary_accession"] for r in m["records"]}


def test_host_label_never_claimed_as_verified(cache):
    m = load(cache)
    for r in m["records"]:
        assert r["host_evidence"] in (pl.HOST_EVIDENCE_NOTE, pl.HOST_EVIDENCE_GB)
        assert r["host_reported_source"] != "" or r["host_reported"] == ""


# ---------------------------------------------------------------------------
# emit tests

def test_panel_ids_deterministic_and_unique(cache):
    m1 = load(cache)
    m2 = load(cache)
    for m in (m1, m2):
        pl.assign_panel_ids(m["records"])
    ids1 = [r["panel_id"] for r in m1["records"]]
    ids2 = [r["panel_id"] for r in m2["records"]]
    assert ids1 == ids2
    assert len(set(ids1)) == len(ids1)
    assert all(i.startswith("PHPUB-") for i in ids1)
    # mycobacteriophages sort first
    scopes = [r["scope"] for r in sorted(m1["records"], key=lambda r: r["panel_id"])]
    first_other = scopes.index("other_actinobacteriophage")
    assert "mycobacteriophage" not in scopes[first_other:]


def test_manifest_roundtrip_and_fasta_prefix(cache, tmp_path):
    m = load(cache)
    pl.assign_panel_ids(m["records"])
    rows = pl.manifest_rows(m["records"])
    by_id = {r["panel_id"]: r for r in rows}
    # every row has required provenance fields
    for r in rows:
        for c in ("panel_id", "primary_accession", "retrieval_timestamp",
                  "seq_sha256", "fasta_bytes_sha256", "length_bp",
                  "evidence_class", "sources", "source_url"):
            assert r[c], f"{c} empty on {r['panel_id']}"
        assert pl.is_versioned_accession(r["primary_accession"]) or \
            r["primary_source"] != "ncbi"
    # fasta bytes: > prefix, ascii, round-trip sha
    for r in rows:
        fb = pl.fasta_bytes(r)
        assert fb.startswith(b">")
        fb.decode("ascii")
        assert hashlib.sha256(fb).hexdigest() == r["fasta_bytes_sha256"]


def test_full_build_deterministic(cache, tmp_path):
    data_dir = str(tmp_path / "data")
    os.makedirs(data_dir)
    os.symlink(str(cache), os.path.join(data_dir, "cache"))
    out1 = str(tmp_path / "panel1")
    out2 = str(tmp_path / "panel2")
    bp.build(data_dir, out1)
    bp.build(data_dir, out2)
    for name in ["panel.full.fasta", "manifest_panel.tsv", "crosswalk.tsv",
                 "exclusions.tsv", "source_summary.tsv", "host_distribution.tsv",
                 "REPORT.md", "PROVENANCE.md"]:
        assert open(os.path.join(out1, name), "rb").read() == \
            open(os.path.join(out2, name), "rb").read(), name


def test_literature_reference_flag(cache):
    m = load(cache)
    lit = [r for r in m["records"] if r["literature_reference"]]
    # fixture has no reference phages; real L5/D29 hit LITERATURE_REFERENCE_PHAGES
    assert lit == []


def test_gb_remainder_host_qualifier(cache, tmp_path):
    gb_path = str(cache / "ncbi" / "gb_remainder.txt")
    open(gb_path, "w").write(
        "LOCUS       JN777777             300 bp    DNA     linear   BCT 17-AUG-2026\n"
        "ACCESSION   JN777777\nVERSION     JN777777.1\n"
        "FEATURES\n  source\n /host=\"Mycobacterium smegmatis\"\n"
        '                     /note="COMPLETENESS: full length"\n//\n')
    summ_extra = str(cache / "ncbi" / "esummary_00100.json")
    open(summ_extra, "w").write(json.dumps({"result": {
        "uids": ["8"],
        "8": dict(_summ(8, "JN777777.1", "Mycobacterium phage Ghost, complete genome", 15000)),
    }}))
    fa_extra = str(cache / "ncbi" / "efetch_fasta_00100.txt")
    open(fa_extra, "w").write(
        _fasta([("JN777777.1 Mycobacterium phage Ghost, complete genome", "ACGT" * 3750)]))
    m = load(cache)
    g = {r["primary_accession"]: r for r in m["records"]}.get("JN777777.1")
    assert g is not None
    assert g["gb_host_qualifier"] == "Mycobacterium smegmatis"
    assert g["host_evidence"] == pl.HOST_EVIDENCE_GB
    assert g["completeness"] == "full length"
    assert g["topology"] == "linear" and g["topology_source"] == "genbank_locus"
    assert g["scope"] == "mycobacteriophage"


def test_seq_sha256_case_and_whitespace_invariant():
    a = pl.seq_sha256("acgt ACGT\n acgt")
    b = pl.seq_sha256("ACGTACGTACGT")
    assert a == b
    assert a != pl.seq_sha256("ACGTACGTACGA")
