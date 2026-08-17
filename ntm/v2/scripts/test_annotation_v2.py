#!/usr/bin/env python3
"""
test_annotation_v2.py — focused tests for the NTM v2 functional-QC pipeline:

  * key-module detection (terminase/portal/capsid/tail/integrase/lysis) from
    Pharokka `annot` text and PHROG `category` mapping
  * graded candidate-functionality classification boundaries, including that
    missing modules stay 'not detected' rather than 'absent'
  * gene / truncation / hypothetical counting
  * cohort derivation (ml_reconstructed vs ml_singleton vs ancestral)
  * report determinism (byte-identical reruns)
  * driver input preparation: FASTA/index round-trip is exact + duplicate-free
    (tiny synthetic release fixture, plus the real 3,639-record release when
    present on this machine)

Run:  python3 -m pytest ntm/v2/scripts/test_annotation_v2.py -v
"""
from __future__ import annotations

import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import build_annotation_report_v2 as rep  # noqa: E402
import run_annotation_v2 as drv  # noqa: E402

REAL_RELEASE = os.path.isdir(drv.RELEASE_DIR)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

def write_lines(path, lines):
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


@pytest.fixture()
def fixture_root(tmp_path):
    """Synthetic 4-genome v2 analysis root with pharokka + checkv outputs."""
    index = [
        "\t".join(["genome_id", "source", "cohort", "clade_id", "length",
                   "status", "host_clades", "origin_fa"]),
        "\t".join(["ntm2_0_0001_ML", "ntm2_ml", "ntm2_ml_reconstructed",
                   "0_0001", "50000", "ml", "host_clade_0002", "ml.fa"]),
        "\t".join(["ntm2_0_0002_ML", "ntm2_ml", "ntm2_ml_singleton",
                   "0_0002", "8000", "singleton", "host_clade_0009", "ml.fa"]),
        "\t".join(["ntm2_0_0001_ANCESTRAL", "ntm2_anc", "ntm2_anc",
                   "0_0001", "48000", "ancestral", "host_clade_0002", "anc.fa"]),
        "\t".join(["ntm2_0_0003_ML", "ntm2_ml", "ntm2_ml_reconstructed",
                   "0_0003", "30000", "ml", "host_clade_0004", "ml.fa"]),
    ]

    cds_header = ("gene\tstart\tstop\tstrand\tcontig\tscore\tpartial\t"
                  "mmseqs_phrog\tmmseqs_alnScore\tmmseqs_seqIdentity\t"
                  "mmseqs_eVal\tpyhmmer_phrog\tannot\tcategory\ttransl_table")

    def cds(gid, n, annot, category, partial="00"):
        gene = f"{gid}_CDS_{n:04d}"
        return f"{gene}\t1\t100\t+\t{gid}\t66.2\t{partial}\t1724\t307\t0.6\t1e-90\tNone\t{annot}\t{category}\t11"

    cds_rows = [
        cds_header,
        # G1: full structural set + lysis + integrase, no truncation
        cds("ntm2_0_0001_ML", 1, "terminase large subunit", "head and packaging"),
        cds("ntm2_0_0001_ML", 2, "portal protein", "head and packaging"),
        cds("ntm2_0_0001_ML", 3, "major capsid protein", "head and packaging"),
        cds("ntm2_0_0001_ML", 4, "tail assembly chaperone", "tail"),
        cds("ntm2_0_0001_ML", 5, "endolysin", "lysis"),
        cds("ntm2_0_0001_ML", 6, "integrase", "integration and excision"),
        cds("ntm2_0_0001_ML", 7, "hypothetical protein", "unknown function"),
        cds("ntm2_0_0001_ML", 8, "hypothetical protein", "unknown function"),
        # G2 (singleton): only a terminase + 1 partial gene
        cds("ntm2_0_0002_ML", 1, "terminase small subunit", "head and packaging"),
        cds("ntm2_0_0002_ML", 2, "hypothetical protein", "unknown function", partial="10"),
        # G3 (ancestral): no key modules, all hypothetical
        cds("ntm2_0_0001_ANCESTRAL", 1, "hypothetical protein", "unknown function"),
        cds("ntm2_0_0001_ANCESTRAL", 2, "hypothetical protein", "unknown function"),
        # G4: capsid + tail + holin but CheckV low-quality
        cds("ntm2_0_0003_ML", 1, "major capsid protein", "head and packaging"),
        cds("ntm2_0_0003_ML", 2, "tail tube protein", "tail"),
        cds("ntm2_0_0003_ML", 3, "holin", "lysis"),
        cds("ntm2_0_0003_ML", 4, "DNA polymerase", "DNA, replication and repair"),
    ]

    lgd_header = "contig\tlength\tgc_perc\ttransl_table\tcds_coding_density"
    lgd = [lgd_header] + [
        f"{r.split(chr(9))[4]}\t{r.split(chr(9))[0].split('_CDS')[0] and ''}50000\t0.62\t11\t71.0"
        for r in []
    ] + [
        "\t".join(["ntm2_0_0001_ML", "50000", "0.62", "11", "71.0"]),
        "\t".join(["ntm2_0_0002_ML", "8000", "0.65", "11", "69.0"]),
        "\t".join(["ntm2_0_0001_ANCESTRAL", "48000", "0.61", "11", "70.0"]),
        "\t".join(["ntm2_0_0003_ML", "30000", "0.63", "11", "72.0"]),
    ]

    checkv_header = ("contig_id\tcontig_length\tprovirus\tproviral_length\t"
                     "gene_count\tviral_genes\thost_genes\tcheckv_quality\t"
                     "miuvig_quality\tcompleteness\tcompleteness_method\t"
                     "contamination\tkmer_freq\twarnings")
    checkv = [checkv_header,
              "\t".join(["ntm2_0_0001_ML", "50000", "No", "NA", "8", "6", "1",
                         "High-quality", "High-quality", "100.0",
                         "AAI-based (high-confidence)", "0.0", "1.78", ""]),
              "\t".join(["ntm2_0_0002_ML", "8000", "No", "NA", "2", "1", "0",
                         "Low-quality", "Low-quality", "42.1",
                         "AAI-based (low-confidence)", "0.0", "1.78", ""]),
              "\t".join(["ntm2_0_0001_ANCESTRAL", "48000", "No", "NA", "2", "0", "0",
                         "Not-determined", "Not-determined", "NA", "NA", "NA",
                         "NA", ""]),
              "\t".join(["ntm2_0_0003_ML", "30000", "No", "NA", "4", "3", "0",
                         "Low-quality", "Low-quality", "55.0",
                         "AAI-based (low-confidence)", "0.0", "1.78", ""])]

    write_lines(tmp_path / "genome_index.tsv", index)
    write_lines(tmp_path / "pharokka_cds_final_merged_output.tsv", cds_rows)
    write_lines(tmp_path / "pharokka_length_gc_cds_density.tsv", lgd)
    write_lines(tmp_path / "quality_summary.tsv", checkv)
    return tmp_path


def build(root):
    return rep.build_rows(
        str(root / "genome_index.tsv"),
        str(root / "pharokka_cds_final_merged_output.tsv"),
        str(root / "pharokka_length_gc_cds_density.tsv"),
        str(root / "quality_summary.tsv"),
    )


# ---------------------------------------------------------------------------
# key-module detection
# ---------------------------------------------------------------------------

class TestKeyModules:
    def test_full_module_set_detected(self, fixture_root):
        rows, _ = build(fixture_root)
        g1 = next(r for r in rows if r["genome_id"] == "ntm2_0_0001_ML")
        for m in ("terminase", "portal", "capsid", "tail", "integrase", "lysis"):
            assert g1[m] == 1, m
        assert g1["n_key_proteins"] == 6
        assert g1["head_packaging_modules"] == 3

    def test_category_only_detection(self):
        """PHROG categories infer modules even with vague annot text."""
        genes, partial, keys, cats, hypoth = rep.parse_pharokka_cds.__wrapped__ \
            if hasattr(rep.parse_pharokka_cds, "__wrapped__") else (None,) * 5
        # direct check via a temp file
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "cds.tsv")
            write_lines(p, [
                "gene\tstart\tstop\tstrand\tcontig\tscore\tpartial\tannot\tcategory",
                "g1\t1\t99\t+\tntm2_x_ML\t66\t00\tprotein of unknown function\tlysis",
                "g2\t1\t99\t-\tntm2_x_ML\t66\t00\tprotein of unknown function\tintegration and excision",
                "g3\t1\t99\t-\tntm2_x_ML\t66\t00\tprotein of unknown function\ttail",
            ])
            genes, partial, keys, cats, hypoth = rep.parse_pharokka_cds(p)
            assert keys["ntm2_x_ML"] == {"lysis", "integrase", "tail"}

    def test_rna_polymerase_is_not_lysis(self):
        """v1 regex bug: bare 'R' matched 'RNA polymerase'. Must not reappear."""
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "cds.tsv")
            write_lines(p, [
                "gene\tstart\tstop\tstrand\tcontig\tscore\tpartial\tannot\tcategory",
                "g1\t1\t99\t+\tntm2_x_ML\t66\t00\tRNA polymerase\tDNA, replication and repair",
                "g2\t1\t99\t+\tntm2_x_ML\t66\t00\tRz protein\tlysis",
            ])
            _g, _p, keys, _c, _h = rep.parse_pharokka_cds(p)
            assert "lysis" in keys["ntm2_x_ML"]          # via Rz
            # without Rz gene, RNA polymerase alone must NOT call lysis
            write_lines(p, [
                "gene\tstart\tstop\tstrand\tcontig\tscore\tpartial\tannot\tcategory",
                "g1\t1\t99\t+\tntm2_x_ML\t66\t00\tRNA polymerase\tDNA, replication and repair",
            ])
            _g, _p, keys, _c, _h = rep.parse_pharokka_cds(p)
            assert "lysis" not in keys["ntm2_x_ML"]

    def test_hypothetical_and_truncation_counting(self, fixture_root):
        rows, _ = build(fixture_root)
        g1 = next(r for r in rows if r["genome_id"] == "ntm2_0_0001_ML")
        assert g1["gene_count"] == 8
        assert g1["hypothetical"] == 2
        assert g1["hypothetical_frac"] == 0.25
        assert g1["truncated_genes"] == 0
        g2 = next(r for r in rows if r["genome_id"] == "ntm2_0_0002_ML")
        assert g2["gene_count"] == 2
        assert g2["truncated_genes"] == 1
        assert g2["truncated_frac"] == 0.5


# ---------------------------------------------------------------------------
# classification
# ---------------------------------------------------------------------------

class TestClassification:
    def test_tier_A_strong(self, fixture_root):
        rows, _ = build(fixture_root)
        g1 = next(r for r in rows if r["genome_id"] == "ntm2_0_0001_ML")
        assert g1["candidate_functionality"] == "A_strong_candidate"
        assert "all tier-A criteria met" in g1["tier_rationale"]

    def test_tier_C_low_checkv_with_modules(self, fixture_root):
        """capsid+tail+holin but CheckV Low-quality -> C, modules recorded."""
        rows, _ = build(fixture_root)
        g4 = next(r for r in rows if r["genome_id"] == "ntm2_0_0003_ML")
        assert g4["candidate_functionality"] == "C_partial_module_evidence"
        assert g4["capsid"] == 1 and g4["tail"] == 1 and g4["lysis"] == 1
        assert "below Medium" in g4["tier_rationale"]

    def test_D_is_not_detected_not_absent(self, fixture_root):
        rows, _ = build(fixture_root)
        g3 = next(r for r in rows if r["genome_id"] == "ntm2_0_0001_ANCESTRAL")
        assert g3["candidate_functionality"] == "D_modules_not_detected"
        assert "NOT proof of biological absence" in g3["tier_rationale"]
        # explicitly must not use absence language
        assert "non-functional" not in g3["candidate_functionality"]
        assert "absent" not in g3["candidate_functionality"].lower().replace("not_", "")

    def test_integrase_not_required_for_top_tier(self):
        tier, _ = rep.classify(
            gene_count=50,
            modules={"terminase", "portal", "capsid", "tail", "lysis"},  # no integrase
            checkv_grade="Medium-quality", truncated_frac=0.0, contamination=0.0)
        assert tier == "A_strong_candidate"

    def test_contamination_breaks_tier_A(self):
        tier, rat = rep.classify(
            gene_count=50,
            modules={"terminase", "portal", "capsid", "tail", "lysis"},
            checkv_grade="High-quality", truncated_frac=0.0, contamination=1.5)
        assert tier == "C_partial_module_evidence"
        assert "contamination" in rat

    def test_heavy_contamination_is_own_reason(self):
        tier, rat = rep.classify(
            gene_count=50, modules={"terminase", "capsid"},
            checkv_grade="High-quality", truncated_frac=0.0, contamination=20.0)
        assert "unreliable" in rat

    def test_high_truncation_blocks_A(self):
        tier, _ = rep.classify(
            gene_count=50,
            modules={"terminase", "portal", "capsid", "tail", "lysis"},
            checkv_grade="Complete", truncated_frac=0.5, contamination=0.0)
        assert tier != "A_strong_candidate"

    def test_single_head_module_with_tail_or_lysis_is_B(self):
        tier, _ = rep.classify(
            gene_count=40, modules={"capsid", "tail"},
            checkv_grade="High-quality", truncated_frac=0.0, contamination=0.0)
        assert tier == "B_moderate_candidate"

    def test_zero_genes_is_no_annotation(self):
        tier, rat = rep.classify(
            gene_count=0, modules=set(), checkv_grade="High-quality",
            truncated_frac=None, contamination=0.0)
        assert tier == "no_annotation"
        assert "not a biological claim" in rat

    def test_missing_checkv_fields_do_not_crash(self):
        tier, _ = rep.classify(
            gene_count=30, modules={"capsid"},
            checkv_grade="Not-determined", truncated_frac=0.1,
            contamination=None, contamination_known=False)
        assert tier == "C_partial_module_evidence"


# ---------------------------------------------------------------------------
# cohorts, summaries, determinism
# ---------------------------------------------------------------------------

class TestCohortsAndSummaries:
    def test_cohort_split(self, fixture_root):
        rows, _ = build(fixture_root)
        by = {r["genome_id"]: r["cohort"] for r in rows}
        assert by["ntm2_0_0001_ML"] == "ntm2_ml_reconstructed"
        assert by["ntm2_0_0002_ML"] == "ntm2_ml_singleton"
        assert by["ntm2_0_0001_ANCESTRAL"] == "ntm2_anc"
        assert by["ntm2_0_0003_ML"] == "ntm2_ml_reconstructed"

    def test_cohort_of_fallback_on_id_suffix(self):
        assert rep.cohort_of("ntm2_7_0007_ML", {}) == "ntm2_ml_reconstructed"
        assert rep.cohort_of("ntm2_7_0007_ANCESTRAL", {}) == "ntm2_anc"
        with pytest.raises(ValueError):
            rep.cohort_of("ecoli_x", {})

    def test_summaries_and_matrix(self, fixture_root):
        rows, _ = build(fixture_root)
        s = {r["group"]: r for r in rep.summarize(rows, "cohort")}
        assert s["ntm2_ml_reconstructed"]["n_genomes"] == 2
        assert s["ntm2_ml_singleton"]["n_genomes"] == 1
        assert s["ntm2_anc"]["n_genomes"] == 1
        assert s["ntm2_ml_reconstructed"]["pct_tier_A"] == 50.0  # 1 of 2
        mat = rep.tier_matrix(rows)
        assert mat[0] == ["cohort"] + rep.TIERS + ["total"]
        totals = {row[0]: row[-1] for row in mat[1:]}
        assert totals["ntm2_anc"] == 1

    def test_report_deterministic(self, fixture_root, tmp_path):
        outs = []
        for i in (1, 2):
            out = tmp_path / f"rep{i}"
            rows, _ = build(fixture_root)
            rep.write_tsv(str(out / "per_genome.tsv"), rep.PER_GENOME_COLS, rows)
            with open(out / "per_genome.tsv", "rb") as f:
                outs.append(f.read())
        assert outs[0] == outs[1]
        # sorted by genome_id
        text = outs[0].decode()
        ids = [ln.split("\t")[0] for ln in text.splitlines()[1:]]
        assert ids == sorted(ids)

    def test_missing_checkv_aborts_by_default(self, fixture_root):
        (fixture_root / "quality_summary.tsv").write_text(
            "contig_id\tcheckv_quality\nntm2_0_0001_ML\tHigh-quality\n")
        with pytest.raises(RuntimeError, match="missing from CheckV"):
            build(fixture_root)
        rows, missing = rep.build_rows(
            str(fixture_root / "genome_index.tsv"),
            str(fixture_root / "pharokka_cds_final_merged_output.tsv"),
            str(fixture_root / "pharokka_length_gc_cds_density.tsv"),
            str(fixture_root / "quality_summary.tsv"),
            missing_checkv_ok=True)
        assert len(missing) == 3
        assert all(r["checkv_grade"] == "Not-determined" for r in rows
                   if r["genome_id"] != "ntm2_0_0001_ML")


# ---------------------------------------------------------------------------
# driver: input preparation round-trip (synthetic + real release)
# ---------------------------------------------------------------------------

def make_release(root):
    ml = root / "all_ntm_ml_phage_genomes.fa"
    anc = root / "all_ntm_ancestral_phage_genomes.fa"
    write_lines(ml, [
        ">ntm2_0_0001_ML status=ml n_members=5 length=30 host_clades=host_clade_0002",
        "ACGTACGTAC",
        ">ntm2_0_0002_ML status=singleton n_members=1 length=10 host_clades=host_clade_0009",
        "TTTTGGCCCCA",
    ])
    write_lines(anc, [
        ">ntm2_0_0001_ANCESTRAL host_clades=host_clade_0002 species=Mycobacterium",
        "GGGGCCTTTTAAAA",
    ])
    return ml, anc


class TestDriverPrepare:
    def test_roundtrip_exact_and_duplicate_free(self, tmp_path, monkeypatch):
        ml, anc = make_release(tmp_path)
        monkeypatch.setattr(drv, "RELEASE_DIR", str(tmp_path))
        monkeypatch.setattr(drv, "ML_FA", str(ml))
        monkeypatch.setattr(drv, "ANC_FA", str(anc))
        monkeypatch.setattr(drv, "EXPECT_ML", 2)
        monkeypatch.setattr(drv, "EXPECT_ANC", 1)
        monkeypatch.setattr(drv, "EXPECT_TOTAL", 3)
        out = tmp_path / "root"
        info = drv.stage_prepare(str(out))
        assert info["n_total"] == 3 and info["n_ml"] == 2 and info["n_anc"] == 1
        # bare-ID FASTA round-trips to index exactly
        ids = [g for g, _d, _s in drv.iter_fasta(info["input_fa"])]
        assert set(ids) == {"ntm2_0_0001_ML", "ntm2_0_0002_ML",
                            "ntm2_0_0001_ANCESTRAL"}
        with open(info["genome_index"]) as f:
            lines = f.read().splitlines()
        assert lines[0].startswith("genome_id\tsource\tcohort")
        srcs = [ln.split("\t")[1] for ln in lines[1:]]
        assert srcs == ["ntm2_ml", "ntm2_ml", "ntm2_anc"]
        # lengths preserved exactly
        lens = {ln.split("\t")[0]: int(ln.split("\t")[4]) for ln in lines[1:]}
        assert lens == {"ntm2_0_0001_ML": 10, "ntm2_0_0002_ML": 11,
                        "ntm2_0_0001_ANCESTRAL": 14}
        # cohort split incl. singleton detection from release status
        cohorts = {ln.split("\t")[0]: ln.split("\t")[2] for ln in lines[1:]}
        assert cohorts["ntm2_0_0002_ML"] == "ntm2_ml_singleton"
        assert cohorts["ntm2_0_0001_ML"] == "ntm2_ml_reconstructed"

    def test_duplicate_ids_abort(self, tmp_path, monkeypatch):
        ml, anc = make_release(tmp_path)
        with open(ml, "a", newline="\n") as f:  # duplicate of an ancestral ID
            f.write(">ntm2_0_0001_ANCESTRAL status=ml length=5\nACGTA\n")
        monkeypatch.setattr(drv, "RELEASE_DIR", str(tmp_path))
        monkeypatch.setattr(drv, "ML_FA", str(ml))
        monkeypatch.setattr(drv, "ANC_FA", str(anc))
        with pytest.raises(RuntimeError, match="duplicate genome IDs"):
            drv.stage_prepare(str(tmp_path / "root"))

    def test_wrong_counts_abort(self, tmp_path, monkeypatch):
        ml, anc = make_release(tmp_path)
        monkeypatch.setattr(drv, "RELEASE_DIR", str(tmp_path))
        monkeypatch.setattr(drv, "ML_FA", str(ml))
        monkeypatch.setattr(drv, "ANC_FA", str(anc))
        monkeypatch.setattr(drv, "EXPECT_ML", 99)  # expectation not met
        with pytest.raises(RuntimeError, match="expected 99 ML"):
            drv.stage_prepare(str(tmp_path / "root"))

    def test_thread_cap_enforced(self, tmp_path, capsys):
        with pytest.raises(SystemExit):
            drv.main(["--dry-run", "--threads", "128", "--root", str(tmp_path)])

    @pytest.mark.skipif(not REAL_RELEASE, reason="real v2 release not mounted")
    def test_real_release_roundtrip(self, tmp_path):
        """The validation criterion: exactly 2,388 ML + 1,251 ancestral =
        3,639 unique records, round-trip exact, duplicate-free."""
        info = drv.stage_prepare(str(tmp_path))
        assert info["n_ml"] == 2388
        assert info["n_anc"] == 1251
        assert info["n_total"] == 3639
        ids = [g for g, _d, _s in drv.iter_fasta(info["input_fa"])]
        assert len(ids) == 3639 and len(set(ids)) == 3639

    def test_guard_rejects_v1_tree(self, tmp_path):
        with pytest.raises(RuntimeError, match="OUTSIDE the v1"):
            drv.guard_paths(drv.V1_ANNOTATION_DIR)
        with pytest.raises(RuntimeError, match="OUTSIDE the v1"):
            drv.guard_paths(os.path.join(drv.V1_ANNOTATION_DIR, "v2new"))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
