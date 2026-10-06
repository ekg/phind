#!/usr/bin/env python3
"""
test_annotation_v3.py — focused tests for the NTM v3 functional-QC pipeline:

  * driver input preparation: FASTA/index round-trip is exact + duplicate-free
    (tiny synthetic release fixture), with ntm3_* clade-id / cohort derivation
  * duplicate IDs and wrong record counts abort
  * guard rejects any root inside the v1 evidence tree
  * GFF3 split: one file per genome, per-genome ##sequence-region preserved,
    embedded ##FASTA omitted, committed index maps genome_id -> path + count
  * split aborts on a seqid absent from the prepared input FASTA
  * report cohort derivation (ntm3_* labels)

Run:  python3 -m pytest ntm/v3/scripts/test_annotation_v3.py -v
"""
from __future__ import annotations

import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import build_annotation_report_v3 as rep  # noqa: E402
import run_annotation_v3 as drv  # noqa: E402


def write_lines(path, lines):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


@pytest.fixture()
def synthetic_release(tmp_path):
    ml = tmp_path / "all_ntm_v3_ml_phage_genomes.fa"
    anc = tmp_path / "all_ntm_v3_ancestral_phage_genomes.fa"
    write_lines(ml, [
        ">ntm3_0_0000_ML status=ml n_members=100 length=20",
        "ACGTACGTACGTACGTACGT",
        ">ntm3_0_0001_ML status=singleton n_members=1 length=6",
        "TTTTGG",
    ])
    write_lines(anc, [
        ">ntm3_0_0000_ANCESTRAL status=ml n_members=100 length=18",
        "GGGGCCTTTTAAAACCCC",
    ])
    return ml, anc


def patch_release(monkeypatch, ml, anc, n_ml=2, n_anc=1):
    monkeypatch.setattr(drv, "ML_FA", str(ml))
    monkeypatch.setattr(drv, "ANC_FA", str(anc))
    monkeypatch.setattr(drv, "EXPECT_ML", n_ml)
    monkeypatch.setattr(drv, "EXPECT_ANC", n_anc)
    monkeypatch.setattr(drv, "EXPECT_TOTAL", n_ml + n_anc)


class TestPrepare:
    def test_roundtrip_clade_and_cohort(self, tmp_path, monkeypatch, synthetic_release):
        ml, anc = synthetic_release
        patch_release(monkeypatch, ml, anc)
        info = drv.stage_prepare(str(tmp_path / "root"))
        assert info["n_total"] == 3 and info["n_ml"] == 2 and info["n_anc"] == 1
        with open(info["genome_index"]) as f:
            rows = [ln.split("\t") for ln in f.read().splitlines()]
        assert rows[0][:3] == ["genome_id", "source", "cohort"]
        by_id = {r[0]: r for r in rows[1:]}
        assert by_id["ntm3_0_0000_ML"][3] == "0_0000"          # clade_id
        assert by_id["ntm3_0_0000_ANCESTRAL"][3] == "0_0000"
        assert by_id["ntm3_0_0001_ML"][2] == "ntm3_ml_singleton"
        assert by_id["ntm3_0_0000_ML"][2] == "ntm3_ml_reconstructed"
        assert by_id["ntm3_0_0000_ANCESTRAL"][2] == "ntm3_anc"
        assert by_id["ntm3_0_0000_ML"][4] == "20"              # length
        assert by_id["ntm3_0_0000_ANCESTRAL"][5] == "ancestral"

    def test_duplicate_ids_abort(self, tmp_path, monkeypatch, synthetic_release):
        ml, anc = synthetic_release
        with open(ml, "a", newline="\n") as f:
            f.write(">ntm3_0_0000_ANCESTRAL status=ml length=18\nACGT\n")
        patch_release(monkeypatch, ml, anc)
        with pytest.raises(RuntimeError, match="duplicate genome IDs"):
            drv.stage_prepare(str(tmp_path / "root"))

    def test_wrong_counts_abort(self, tmp_path, monkeypatch, synthetic_release):
        ml, anc = synthetic_release
        patch_release(monkeypatch, ml, anc, n_ml=99)
        with pytest.raises(RuntimeError, match="expected 99 ML"):
            drv.stage_prepare(str(tmp_path / "root"))

    def test_thread_cap_enforced(self, tmp_path):
        with pytest.raises(SystemExit):
            drv.main(["--dry-run", "--threads", "128", "--root", str(tmp_path)])

    def test_guard_rejects_v1_tree(self, tmp_path):
        with pytest.raises(RuntimeError, match="OUTSIDE the v1"):
            drv.guard_paths(drv.V1_ANNOTATION_DIR)
        with pytest.raises(RuntimeError, match="OUTSIDE the v1"):
            drv.guard_paths(os.path.join(drv.V1_ANNOTATION_DIR, "v3new"))


class TestSplit:
    def _root(self, tmp_path):
        root = tmp_path / "root"
        write_lines(root / "input" / "all_v3_phage_genomes.fa", [
            ">ntm3_0_0000_ML", "ACGTACGTACGTACGTACGT",
            ">ntm3_0_0000_ANCESTRAL", "GGGGCCTTTTAAAACCCC",
            ">ntm3_0_0001_ML", "TTTTGG",
        ])
        write_lines(root / "pharokka_out" / "pharokka.gff", [
            "##gff-version 3",
            "##sequence-region ntm3_0_0000_ML 1 20",
            "##sequence-region ntm3_0_0000_ANCESTRAL 1 18",
            "##sequence-region ntm3_0_0001_ML 1 6",
            "\t".join(["ntm3_0_0000_ML", "pyrodigal-gv", "CDS", "2", "10",
                       ".", "+", "0", "ID=ntm3_0_0000_ML_CDS_0001"]),
            "\t".join(["ntm3_0_0000_ML", "pyrodigal-gv", "CDS", "12", "20",
                       ".", "-", "0", "ID=ntm3_0_0000_ML_CDS_0002"]),
            "\t".join(["ntm3_0_0000_ANCESTRAL", "pyrodigal-gv", "CDS", "1", "18",
                       ".", "+", "0", "ID=ntm3_0_0000_ANCESTRAL_CDS_0001"]),
            "##FASTA",
            ">ntm3_0_0000_ML", "ACGTACGTACGTACGTACGT",
        ])
        return root

    def test_split_three_genomes(self, tmp_path):
        root = self._root(tmp_path)
        idx = tmp_path / "annotation_gff_index.tsv"
        info = drv.split_gff(str(root), str(idx))
        assert info["n_genomes"] == 3
        assert info["n_features"] == 3
        assert (root / "pharokka_out" / "per_genome_gff" /
                "ntm3_0_0000_ML.gff").exists()
        # empty genome still gets a header-only GFF + a 0 count
        g1 = (root / "pharokka_out" / "per_genome_gff" / "ntm3_0_0001_ML.gff")
        assert g1.read_text() == "##gff-version 3\n" + \
            "##sequence-region ntm3_0_0001_ML 1 6\n"
        # embedded ##FASTA omitted, per-genome region preserved
        body = (root / "pharokka_out" / "per_genome_gff" /
                "ntm3_0_0000_ML.gff").read_text()
        assert body.startswith("##gff-version 3\n##sequence-region ntm3_0_0000_ML 1 20\n")
        assert "##FASTA" not in body
        assert "ACGTACGTACGTACGTACGT" not in body
        # committed index maps genome_id -> gff path + gene count
        with open(idx) as f:
            lines = f.read().splitlines()
        assert lines[0] == "genome_id\tgff_path\tgene_count"
        n = {ln.split("\t")[0]: ln.split("\t")[2] for ln in lines[1:]}
        assert n == {"ntm3_0_0000_ML": "2", "ntm3_0_0000_ANCESTRAL": "1",
                     "ntm3_0_0001_ML": "0"}
        # NVMe manifest mirror is identical
        assert (root / "pharokka_out" / "per_genome_gff" /
                "manifest.tsv").read_text() == idx.read_text()

    def test_split_unknown_seqid_aborts(self, tmp_path):
        root = self._root(tmp_path)
        gff = root / "pharokka_out" / "pharokka.gff"
        text = gff.read_text()
        bad = "\t".join(["ntm3_9_9999_ML", "pyrodigal-gv", "CDS", "1",
                         "3", ".", "+", "0", "ID=x"])
        gff.write_text(text.replace("##FASTA", bad + "\n##FASTA"))
        with pytest.raises(RuntimeError, match="not in the prepared input"):
            drv.split_gff(str(root), str(tmp_path / "idx.tsv"))


class TestReportCohorts:
    def test_ntm3_cohort_labels(self):
        assert rep.cohort_of("ntm3_7_0007_ML", {}) == "ntm3_ml_reconstructed"
        assert rep.cohort_of("ntm3_7_0007_ANCESTRAL", {}) == "ntm3_anc"
        assert rep.cohort_of("ntm3_x_ML", {"status": "singleton"}) == \
            "ntm3_ml_singleton"
        with pytest.raises(ValueError):
            rep.cohort_of("ntm2_7_0007_ML", {})


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
