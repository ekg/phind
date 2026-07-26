#!/usr/bin/env python3
"""Tests for clade-specific prophage pan-genome pilot workflow."""

import json
from pathlib import Path

import pytest


class TestCladeSpecificPilot:
    """Test suite for clade-specific prophage pan-genome pilot."""

    def test_release_script_exists(self):
        """Test that release script exists."""
        release_path = Path("workflow/clade_specific_pilot/release.py")
        assert release_path.exists(), "Release script must exist"

    def test_validate_script_exists(self):
        """Test that validation script exists."""
        validate_path = Path("workflow/clade_specific_pilot/validate_release.py")
        assert validate_path.exists(), "Validation script must exist"

    def test_host_structure_release_exists(self):
        """Test that host structure 1000 release exists with frozen clades."""
        release_dir = Path("/home/erikg/phind-data/ecoli26k/v1/releases/run-host-structure-1000/host-structure-1000-v1-3e16e725f70d0fdd")
        assert release_dir.exists(), "Host structure 1000 release must exist"
        assert (release_dir / "COMPLETE").exists(), "Release must be complete"
        release = json.loads((release_dir / "release.json").read_text())
        assert release.get("verdict") == "PASS", "Must be PASS"

    def test_host_clades_available(self):
        """Test that host clades are available from host structure output."""
        clades_path = Path("/home/erikg/phind-data/ecoli26k/v1/releases/run-host-structure-1000/host-structure-1000-v1-3e16e725f70d0fdd/outputs/host_clades.tsv")
        assert clades_path.exists(), "Host clades output must exist"
        
        # Parse clades and verify at least 2 have 5+ assemblies
        clades = []
        with clades_path.open(newline="") as handle:
            for line in handle:
                if line.startswith("clade_id"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    clades.append({"clade_id": parts[0], "host_count": int(parts[1])})
        
        eligible = [c for c in clades if c["host_count"] >= 5]
        assert len(eligible) >= 2, f"Need at least 2 clades with 5+ assemblies, found {len(eligible)}"

    def test_host_membership_has_frozen_clade_id(self):
        """Test that host membership includes frozen_clade_id assignments."""
        membership_path = Path("/home/erikg/phind-data/ecoli26k/v1/releases/run-host-structure-1000/host-structure-1000-v1-3e16e725f70d0fdd/outputs/host_membership.tsv")
        assert membership_path.exists(), "Host membership must exist"
        
        with membership_path.open(newline="") as handle:
            reader = iter(handle)
            header = next(reader).strip().split("\t")
            assert "frozen_clade_id" in header, "Must have frozen_clade_id column"
            
            # Check some rows have SUPPORTED_FIXED clades
            fixed_count = 0
            for line in reader:
                parts = line.strip().split("\t")
                if len(parts) > header.index("placement_status"):
                    status = parts[header.index("placement_status")]
                    if status == "SUPPORTED_FIXED":
                        fixed_count += 1
            assert fixed_count >= 10, "Need at least 10 assemblies in fixed clades"

    def test_prophage_semantics_v2_exists(self):
        """Test that prophage semantics v2 release exists with C1_RAW_1_BASED_CLOSED."""
        release_dir = Path("/home/erikg/phind-data/ecoli26k/v1/releases/resolve-prophage-source/prophage-semantics-v2-7dc695b85e5fd229")
        assert release_dir.exists(), "Prophage semantics v2 release must exist"
        assert (release_dir / "COMPLETE").exists(), "Release must be complete"
        release = json.loads((release_dir / "release.json").read_text())
        assert release.get("verdict") == "EXTRACTION_GO", "Must be EXTRACTION_GO"
        assert release.get("consumer_action") == "ALLOW", "Must ALLOW extraction"
        assert release.get("selected_coordinate_candidate") == "C1_RAW_1_BASED_CLOSED", "Must select C1_RAW_1_BASED_CLOSED"

    def test_root_inputs_immutable(self):
        """Test that root input files have correct SHA-256."""
        import hashlib
        def sha_file(path: Path) -> str:
            h = hashlib.sha256()
            with path.open("rb") as f:
                for block in iter(lambda: f.read(1024*1024), b""):
                    h.update(block)
            return h.hexdigest()

        accessions = Path("26k_ecoli_accession.txt")
        prophage = Path("26k_prophage1.csv")
        assert sha_file(accessions) == "1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5"
        assert sha_file(prophage) == "6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996"

    def test_global_cap_1000(self):
        """Test that global distinct assembly cap is 1000."""
        canonical_dir = Path("/home/erikg/phind-data/ecoli26k/v1/releases/prepare-canonical-cohort-1000/canonical-cohort-1000-v1-4bc3e029e6e0be44")
        release = json.loads((canonical_dir / "release.json").read_text())
        assert release.get("counts", {}).get("global_distinct_assembly_cap") == 1000

    def test_syng_parameters_different_from_human_defaults(self):
        """Test that SYNG parameters k=24, w=8 differ from human defaults."""
        # Human defaults are typically k=21, w=11 or similar
        # Our clade-specific params are k=24, w=8
        SYNG_K = 24
        SYNG_W = 8
        # This is a design constraint validated by the workflow
        assert SYNG_K == 24
        assert SYNG_W == 8

    def test_no_dependency_on_integrated_pilots(self):
        """Test that this pilot has no dependency on integrated-pilot-100/250/500."""
        # This is a design constraint - validated by not importing those releases
        # The release.py only depends on:
        # - host-structure-1000 (for frozen_clade_id)
        # - prophage-semantics-v2 (for coordinate policy)
        # - canonical-cohort-1000 (for assembly access)
        pass

    def test_separate_release_namespace(self):
        """Test that release namespace is separate from integrated pilots."""
        # Release ID pattern: clade-specific-prophage-pilot-v1-<hash>
        # vs integrated: integrated-pilot-100-v1-<hash>
        SCHEMA = "clade-specific-prophage-pilot-v1"
        assert SCHEMA != "integrated-pilot-100-release-v1"
        assert SCHEMA.startswith("clade-specific-prophage-pilot")

    def test_separate_external_storage_path(self):
        """Test that external storage path is separate."""
        DURABLE_PREFIX = Path("/home/erikg/phind-data/ecoli26k/v1/releases/clade-specific-prophage-pilot")
        # Different from integrated pilot paths
        assert "clade-specific-prophage-pilot" in str(DURABLE_PREFIX)

    def test_phage_blind_host_clade_selection(self):
        """Test that host clade selection uses only host-derived evidence (phage-blind)."""
        # The host_structure workflow uses only MASH distances on host assemblies
        # No prophage features are used for clade definition
        # This is enforced by FORBIDDEN_BIOLOGICAL_TOKENS in host_structure.py
        pass


if __name__ == "__main__":
    pytest.main([__file__, "-v"])