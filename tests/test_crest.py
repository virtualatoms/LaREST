"""Tests for larest.crest CREST output parsing and execution."""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import pytest

from larest.constants import CALMOL_TO_JMOL
from larest.crest import (
    parse_crest_entropy_output,
    run_crest_confgen,
    run_crest_entropy,
)

# Expected values from tests/data/crest_entropy_output.txt
_S_CONF_CALMOL = 12.3456
_S_RRHO_CALMOL = 5.6789
_S_TOTAL_CALMOL = 18.0245


class TestParseCrestEntropyOutput:
    def test_parses_s_conf(self, crest_entropy_output_file):
        result = parse_crest_entropy_output(crest_entropy_output_file)
        assert result["S_conf"] == pytest.approx(
            _S_CONF_CALMOL * CALMOL_TO_JMOL,
            rel=1e-5,
        )

    def test_parses_s_rrho(self, crest_entropy_output_file):
        result = parse_crest_entropy_output(crest_entropy_output_file)
        assert result["S_rrho"] == pytest.approx(
            _S_RRHO_CALMOL * CALMOL_TO_JMOL,
            rel=1e-5,
        )

    def test_parses_s_total(self, crest_entropy_output_file):
        result = parse_crest_entropy_output(crest_entropy_output_file)
        assert result["S_total"] == pytest.approx(
            _S_TOTAL_CALMOL * CALMOL_TO_JMOL,
            rel=1e-5,
        )

    def test_returns_none_for_missing(self, tmp_path):
        empty = tmp_path / "empty.txt"
        empty.write_text("nothing here\n")
        result = parse_crest_entropy_output(empty)
        assert result["S_conf"] is None
        assert result["S_rrho"] is None
        assert result["S_total"] is None

    def test_keys_present(self, crest_entropy_output_file):
        result = parse_crest_entropy_output(crest_entropy_output_file)
        assert set(result.keys()) == {"S_conf", "S_rrho", "S_total"}


class TestRunCrestConfgenMocked:
    def _setup_rdkit_checkpoint(self, dir_path):
        """Create the minimal checkpoint files that run_crest_confgen needs."""
        xtb_rdkit_dir = dir_path / "xtb" / "rdkit"
        xtb_rdkit_dir.mkdir(parents=True)
        results_csv = xtb_rdkit_dir / "results.csv"
        results_csv.write_text("conformer_id,H,S,G\n0,-324567.89,-85.12,-350000.00\n")

        # Create a fake optimised xyz (crest reads this path)
        conformer_dir = xtb_rdkit_dir / "conformer_0"
        conformer_dir.mkdir()
        xyz = conformer_dir / "conformer_0.xtbopt.xyz"
        xyz.write_text("1\ntest\nC 0.0 0.0 0.0\n")
        return dir_path

    def test_runs_crest_and_xtb(self, tmp_path, minimal_config):
        self._setup_rdkit_checkpoint(tmp_path)

        fake_xtb_results = {"H": -100000.0, "S": -50.0, "G": -115000.0}

        with (
            patch("larest.crest.subprocess.run") as mock_crest,
            patch("larest.crest.run_xtb", return_value=fake_xtb_results) as mock_xtb,
        ):
            mock_crest.return_value = MagicMock(returncode=0)
            result = run_crest_confgen(tmp_path, minimal_config)

        mock_crest.assert_called_once()
        mock_xtb.assert_called_once()
        assert result == {"crest": fake_xtb_results}

    def test_skips_xtb_when_disabled(self, tmp_path, minimal_config):
        cfg = {**minimal_config, "steps": {**minimal_config["steps"], "xtb": False}}
        self._setup_rdkit_checkpoint(tmp_path)

        with (
            patch("larest.crest.subprocess.run"),
            patch("larest.crest.run_xtb") as mock_xtb,
        ):
            result = run_crest_confgen(tmp_path, cfg)

        mock_xtb.assert_not_called()
        assert result == {}

    def test_raises_when_no_rdkit_results(self, tmp_path, minimal_config):
        # No checkpoint files present
        with pytest.raises(FileNotFoundError):
            run_crest_confgen(tmp_path, minimal_config)


class TestRunCrestEntropyMocked:
    def _setup_censo_checkpoint(self, dir_path):
        censo_dir = dir_path / "censo"
        censo_dir.mkdir(parents=True)

        # censo.txt with a Highest ranked conformer line
        censo_txt = censo_dir / "censo.txt"
        censo_txt.write_text(
            "part0 -0.12 -0.11 x\n  Highest ranked conformer CONF1\n"
            "part1 -0.22 -0.21 x\n  Highest ranked conformer CONF1\n"
            "part2 -0.32 -0.31 x\n  Highest ranked conformer CONF1\n"
            "part3 -0.42 -0.41 x\n  Highest ranked conformer CONF1\n",
        )

        # 3_REFINEMENT.xyz with CONF1
        xyz = censo_dir / "3_REFINEMENT.xyz"
        xyz.write_text(
            "3\nCONF1 energy= -0.12345\n"
            "O   0.0 0.0 0.0\nH   1.0 0.0 0.0\nH   0.0 1.0 0.0\n",
        )

    def test_runs_crest_entropy(self, tmp_path, minimal_config):
        self._setup_censo_checkpoint(tmp_path)

        fake_results = {"S_conf": 51.7, "S_rrho": 23.8, "S_total": 75.5}

        with (
            patch("larest.crest.subprocess.run"),
            patch("larest.crest.parse_crest_entropy_output", return_value=fake_results),
        ):
            result = run_crest_entropy(tmp_path, minimal_config)

        assert result == fake_results

    def test_writes_results_json(self, tmp_path, minimal_config):
        self._setup_censo_checkpoint(tmp_path)
        fake_results = {"S_conf": 51.7, "S_rrho": 23.8, "S_total": 75.5}

        with (
            patch("larest.crest.subprocess.run"),
            patch("larest.crest.parse_crest_entropy_output", return_value=fake_results),
        ):
            run_crest_entropy(tmp_path, minimal_config)

        results_json = tmp_path / "crest_entropy" / "results.json"
        assert results_json.exists()
        data = json.loads(results_json.read_text())
        assert data["S_total"] == pytest.approx(75.5)

    def test_raises_when_no_censo_results(self, tmp_path, minimal_config):
        with pytest.raises(ValueError, match="Failed to extract best CENSO conformer"):
            run_crest_entropy(tmp_path, minimal_config)


# ---------------------------------------------------------------------------
# Integration tests (require crest installed)
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestRunCrestIntegration:
    """Runs the full CREST confgen stage on a tiny molecule."""

    def test_crest_confgen_produces_best_xyz(self, tmp_path, minimal_config):
        """Run rdkit first, then run crest confgen — check crest_best.xyz exists."""
        from larest.rdkit import run_rdkit

        run_rdkit("C1CC(=O)O1", tmp_path, minimal_config)
        run_crest_confgen(tmp_path, minimal_config)

        assert (tmp_path / "crest_confgen" / "crest_best.xyz").exists()
