"""Tests for larest.main pipeline orchestration."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from larest.data import Initiator, MolResults, Monomer, Polymer
from larest.main import compile_results, run_pipeline


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_MONOMER_SMILES = "C1CC(=O)O1"
_INITIATOR_SMILES = "CCO"

_SECTIONS = [
    "rdkit",
    "crest",
    "censo_prescreening",
    "censo_screening",
    "censo_optimization",
    "censo_refinement",
]


def _make_mol_results(smiles: str, h: float = -100000.0, s: float = -50.0, g: float = -115000.0):
    sections = {section: {"H": h, "S": s, "G": g} for section in _SECTIONS}
    return MolResults(smiles=smiles, sections=sections)


# ---------------------------------------------------------------------------
# compile_results
# ---------------------------------------------------------------------------


class TestCompileResults:
    def test_creates_summary_directory(self, tmp_path):
        monomer_results = _make_mol_results(_MONOMER_SMILES)
        polymer_results = [(2, _make_mol_results("OCC(=O)OCC(=O)O", h=-200000.0, s=-100.0, g=-230000.0))]

        compile_results(
            monomer_smiles=_MONOMER_SMILES,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=None,
            output_dir=tmp_path,
            reaction_type="RER",
        )

        from larest.output import slugify

        summary_dir = tmp_path / "Monomer" / slugify(_MONOMER_SMILES) / "summary"
        assert summary_dir.exists()

    def test_creates_csv_per_section(self, tmp_path):
        monomer_results = _make_mol_results(_MONOMER_SMILES)
        polymer_results = [(2, _make_mol_results("OCC(=O)OCC(=O)O", h=-200000.0, s=-100.0, g=-230000.0))]

        compile_results(
            monomer_smiles=_MONOMER_SMILES,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=None,
            output_dir=tmp_path,
            reaction_type="RER",
        )

        from larest.output import slugify

        summary_dir = tmp_path / "Monomer" / slugify(_MONOMER_SMILES) / "summary"
        for section in _SECTIONS:
            assert (summary_dir / f"{section}.csv").exists()

    def test_delta_computation_rer(self, tmp_path):
        """delta_param = (poly - n * mono - 0) / n"""
        monomer_results = _make_mol_results(_MONOMER_SMILES, h=-100000.0, s=-50.0, g=-115000.0)
        n = 2
        polymer_results = [(n, _make_mol_results("P", h=-210000.0, s=-105.0, g=-241500.0))]

        compile_results(
            monomer_smiles=_MONOMER_SMILES,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=None,
            output_dir=tmp_path,
            reaction_type="RER",
        )

        import pandas as pd
        from larest.output import slugify

        summary_dir = tmp_path / "Monomer" / slugify(_MONOMER_SMILES) / "summary"
        df = pd.read_csv(summary_dir / "rdkit.csv")
        row = df[df["polymer_length"] == n].iloc[0]

        expected_delta_h = (-210000.0 - n * -100000.0 - 0) / n
        assert row["delta_H"] == pytest.approx(expected_delta_h, rel=1e-6)

    def test_ror_raises_without_initiator(self, tmp_path):
        monomer_results = _make_mol_results(_MONOMER_SMILES)
        polymer_results = [(1, _make_mol_results("P"))]

        with pytest.raises(ValueError, match="Initiator results are required"):
            compile_results(
                monomer_smiles=_MONOMER_SMILES,
                monomer_results=monomer_results,
                polymer_results=polymer_results,
                initiator_results=None,
                output_dir=tmp_path,
                reaction_type="ROR",
            )

    def test_ror_includes_initiator_in_delta(self, tmp_path):
        monomer_results = _make_mol_results(_MONOMER_SMILES, h=-100000.0, s=-50.0, g=-115000.0)
        initiator_results = _make_mol_results(_INITIATOR_SMILES, h=-50000.0, s=-30.0, g=-59000.0)
        n = 1
        polymer_results = [(n, _make_mol_results("P", h=-160000.0, s=-85.0, g=-180000.0))]

        compile_results(
            monomer_smiles=_MONOMER_SMILES,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=initiator_results,
            output_dir=tmp_path,
            reaction_type="ROR",
        )

        import pandas as pd
        from larest.output import slugify

        summary_dir = tmp_path / "Monomer" / slugify(_MONOMER_SMILES) / "summary"
        df = pd.read_csv(summary_dir / "rdkit.csv")
        row = df[df["polymer_length"] == n].iloc[0]

        expected_delta_h = (-160000.0 - n * -100000.0 - -50000.0) / n
        assert row["delta_H"] == pytest.approx(expected_delta_h, rel=1e-6)

    def test_multiple_polymer_lengths(self, tmp_path):
        monomer_results = _make_mol_results(_MONOMER_SMILES)
        polymer_results = [
            (2, _make_mol_results("P2")),
            (3, _make_mol_results("P3")),
            (4, _make_mol_results("P4")),
        ]

        compile_results(
            monomer_smiles=_MONOMER_SMILES,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=None,
            output_dir=tmp_path,
            reaction_type="RER",
        )

        import pandas as pd
        from larest.output import slugify

        summary_dir = tmp_path / "Monomer" / slugify(_MONOMER_SMILES) / "summary"
        df = pd.read_csv(summary_dir / "rdkit.csv")
        assert len(df) == 3
        assert set(df["polymer_length"].astype(int).tolist()) == {2, 3, 4}


# ---------------------------------------------------------------------------
# run_pipeline (mocked)
# ---------------------------------------------------------------------------


class TestRunPipelineMocked:
    def _make_config(self):
        return {
            "steps": {
                "rdkit": True,
                "crest_confgen": False,
                "censo": False,
                "crest_entropy": False,
                "xtb": False,
            },
            "rdkit": {
                "n_conformers": 1,
                "random_seed": 42,
                "conformer_box_size": 2.0,
                "mmff": "MMFF94",
                "mmff_iters": 10,
                "align_iters": 5,
                "precision": 6,
                "n_cores": 1,
            },
            "xtb": {"parallel": 1, "gfn": 2, "ohess": "vtight", "alpb": "toluene", "etemp": 298.15},
            "crest": {"confgen": {}, "entropy": {}},
            "censo": {"cli": {}, "general": {"temperature": 298.15}},
            "reaction": {"type": "RER", "lengths": [2], "initiator": "", "monomers": [_MONOMER_SMILES]},
        }

    def _mol_dir(self, tmp_path, mol):
        """Return and pre-create the molecule output directory.

        run_rdkit normally creates this via create_dir(rdkit_dir, parents=True).
        When run_rdkit is mocked we must create it ourselves so results.json
        can be written.
        """
        from larest.output import create_dir, slugify

        match mol:
            case Polymer():
                d = tmp_path / "Polymer" / f"{slugify(mol.monomer_smiles)}_{mol.length}"
            case _:
                d = tmp_path / type(mol).__name__ / slugify(mol.smiles)
        create_dir(d)
        return d

    def test_run_pipeline_monomer(self, tmp_path):
        config = self._make_config()
        fake_rdkit = {"rdkit": {"H": -100000.0, "S": -50.0, "G": -115000.0}}
        mol = Monomer(_MONOMER_SMILES)
        self._mol_dir(tmp_path, mol)

        with patch("larest.main.run_rdkit", return_value=fake_rdkit):
            result = run_pipeline(mol, tmp_path, config)

        assert isinstance(result, MolResults)
        assert result.smiles == _MONOMER_SMILES
        assert result.sections["rdkit"]["H"] == pytest.approx(-100000.0)

    def test_run_pipeline_writes_results_json(self, tmp_path):
        config = self._make_config()
        fake_rdkit = {"rdkit": {"H": -100000.0, "S": -50.0, "G": -115000.0}}
        mol = Monomer(_MONOMER_SMILES)
        mol_dir = self._mol_dir(tmp_path, mol)

        with patch("larest.main.run_rdkit", return_value=fake_rdkit):
            run_pipeline(mol, tmp_path, config)

        assert (mol_dir / "results.json").exists()

    def test_run_pipeline_polymer_dir_name(self, tmp_path):
        config = self._make_config()
        fake_rdkit = {"rdkit": {"H": -100000.0, "S": -50.0, "G": -115000.0}}
        polymer = Polymer(smiles="OCC(=O)OCC(=O)O", monomer_smiles=_MONOMER_SMILES, length=2)
        expected_dir = self._mol_dir(tmp_path, polymer)

        with patch("larest.main.run_rdkit", return_value=fake_rdkit):
            run_pipeline(polymer, tmp_path, config)

        assert expected_dir.exists()

    def test_run_pipeline_skips_disabled_steps(self, tmp_path):
        config = self._make_config()
        config["steps"]["rdkit"] = False
        mol = Monomer(_MONOMER_SMILES)
        self._mol_dir(tmp_path, mol)

        with (
            patch("larest.main.run_rdkit") as mock_rdkit,
            patch("larest.main.restore_results", return_value=({}, MagicMock(value=5))),
        ):
            run_pipeline(mol, tmp_path, config)

        mock_rdkit.assert_not_called()
