"""Tests for larest.data dataclasses."""

from larest.data import Initiator, MolResults, Monomer, Polymer


class TestMonomer:
    def test_creation(self):
        m = Monomer(smiles="C1CC(=O)O1")
        assert m.smiles == "C1CC(=O)O1"

    def test_equality(self):
        assert Monomer("C1CC(=O)O1") == Monomer("C1CC(=O)O1")
        assert Monomer("C1CC(=O)O1") != Monomer("C1CCCC(=O)O1")


class TestInitiator:
    def test_creation(self):
        i = Initiator(smiles="CCO")
        assert i.smiles == "CCO"


class TestPolymer:
    def test_creation(self):
        p = Polymer(smiles="OCC(=O)OCC(=O)O", monomer_smiles="C1CC(=O)O1", length=2)
        assert p.smiles == "OCC(=O)OCC(=O)O"
        assert p.monomer_smiles == "C1CC(=O)O1"
        assert p.length == 2

    def test_equality(self):
        p1 = Polymer("ABC", "C1CC(=O)O1", 2)
        p2 = Polymer("ABC", "C1CC(=O)O1", 2)
        assert p1 == p2

    def test_inequality_length(self):
        p1 = Polymer("ABC", "C1CC(=O)O1", 2)
        p2 = Polymer("ABC", "C1CC(=O)O1", 3)
        assert p1 != p2


class TestMolResults:
    def test_creation_empty(self):
        r = MolResults(smiles="C1CC(=O)O1")
        assert r.smiles == "C1CC(=O)O1"
        assert r.sections == {}

    def test_creation_with_sections(self):
        sections: dict[str, dict[str, float | None]] = {
            "rdkit": {"H": 1.0, "S": 2.0, "G": 3.0},
        }
        r = MolResults(smiles="C1CC(=O)O1", sections=sections)
        assert r.sections["rdkit"]["H"] == 1.0

    def test_sections_default_is_independent(self):
        """Each instance should have its own dict."""
        r1 = MolResults("A")
        r2 = MolResults("B")
        r1.sections["rdkit"] = {"H": 1.0}
        assert "rdkit" not in r2.sections

    def test_none_values_allowed(self):
        r = MolResults("C", sections={"rdkit": {"H": None, "S": None, "G": None}})
        assert r.sections["rdkit"]["H"] is None
