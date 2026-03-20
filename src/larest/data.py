"""Dataclasses representing molecules and their computed results in the LaREST pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Monomer:
    """A lactone monomer molecule.

    Attributes
    ----------
    smiles : str
        SMILES string of the monomer.
    """

    smiles: str


@dataclass
class Initiator:
    """An initiating alcohol molecule for ring-opening polymerization.

    Attributes
    ----------
    smiles : str
        SMILES string of the initiator.
    """

    smiles: str


@dataclass
class Polymer:
    """A polymer chain at a given repeat unit length.

    Attributes
    ----------
    smiles : str
        SMILES string of the polymer chain.
    monomer_smiles : str
        SMILES string of the monomer used to build the chain.
    length : int
        Number of monomer repeat units in the chain.
    """

    smiles: str
    monomer_smiles: str
    length: int


@dataclass
class MolResults:
    """Thermodynamic results for a single molecule across pipeline sections.

    Attributes
    ----------
    smiles : str
        SMILES string of the molecule.
    sections : dict[str, dict[str, float | None]]
        Mapping of pipeline section name (e.g. ``"rdkit"``, ``"3_REFINEMENT"``)
        to a dict of thermodynamic parameters (``"H"``, ``"S"``, ``"G"``).
        Values are ``None`` where a parameter could not be computed.
    """

    smiles: str
    sections: dict[str, dict[str, float | None]] = field(default_factory=dict)
