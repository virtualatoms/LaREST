from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Monomer:
    smiles: str


@dataclass
class Initiator:
    smiles: str


@dataclass
class Polymer:
    smiles: str
    monomer_smiles: str
    length: int


@dataclass
class MolResults:
    smiles: str
    sections: dict[str, dict[str, float | None]] = field(default_factory=dict)
