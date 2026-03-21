# LaREST

**LaREST** (Lactone Ring-opening Energetics Sorting Tool) is a computational chemistry pipeline for predicting the thermodynamics of lactone ring-opening polymerization reactions.

Given a set of monomer SMILES strings and polymer chain lengths, it computes Boltzmann-averaged enthalpies (H), entropies (S), and free energies (G) by chaining together four levels of theory:

| Stage | Tool | Description |
|---|---|---|
| 1 | RDKit + xTB | MMFF conformer generation, ranked by xTB free energy |
| 2 | CREST | Conformer/rotamer ensemble exploration (iMTD-GC) |
| 3 | CENSO + ORCA | DFT refinement of the ensemble (four sub-stages) |
| 4 | CREST entropy | Conformational entropy correction |

Two reaction types are supported:

- **RER** (Ring Equilibrium Reaction) — no initiator; polymer chain only
- **ROR** (Ring-Opening polymerization Reaction) — includes an initiating alcohol

Final reaction thermodynamics are computed as:

$$\Delta X = \frac{X_\text{polymer} - n \cdot X_\text{monomer} - X_\text{initiator}}{n}$$

where *n* is the polymer chain length.

```{toctree}
:maxdepth: 2
:caption: Getting started
:hidden:

installation
usage
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

configuration
pipeline
output
testing
```
