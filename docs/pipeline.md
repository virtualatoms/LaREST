# Pipeline architecture

LaREST chains four computational stages to produce thermodynamic parameters (H, S, G) for lactone ring-opening polymerization reactions. Each stage feeds into the next, progressively increasing the level of theory.

## Stages

### 1. RDKit + xTB

**Purpose:** Generate a set of initial molecular conformers and rank them by free energy.

1. RDKit generates MMFF conformers (default: 50 per molecule).
2. Each conformer is optimised at the xTB GFN2 level and a Hessian is computed to obtain H, S, and G.
3. The lowest-free-energy conformer is passed to the next stage.

Controlled by `[rdkit]` and `[xtb]` in `config.toml`. Set `steps.rdkit = false` to skip.

### 2. CREST conformer ensemble

**Purpose:** Explore the conformational space more thoroughly using metadynamics.

1. CREST's iMTD-GC algorithm generates a conformer/rotamer ensemble from the best RDKit conformer.
2. The ensemble is deduplicated and sorted by energy (CREGEN).
3. xTB re-ranks the ensemble by free energy; the lowest-energy conformer is passed forward.

Controlled by `[crest.confgen]` in `config.toml`. Set `steps.crest_confgen = false` to skip.

### 3. CENSO + ORCA (DFT refinement)

**Purpose:** Refine the CREST ensemble with density functional theory using four sub-stages of increasing accuracy.

| Sub-stage | Label | Default functional | Default basis |
|---|---|---|---|
| Prescreening | `censo_prescreening` | PBE-D4 | def2-SV(P) |
| Screening | `censo_screening` | r2SCAN-3c | def2-TZVP |
| Optimisation | `censo_optimization` | r2SCAN-3c | def2-TZVP |
| Refinement | `censo_refinement` | wB97X-V | def2-TZVP |

Each sub-stage applies an energy window threshold to prune the ensemble before passing it to the next stage. ORCA is used as the QM backend throughout.

Controlled by `[censo.*]` in `config.toml`. Set `steps.censo = false` to skip.

### 4. CREST entropy

**Purpose:** Compute the conformational entropy correction using CREST's entropy mode.

1. CREST re-explores the conformational space using GFN-FF (fast) to obtain a well-converged entropy estimate.
2. The resulting S_conf is added to the CENSO refinement results to produce the `censo_corrected` section.

Controlled by `[crest.entropy]` in `config.toml`. Set `steps.crest_entropy = false` to skip.

## Checkpointing

At the start of each molecule's run, LaREST walks through the output directory and identifies the first missing result file. All stages up to that point are skipped; execution resumes from there. This means interrupted runs can be restarted without any manual intervention.

## Molecules processed

For each monomer SMILES in the config, LaREST runs the full pipeline for:

1. The **monomer** itself.
2. The **initiator** (ROR only).
3. Each **polymer** at every requested chain length.

Polymer SMILES are constructed automatically from the monomer and chain length. Final reaction thermodynamics (ΔH, ΔS, ΔG) are computed in `compile_results` and written to `summary/` CSVs.
