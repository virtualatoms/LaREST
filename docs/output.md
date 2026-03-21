# Output

All results are written under the directory specified by `-o` (default: `./output`).

## Directory structure

```
output/
  Monomer/<smiles_slug>/
    rdkit/
      conformers.sdf              # all MMFF conformers
      xtb/rdkit/results.csv       # xTB H, S, G for each conformer
    crest_confgen/
      crest_best.xyz              # lowest-energy CREST conformer
      crest_conformers.xyz        # full CREST ensemble
    censo/
      censo.txt                   # raw CENSO output
      results.json                # parsed H, S, G for each sub-stage
      3_REFINEMENT.xyz            # refined conformer ensemble
    crest_entropy/
      crest.txt                   # raw CREST entropy output
      results.json                # S_conf, S_rrho, S_total
    results.json                  # merged results across all sections
    summary/
      rdkit.csv                   # ΔH, ΔS, ΔG vs polymer length (RDKit level)
      crest.csv                   # ΔH, ΔS, ΔG vs polymer length (CREST level)
      censo_refinement.csv        # ΔH, ΔS, ΔG vs polymer length (DFT level)
      censo_corrected.csv         # ΔH, ΔS, ΔG with CREST entropy correction

  Polymer/<smiles_slug>_<length>/
    ...                           # same structure as Monomer/

  Initiator/<smiles_slug>/
    ...                           # same structure as Monomer/

  larest.log                      # full pipeline log
```

## Summary CSVs

The `summary/` CSVs are the primary output. Each file corresponds to one pipeline section and contains one row per polymer length with the following columns:

| Column | Description |
|---|---|
| `polymer_length` | Chain length *n* |
| `monomer_H/S/G` | Monomer thermodynamic parameters (J/mol) |
| `initiator_H/S/G` | Initiator parameters (J/mol); zero for RER |
| `polymer_H/S/G` | Polymer parameters (J/mol) |
| `delta_H` | ΔH of ring-opening (J/mol) |
| `delta_S` | ΔS of ring-opening (J/mol/K) |
| `delta_G` | ΔG of ring-opening (J/mol) |

Delta values are computed as:

```
delta_param = (polymer_param - n * monomer_param - initiator_param) / n
```

## Sections

| Section | Description |
|---|---|
| `rdkit` | RDKit + xTB level |
| `crest` | CREST + xTB level |
| `censo_prescreening` | CENSO prescreening sub-stage |
| `censo_screening` | CENSO screening sub-stage |
| `censo_optimization` | CENSO optimisation sub-stage |
| `censo_refinement` | CENSO refinement sub-stage (highest accuracy) |
| `censo_corrected` | `censo_refinement` + CREST conformational entropy |

## Log file

`larest.log` is written to the output directory and contains the full pipeline log at DEBUG level. Console output respects the `--verbose` flag.
