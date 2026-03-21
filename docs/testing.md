# Testing

The test suite lives in `tests/` and uses [pytest](https://pytest.org).

## Running tests

**Unit tests only** (no external tools required):

```bash
pytest tests/
```

**Including integration tests** (requires `xtb` and `crest` on `PATH`):

```bash
pytest tests/ --integration
```

Integration tests are marked `@pytest.mark.integration` and skipped by default. CENSO/ORCA tests are not included, as they require a full DFT installation.

## Test modules

| Module | Coverage |
|---|---|
| `test_data.py` | `Monomer`, `Polymer`, `Initiator`, `MolResults` dataclasses |
| `test_output.py` | `slugify`, `create_dir`, `remove_dir` |
| `test_setup.py` | Config loading, deep-merge, parallelisation propagation, `parse_command_args` |
| `test_chem.py` | SMILES parsing, ring-size detection, polymer construction (RER and ROR) |
| `test_xtb.py` | xTB output parsing and `run_xtb` (mocked + integration) |
| `test_crest.py` | CREST entropy output parsing, `run_crest_confgen`, `run_crest_entropy` (mocked + integration) |
| `test_censo.py` | CENSO output parsing, conformer extraction, `.censo2rc` creation (no ORCA required) |
| `test_checkpoint.py` | `PipelineStage` ordering, `restore_results`, `apply_entropy_correction` |
| `test_rdkit_stage.py` | `parse_best_rdkit_conformer`, `run_rdkit` (mocked + integration) |
| `test_main.py` | `compile_results` delta calculations, `run_pipeline` orchestration |

## Test data

Sample fixture files live in `tests/data/` and include xTB output, CREST entropy output, an RDKit results CSV, CENSO output, and a multi-conformer XYZ file. These are used by the unit tests so no external tools need to be running.

## Continuous integration

Tests run automatically on every push and pull request to `main` via GitHub Actions (see `.github/workflows/test.yml`). The workflow sets up the conda environment from `env.yaml` and runs the full integration test suite.
