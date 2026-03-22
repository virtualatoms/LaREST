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

## Test data

Sample fixture files live in `tests/data/` and include xTB output, CREST entropy output, an RDKit results CSV, CENSO output, and a multi-conformer XYZ file. These are used by the unit tests so no external tools need to be running.

## Continuous integration

Tests run automatically on every push and pull request to `main` via GitHub Actions (see `.github/workflows/test.yml`). The workflow sets up the conda environment from `env.yaml` and runs the full integration test suite.
