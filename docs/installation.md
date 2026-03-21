# Installation

## Prerequisites

LaREST requires [conda](https://docs.conda.io) to manage its environment. The `env.yaml` file installs xTB, CREST, and all Python dependencies automatically.

**ORCA** must be installed separately — it is not available via conda. Download it from the [ORCA forum](https://orcaforum.kofo.mpg.de) and ensure the `orca` binary is on your `PATH`.

## Steps

```bash
git clone https://github.com/Ryan-Reese/LaREST.git
cd LaREST
conda env create -f env.yaml
conda activate larest
```

This installs xTB, CREST, and the `larest` Python package (plus all Python dependencies) into the `larest` conda environment.

## Post-installation

After installing ORCA, update the `orcaversion` field in your config file to match your installed version:

```toml
[censo.paths]
orcaversion = "6.1.0"   # must match your ORCA installation exactly
```

## Tested versions

| Dependency | Version |
|---|---|
| xTB | 6.7.1 |
| CREST | 3.0.2 |
| CENSO | 2.1.4 |
| ORCA | 6.1.0 |
