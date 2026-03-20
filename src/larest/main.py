"""Top-level pipeline orchestration for LaREST.

This module contains the three main functions that drive a full LaREST run:

* :func:`run_pipeline` — executes the four stages (RDKit, CREST confgen,
  CENSO, CREST entropy) for a single molecule and writes ``results.json``.
* :func:`compile_results` — combines monomer, polymer, and (for ROR)
  initiator results into per-section ``delta_H / delta_S / delta_G`` CSV
  summaries inside ``<output>/Monomer/<slug>/summary/``.
* :func:`main` — iterates over all configured monomers and polymer lengths,
  calling the two functions above for each combination.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from tqdm import tqdm

from larest.censo import run_censo
from larest.checkpoint import PipelineStage, apply_entropy_correction, restore_results
from larest.chem import build_polymer
from larest.constants import THERMODYNAMIC_PARAMS
from larest.crest import run_crest_confgen, run_crest_entropy
from larest.data import Initiator, MolResults, Monomer, Polymer
from larest.output import create_dir, remove_dir, slugify
from larest.rdkit import run_rdkit

logger = logging.getLogger(__name__)


def run_pipeline(
    mol: Monomer | Polymer | Initiator,
    output_dir: Path,
    config: dict[str, Any],
) -> MolResults:
    """Run the full four-stage pipeline for a single molecule.

    Restores any existing checkpoint results, then conditionally executes
    each enabled stage (controlled by the ``[steps]`` config section) before
    writing a final ``results.json`` file.

    Parameters
    ----------
    mol : Monomer | Polymer | Initiator
        The molecule to process.  Output is placed under
        ``output_dir/Polymer/<slug>_<length>/`` for polymers or
        ``output_dir/<MolType>/<slug>/`` for monomers and initiators.
    output_dir : Path
        Root output directory for the run.
    config : dict[str, Any]
        Full pipeline configuration dict loaded from ``config.toml``.

    Returns
    -------
    MolResults
        Dataclass containing the molecule's SMILES and a mapping of pipeline
        section names to thermodynamic parameter dicts.
    """
    match mol:
        case Polymer():
            dir_path = output_dir / "Polymer" / f"{slugify(mol.monomer_smiles)}_{mol.length}"
        case _:
            dir_path = output_dir / type(mol).__name__ / slugify(mol.smiles)

    create_dir(output_dir / "temp")
    results, stage = restore_results(dir_path)

    if config["steps"]["rdkit"] and stage <= PipelineStage.RDKIT:
        results |= run_rdkit(mol.smiles, dir_path, config)
    if config["steps"]["crest_confgen"] and stage <= PipelineStage.CREST_CONFGEN:
        results |= run_crest_confgen(dir_path, config)
    if config["steps"]["censo"] and stage <= PipelineStage.CENSO:
        results |= run_censo(dir_path, output_dir, config)
    if config["steps"]["crest_entropy"] and stage <= PipelineStage.CREST_ENTROPY:
        crest_entropy_results = run_crest_entropy(dir_path, config)
        if config["steps"]["censo"]:
            results |= apply_entropy_correction(results["censo_refinement"], crest_entropy_results)

    results_file = dir_path / "results.json"
    logger.debug(f"Writing final results to {results_file}")
    with open(results_file, "w") as fstream:
        json.dump(results, fstream, sort_keys=True, indent=4, allow_nan=True)

    remove_dir(output_dir / "temp")
    return MolResults(smiles=mol.smiles, sections=results)


def compile_results(
    monomer_smiles: str,
    monomer_results: MolResults,
    polymer_results: list[tuple[int, MolResults]],
    initiator_results: MolResults | None,
    output_dir: Path,
    reaction_type: str,
) -> None:
    """Compute reaction thermodynamics and write per-section CSV summaries.

    For each pipeline section present in *monomer_results*, constructs a
    DataFrame containing the raw monomer, initiator, and polymer thermodynamic
    values at every polymer length, then appends columns for
    ``delta_H``, ``delta_S``, and ``delta_G`` computed as::

        delta_param = (polymer_param - n * monomer_param - initiator_param) / n

    where *n* is the polymer length.  One CSV file per section is written to
    ``<output_dir>/Monomer/<slug>/summary/``.

    Parameters
    ----------
    monomer_smiles : str
        SMILES string of the monomer (used to locate the summary directory).
    monomer_results : MolResults
        Pipeline results for the monomer molecule.
    polymer_results : list[tuple[int, MolResults]]
        Pairs of (polymer_length, MolResults) for each requested chain length.
    initiator_results : MolResults | None
        Pipeline results for the initiator.  Required when *reaction_type* is
        ``"ROR"``; pass ``None`` for ``"RER"`` reactions.
    output_dir : Path
        Root output directory for the run.
    reaction_type : str
        Either ``"ROR"`` or ``"RER"``.  Determines whether initiator
        thermodynamics are included in the delta calculation.

    Raises
    ------
    ValueError
        If *reaction_type* is ``"ROR"`` but *initiator_results* is ``None``.
    """
    if reaction_type == "ROR" and initiator_results is None:
        raise ValueError(f"Initiator results are required for ROR reactions but are missing (monomer: {monomer_smiles})")

    summary_dir = output_dir / "Monomer" / slugify(monomer_smiles) / "summary"
    create_dir(summary_dir)

    summary_headings: list[str] = ["polymer_length"]
    for mol_type in ["monomer", "initiator", "polymer"]:
        summary_headings.extend([f"{mol_type}_{param}" for param in THERMODYNAMIC_PARAMS])

    for section in monomer_results.sections:
        summary: dict[str, list[float | None]] = {heading: [] for heading in summary_headings}

        for length, poly_results in polymer_results:
            summary["polymer_length"].append(length)
            for param in THERMODYNAMIC_PARAMS:
                summary[f"polymer_{param}"].append(poly_results.sections[section][param])
                summary[f"monomer_{param}"].append(monomer_results.sections[section][param])
                summary[f"initiator_{param}"].append(
                    initiator_results.sections[section][param] if reaction_type == "ROR" else 0,
                )

        summary_df = pd.DataFrame(
            summary,
            index=None,
            dtype=np.float64,
        ).sort_values("polymer_length", ascending=True)

        for param in THERMODYNAMIC_PARAMS:
            summary_df[f"delta_{param}"] = (
                summary_df[f"polymer_{param}"]
                - (summary_df["polymer_length"] * summary_df[f"monomer_{param}"])
                - summary_df[f"initiator_{param}"]
            ) / summary_df["polymer_length"]

        summary_df.to_csv(summary_dir / f"{section}.csv", header=True, index=False)


def main(output_dir: Path, config: dict[str, Any]) -> None:
    """Run the complete LaREST pipeline for all configured monomers and lengths.

    Iterates over every monomer SMILES in ``config["reaction"]["monomers"]``,
    runs the pipeline for the monomer (and initiator if ROR), builds all
    requested polymer lengths, runs the pipeline for each polymer, then calls
    :func:`compile_results` to produce the summary CSVs.  Failures for
    individual molecules or polymer lengths are logged and skipped so that the
    remaining molecules can still be processed.

    Parameters
    ----------
    output_dir : Path
        Root output directory; sub-directories are created automatically.
    config : dict[str, Any]
        Full pipeline configuration dict loaded from ``config.toml``.
    """
    reaction_type = config["reaction"]["type"]

    logger.info("Running pipeline for monomers")
    for monomer_smiles in tqdm(config["reaction"]["monomers"], desc="Running pipeline for monomers"):
        monomer = Monomer(smiles=monomer_smiles)

        try:
            monomer_results = run_pipeline(monomer, output_dir, config)
        except Exception:
            logger.exception(f"Failed to run pipeline for monomer {monomer_smiles}")
            continue

        initiator_results = None
        if reaction_type == "ROR":
            initiator = Initiator(smiles=config["reaction"]["initiator"])
            try:
                initiator_results = run_pipeline(initiator, output_dir, config)
            except Exception:
                logger.exception(f"Failed to run pipeline for initiator {config['reaction']['initiator']}")
                continue

        polymers = []
        for length in config["reaction"]["lengths"]:
            try:
                polymer_smiles = build_polymer(
                    monomer_smiles=monomer_smiles,
                    polymer_length=length,
                    reaction_type=reaction_type,
                    config=config,
                )
            except Exception:
                logger.exception(f"Failed to build polymer for monomer {monomer_smiles} (length: {length})")
                continue
            polymers.append(Polymer(smiles=polymer_smiles, monomer_smiles=monomer_smiles, length=length))

        polymer_results: list[tuple[int, MolResults]] = []
        for polymer in tqdm(polymers, desc="Running pipeline for each polymer length", leave=False):
            try:
                result = run_pipeline(polymer, output_dir, config)
                polymer_results.append((polymer.length, result))
            except Exception:
                logger.exception(f"Failed to run pipeline for polymer {monomer_smiles} (length: {polymer.length})")
                continue

        compile_results(
            monomer_smiles=monomer_smiles,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=initiator_results,
            output_dir=output_dir,
            reaction_type=reaction_type,
        )
