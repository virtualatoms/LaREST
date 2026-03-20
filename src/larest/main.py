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
            results |= apply_entropy_correction(results["3_REFINEMENT"], crest_entropy_results)

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
        for polymer in tqdm(polymers, desc="Running pipeline for each polymer length"):
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
