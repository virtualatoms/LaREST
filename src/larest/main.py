from argparse import Namespace
from logging import Logger
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from tqdm import tqdm

from larest.chem import build_polymer
from larest.constants import THERMODYNAMIC_PARAMS
from larest.data import Initiator, MolResults, Monomer, Polymer
from larest.output import create_dir, slugify
from larest.parsers import LarestArgumentParser
from larest.pipeline import MolPipeline
from larest.setup import get_config, get_logger


def compile_results(
    monomer_smiles: str,
    monomer_results: MolResults,
    polymer_results: list[tuple[int, MolResults]],
    initiator_results: MolResults | None,
    output_dir: Path,
    reaction_type: str,
) -> None:
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
                    initiator_results.sections[section][param]
                    if reaction_type == "ROR" and initiator_results is not None
                    else 0,
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


def main(
    args: Namespace,
    config: dict[str, Any],
    logger: Logger,
) -> None:
    output_dir = Path(args.output)
    reaction_type = config["reaction"]["type"]

    logger.info("Running pipeline for monomers")
    for monomer_smiles in tqdm(
        config["reaction"]["monomers"],
        desc="Running pipeline for monomers",
    ):
        monomer = Monomer(smiles=monomer_smiles)

        try:
            monomer_results = MolPipeline(monomer, output_dir, config).run()
        except Exception:
            logger.exception(f"Failed to run pipeline for monomer {monomer_smiles}")
            continue

        initiator_results = None
        if reaction_type == "ROR":
            initiator = Initiator(smiles=config["reaction"]["initiator"])
            try:
                initiator_results = MolPipeline(initiator, output_dir, config).run()
            except Exception:
                logger.exception(
                    f"Failed to run pipeline for initiator {config['reaction']['initiator']}",
                )
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
                logger.exception(
                    f"Failed to build polymer for monomer {monomer_smiles} (length: {length})",
                )
                continue
            polymers.append(
                Polymer(smiles=polymer_smiles, monomer_smiles=monomer_smiles, length=length)
            )

        polymer_results: list[tuple[int, MolResults]] = []
        for polymer in tqdm(polymers, desc="Running pipeline for each polymer length"):
            try:
                result = MolPipeline(polymer, output_dir, config).run()
                polymer_results.append((polymer.length, result))
            except Exception:
                logger.exception(
                    f"Failed to run pipeline for polymer {monomer_smiles} (length: {polymer.length})",
                )
                continue

        compile_results(
            monomer_smiles=monomer_smiles,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=initiator_results,
            output_dir=output_dir,
            reaction_type=reaction_type,
        )


def entry_point() -> None:
    """Entry point of the LaREST package script"""
    parser: LarestArgumentParser = LarestArgumentParser()
    args: Namespace = parser.parse_args()

    try:
        config: dict[str, Any] = get_config(args=args)
    except Exception as err:
        raise SystemExit(1) from err

    try:
        logger: Logger = get_logger(
            name=__name__,
            args=args,
            config=config,
        )
    except Exception as err:
        raise SystemExit(1) from err
    else:
        logger.info("LaREST Initialised")
        for config_key, config_value in config.items():
            logger.debug(f"{config_key} config:\n{config_value}")

    main(args=args, config=config, logger=logger)
