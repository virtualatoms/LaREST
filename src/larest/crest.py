from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path
from typing import Any

from larest.censo import extract_best_conformer_xyz, parse_best_censo_conformers
from larest.constants import CALMOL_TO_JMOL, CREST_ENTROPY_OUTPUT_PARAMS
from larest.output import create_dir
from larest.parsers import parse_command_args
from larest.rdkit import parse_best_rdkit_conformer
from larest.xtb import run_xtb

logger = logging.getLogger(__name__)


def run_crest_confgen(
    dir_path: Path,
    config: dict[str, Any],
) -> dict[str, dict[str, float | None]]:
    crest_dir = dir_path / "crest_confgen"
    create_dir(crest_dir)

    xtb_rdkit_results_file = dir_path / "xtb" / "rdkit" / "results.csv"
    best_rdkit_conformer_results = parse_best_rdkit_conformer(xtb_rdkit_results_file)
    if best_rdkit_conformer_results["conformer_id"] is None:
        raise ValueError("Failed to extract best RDKit conformer for crest_confgen")

    best_rdkit_conformer_id = int(best_rdkit_conformer_results["conformer_id"])
    best_rdkit_conformer_xyz_file = (
        dir_path
        / "xtb"
        / "rdkit"
        / f"conformer_{best_rdkit_conformer_id}"
        / f"conformer_{best_rdkit_conformer_id}.xtbopt.xyz"
    )

    crest_output_file = crest_dir / "crest.txt"
    crest_args: list[str] = [
        "crest",
        str(best_rdkit_conformer_xyz_file.absolute()),
    ] + parse_command_args(sub_config=["crest", "confgen"], config=config)

    try:
        with open(crest_output_file, "w") as fstream:
            subprocess.run(
                crest_args,
                stdout=fstream,
                stderr=subprocess.STDOUT,
                cwd=crest_dir,
                check=True,
            )
    except Exception:
        logger.exception("Failed to run CREST confgen")
        raise

    if config["steps"]["xtb"]:
        best_crest_conformer_xyz_file = crest_dir / "crest_best.xyz"
        xtb_dir = dir_path / "xtb" / "crest"
        create_dir(xtb_dir)
        xtb_results = run_xtb(
            xtb_input_file=best_crest_conformer_xyz_file,
            xtb_dir=xtb_dir,
            config=config,
        )
        return {"crest": xtb_results}

    return {}


def run_crest_entropy(
    dir_path: Path,
    config: dict[str, Any],
) -> dict[str, float | None]:
    censo_dir = dir_path / "censo"
    best_censo_conformer_xyz_file = censo_dir / "censo_best.xyz"

    try:
        best_censo_conformer = parse_best_censo_conformers(
            censo_output_file=censo_dir / "censo.txt",
        )["3_REFINEMENT"]
        extract_best_conformer_xyz(
            censo_conformers_xyz_file=censo_dir / "3_REFINEMENT.xyz",
            best_conformer_id=best_censo_conformer,
            output_xyz_file=best_censo_conformer_xyz_file,
        )
    except Exception as err:
        raise ValueError(
            "Failed to extract best CENSO conformer for crest_entropy",
        ) from err

    crest_dir = dir_path / "crest_entropy"
    create_dir(crest_dir)

    crest_output_file = crest_dir / "crest.txt"
    crest_args: list[str] = [
        "crest",
        str(best_censo_conformer_xyz_file.absolute()),
    ] + parse_command_args(sub_config=["crest", "entropy"], config=config)

    try:
        with open(crest_output_file, "w") as fstream:
            subprocess.run(
                crest_args,
                stdout=fstream,
                stderr=subprocess.STDOUT,
                cwd=crest_dir,
                check=True,
            )
    except Exception:
        logger.exception("Failed to run CREST entropy")
        raise

    crest_results: dict[str, float | None] = parse_crest_entropy_output(
        crest_output_file=crest_output_file,
    )

    crest_results_file = crest_dir / "results.json"
    logger.debug(f"Writing results to {crest_results_file}")
    with open(crest_results_file, "w") as fstream:
        json.dump(crest_results, fstream, sort_keys=True, allow_nan=True, indent=4)

    return crest_results


def parse_crest_entropy_output(
    crest_output_file: Path,
) -> dict[str, float | None]:
    crest_output: dict[str, float | None] = dict.fromkeys(
        CREST_ENTROPY_OUTPUT_PARAMS,
        None,
    )

    logger.debug(f"Searching for results in file {crest_output_file}")
    try:
        with open(crest_output_file) as fstream:
            for i, line in enumerate(fstream):
                if "Sconf" in line:
                    try:
                        crest_output["S_conf"] = (
                            float(line.split()[-1]) * CALMOL_TO_JMOL
                        )
                    except Exception:
                        logger.exception(
                            f"Failed to extract S_conf from line {i}: {line}",
                        )
                elif ("+" in line) and ("δSrrho" in line) and (len(line.split()) == 4):
                    try:
                        crest_output["S_rrho"] = (
                            float(line.split()[-1]) * CALMOL_TO_JMOL
                        )
                    except Exception:
                        logger.exception(
                            f"Failed to extract S_rrho from line {i}: {line}",
                        )
                elif ("S(total)" in line) and ("cal" in line):
                    try:
                        crest_output["S_total"] = (
                            float(line.split()[3]) * CALMOL_TO_JMOL
                        )
                    except Exception:
                        logger.exception(
                            f"Failed to extract S_total from line {i}: {line}",
                        )
    except Exception:
        logger.exception(
            f"Failed to parse crest entropy results from {crest_output_file}",
        )
        raise

    if not all(value is not None for value in crest_output.values()):
        logger.warning(f"Failed to extract necessary data from {crest_output_file}")
        logger.warning("Missing data will be assigned None")
    logger.debug(
        f"Found S_conf: {crest_output['S_conf']}, S_rrho: {crest_output['S_rrho']}",
    )
    logger.debug(f"S_total {crest_output['S_total']}")

    return crest_output
