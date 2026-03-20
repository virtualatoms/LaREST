from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path
from typing import Any

from larest.constants import CENSO_SECTIONS, HARTTREE_TO_JMOL, THERMODYNAMIC_PARAMS
from larest.output import create_dir
from larest.setup import parse_command_args
from larest.setup import create_censorc

logger = logging.getLogger(__name__)


def run_censo(
    dir_path: Path,
    output_dir: Path,
    config: dict[str, Any],
) -> dict[str, dict[str, float | None]]:
    temp_dir = output_dir / "temp"
    create_censorc(config=config, temp_dir=temp_dir)

    censo_dir = dir_path / "censo"
    create_dir(censo_dir)

    censo_output_file = censo_dir / "censo.txt"
    censo_config_file = temp_dir / ".censo2rc"
    crest_conformers_file = dir_path / "crest_confgen" / "crest_conformers.xyz"

    censo_args: list[str] = [
        "censo",
        "--input",
        str(crest_conformers_file.absolute()),
        "--inprc",
        str(censo_config_file.absolute()),
    ] + parse_command_args(sub_config=["censo", "cli"], config=config)

    with open(censo_output_file, "w") as fstream:
        subprocess.run(
            censo_args,
            stdout=fstream,
            stderr=subprocess.STDOUT,
            cwd=censo_dir,
            check=True,
        )

    censo_results: dict[str, dict[str, float | None]] = parse_censo_output(
        censo_output_file=censo_output_file,
        temperature=config["censo"]["general"]["temperature"],
    )

    censo_results_file = censo_dir / "results.json"
    logger.debug(f"Writing results to {censo_results_file}")
    with open(censo_results_file, "w") as fstream:
        json.dump(censo_results, fstream, sort_keys=True, allow_nan=True, indent=4)

    return censo_results


def parse_censo_output(
    censo_output_file: Path,
    temperature: float,
) -> dict[str, dict[str, float | None]]:
    # WARN: cannot use fromkeys, otherwise they all point to the same mutable dict
    censo_output: dict[str, dict[str, float | None]] = {
        section: dict.fromkeys(THERMODYNAMIC_PARAMS, None) for section in CENSO_SECTIONS
    }

    logger.debug(f"Searching for CENSO results in file {censo_output_file}")
    with open(censo_output_file) as fstream:
        section_no: int = 0
        for i, line in enumerate(fstream):
            if f"part{section_no}" in line:
                try:
                    censo_output[CENSO_SECTIONS[section_no]]["H"] = (
                        float(line.split()[1]) * HARTTREE_TO_JMOL
                    )
                    censo_output[CENSO_SECTIONS[section_no]]["G"] = (
                        float(line.split()[2]) * HARTTREE_TO_JMOL
                    )
                except Exception:
                    logger.exception(
                        f"Failed to extract H and G from line {i}: {line}",
                    )
                else:
                    section_no += 1

    for params in censo_output.values():
        if params["H"] is not None and params["G"] is not None:
            params["S"] = (params["H"] - params["G"]) / temperature
        if not all(v is not None for v in params.values()):
            logger.warning(
                f"Failed to extract necessary data from {censo_output_file}",
            )
            logger.warning("Missing data will be assigned None")
        logger.debug(
            f"Found enthalpy: {params['H']}, free energy: {params['G']}, entropy {params['S']}",
        )

    return censo_output


def parse_best_censo_conformers(
    censo_output_file: Path,
) -> dict[str, str]:
    best_censo_conformers: dict[str, str] = dict.fromkeys(CENSO_SECTIONS, "CONF0")

    logger.debug(f"Searching for results in file {censo_output_file}")
    section_no: int = 0
    with open(censo_output_file) as fstream:
        for i, line in enumerate(fstream):
            if "Highest ranked conformer" in line:
                try:
                    best_censo_conformers[CENSO_SECTIONS[section_no]] = (
                        line.split()[-1]
                    )
                except Exception:
                    logger.exception(
                        f"Failed to extract best conformer from line {i}: {line}",
                    )
                else:
                    section_no += 1

    if not all(best_censo_conformers.values()):
        logger.warning(
            f"Failed to extract best conformers from {censo_output_file}",
        )
        logger.warning("Missing sections will be assigned first conformer (CONF0)")
    for section in CENSO_SECTIONS:
        logger.debug(
            f"Best conformer in {section}: {best_censo_conformers[section]}",
        )

    return best_censo_conformers


def extract_best_conformer_xyz(
    censo_conformers_xyz_file: Path,
    best_conformer_id: str,
    output_xyz_file: Path,
) -> None:
    logger.debug(
        f"Extracting best conformer ({best_conformer_id}) .xyz from {censo_conformers_xyz_file}",
    )
    with open(censo_conformers_xyz_file) as fin:
        conformers_xyz: list[str] = fin.readlines()
        n_atoms: int = int(conformers_xyz[0])
        for i, line in enumerate(conformers_xyz):
            if best_conformer_id in line.split():
                with open(output_xyz_file, "w") as fout:
                    fout.writelines(conformers_xyz[i - 1 : i + n_atoms + 1])
                break

    logger.debug(f"Finished extracting best conformer xyz to {output_xyz_file}")
