from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path
from typing import Any

from larest.constants import HARTTREE_TO_JMOL, THERMODYNAMIC_PARAMS
from larest.setup import parse_command_args

logger = logging.getLogger(__name__)


def run_xtb(xtb_input_file: Path, xtb_dir: Path, config: dict[str, Any]) -> dict[str, float | None]:
    xtb_args: list[str] = [
        "xtb",
        str(xtb_input_file.absolute()),
        "--namespace",
        xtb_input_file.name.split(".")[0],
    ] + parse_command_args(sub_config=["xtb"], config=config)

    xtb_output_file = xtb_dir / "xtb.txt"
    with open(xtb_output_file, "w") as fstream:
        subprocess.run(
            xtb_args,
            stdout=fstream,
            stderr=subprocess.STDOUT,
            cwd=xtb_dir,
            check=True,
        )

    xtb_results: dict[str, float | None] = parse_xtb_output(
        xtb_output_file=xtb_output_file,
        temperature=config["xtb"]["etemp"],
    )

    xtb_results_file = xtb_dir / "results.json"
    logger.debug(f"Writing results to {xtb_results_file}")
    with open(xtb_results_file, "w") as fstream:
        json.dump(xtb_results, fstream, sort_keys=True, allow_nan=True, indent=4)

    return xtb_results


def parse_xtb_output(
    xtb_output_file: Path,
    temperature: float,
) -> dict[str, float | None]:
    xtb_output: dict[str, float | None] = dict.fromkeys(THERMODYNAMIC_PARAMS, None)

    logger.debug(f"Searching for xTB results in file {xtb_output_file}")
    with open(xtb_output_file) as fstream:
        for i, line in enumerate(fstream):
            if "TOTAL ENTHALPY" in line:
                try:
                    xtb_output["H"] = float(line.split()[3]) * HARTTREE_TO_JMOL
                except Exception:
                    logger.exception(
                        f"Failed to extract total enthalpy from line {i}: {line}",
                    )
            elif "TOTAL FREE ENERGY" in line:
                try:
                    xtb_output["G"] = float(line.split()[4]) * HARTTREE_TO_JMOL
                except Exception:
                    logger.exception(
                        f"Failed to extract total free energy from line {i}: {line}",
                    )

    if xtb_output["H"] is not None and xtb_output["G"] is not None:
        xtb_output["S"] = (xtb_output["H"] - xtb_output["G"]) / temperature
    if not all(v is not None for v in xtb_output.values()):
        logger.warning(f"Failed to extract necessary data from {xtb_output_file}")
        logger.warning("Missing data will be assigned None")
    logger.debug(
        f"Found enthalpy: {xtb_output['H']}, entropy: {xtb_output['S']}, free energy {xtb_output['G']}",
    )

    return xtb_output
