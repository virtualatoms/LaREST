from __future__ import annotations

import copy
import logging
import logging.config
import tomllib
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def get_config(config_dir: Path) -> dict[str, Any]:
    config_file = config_dir / "config.toml"
    try:
        with open(config_file, "rb") as fstream:
            return tomllib.load(fstream)
    except Exception:
        print(f"Failed to load config from {config_file}")
        raise


def get_logger(output_dir: Path, config: dict[str, Any]) -> None:
    try:
        log_config: dict[str, Any] = copy.deepcopy(config["logging"])
        log_config["handlers"]["file"]["filename"] = (
            output_dir / log_config["handlers"]["file"]["filename"]
        ).resolve()
        logging.config.dictConfig(log_config)
    except Exception:
        print(f"Failed to setup logging config from {config}")
        raise


def parse_command_args(sub_config: list[str], config: dict[str, Any]) -> list[str]:
    cfg = config
    try:
        for key in sub_config:
            cfg = cfg[key]
    except KeyError:
        logger.warning(f"Failed to find sub-config {sub_config} in config, using no arguments")
        return []

    args = []
    for k, v in cfg.items():
        if v is False:
            continue
        args.append(f"--{k}")
        if v is not True:
            args.append(str(v))
    logger.debug(f"Returning {sub_config} args: {args}")
    return args


def create_censorc(config: dict[str, Any], temp_dir: Path) -> None:
    censorc_file = temp_dir / ".censo2rc"
    censo_config: dict[str, Any] = config["censo"]

    with open(censorc_file, "w") as fstream:
        for header, sub_config in censo_config.items():
            fstream.write(f"[{header}]\n")
            fstream.writelines(f"{key} = {value}\n" for key, value in sub_config.items())
            fstream.write("\n")
    logger.debug(f"Created censo config file at {censorc_file}")
