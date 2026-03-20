import argparse
import logging
from typing import Any

logger = logging.getLogger(__name__)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="LaREST")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="./output",
        help="Output directory",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default="./config",
        help="Config directory",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="increase output verbosity",
        action="store_true",
    )
    return parser


def parse_command_args(
    sub_config: list[str],
    config: dict[str, Any],
) -> list[str]:
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
