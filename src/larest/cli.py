from __future__ import annotations

import argparse
import logging
from pathlib import Path

from larest.main import main
from larest.setup import get_config, get_logger

logger = logging.getLogger(__name__)


def entry_point() -> None:
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
    args = parser.parse_args()

    output_dir = Path(args.output)
    config_dir = Path(args.config)

    try:
        config = get_config(config_dir)
    except Exception as err:
        raise SystemExit(1) from err

    try:
        get_logger(output_dir=output_dir, config=config)
    except Exception as err:
        raise SystemExit(1) from err

    logger.info("LaREST Initialised")
    for config_key, config_value in config.items():
        logger.debug(f"{config_key} config:\n{config_value}")

    main(output_dir=output_dir, config=config)
