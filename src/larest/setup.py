"""Configuration loading and logging setup for the LaREST pipeline.

:func:`get_config` reads the user's ``config.toml`` directly into a plain
dict with no defaults merging — the file must be complete.

:func:`get_logger` applies the ``[logging]`` section of that dict via
:func:`logging.config.dictConfig`, substituting the resolved path for the
log file handler.

:func:`parse_command_args` converts a nested sub-config dict into a flat
list of CLI flags suitable for passing directly to ``subprocess.run``.
"""

from __future__ import annotations

import copy
import logging
import logging.config
import tomllib
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def get_config(config_dir: Path) -> dict[str, Any]:
    """Load the pipeline configuration from ``config.toml``.

    Parameters
    ----------
    config_dir : Path
        Directory containing ``config.toml``.

    Returns
    -------
    dict[str, Any]
        Parsed TOML configuration as a nested dict.

    Raises
    ------
    Exception
        Re-raises any exception raised by :func:`tomllib.load` after printing
        a human-readable error message.
    """
    config_file = config_dir / "config.toml"
    try:
        with open(config_file, "rb") as fstream:
            return tomllib.load(fstream)
    except Exception:
        print(f"Failed to load config from {config_file}")
        raise


def get_logger(output_dir: Path, config: dict[str, Any]) -> None:
    """Configure the Python logging system from the ``[logging]`` config section.

    Deep-copies the logging sub-config before patching the file handler's
    ``filename`` with an absolute path inside *output_dir*, then applies it
    with :func:`logging.config.dictConfig`.

    Parameters
    ----------
    output_dir : Path
        Root output directory; the log file is written here using the filename
        specified in ``config["logging"]["handlers"]["file"]["filename"]``.
    config : dict[str, Any]
        Full pipeline configuration dict.  Must contain a ``[logging]`` section
        compatible with :func:`logging.config.dictConfig`.

    Raises
    ------
    Exception
        Re-raises any exception raised during logging setup after printing a
        human-readable error message.
    """
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
    """Convert a nested config section into a flat list of CLI flags.

    Traverses *config* using the keys in *sub_config* to locate the target
    section, then converts each key-value pair to ``--key value`` tokens.
    Boolean values are handled specially: ``True`` emits ``--key`` with no
    value; ``False`` is omitted entirely.

    Parameters
    ----------
    sub_config : list[str]
        Sequence of keys used to navigate to the desired section, e.g.
        ``["xtb"]`` or ``["crest", "confgen"]``.
    config : dict[str, Any]
        Full pipeline configuration dict.

    Returns
    -------
    list[str]
        Flat list of CLI argument tokens ready to be appended to a
        ``subprocess.run`` command list.  Returns an empty list if the
        sub-config path cannot be found.
    """
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

