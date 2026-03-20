"""File-system utilities for directory management and SMILES-to-path conversion.

Provides helpers for creating and removing output directories, and a
:func:`slugify` function that converts arbitrary SMILES strings into safe,
filesystem-compatible directory names.
"""

import logging
import re
import shutil
import unicodedata
from pathlib import Path

logger = logging.getLogger(__name__)


def create_dir(dir_path: Path) -> None:
    """Create a directory, including any missing parent directories.

    Equivalent to ``mkdir -p``; no error is raised if the directory already
    exists.

    Parameters
    ----------
    dir_path : Path
        Path of the directory to create.
    """
    logger.debug(f"Creating directory: {dir_path}")
    dir_path.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Directory {dir_path} created")


def remove_dir(dir_path: Path) -> None:
    """Remove a directory tree, logging a warning on failure instead of raising.

    Parameters
    ----------
    dir_path : Path
        Root of the directory tree to remove.
    """
    logger.debug(f"Removing directory: {dir_path}")
    try:
        shutil.rmtree(dir_path, ignore_errors=False)
    except Exception:
        logger.warning(f"Failed to remove directory {dir_path}")
    else:
        logger.debug(f"Directory {dir_path} removed")


# modified from django.utils
def slugify(smiles: str) -> str:
    """Convert a SMILES string into a safe filesystem directory name.

    Applies Unicode normalisation, strips non-ASCII characters, replaces
    parentheses with hyphens, removes remaining special characters (except
    ``@`` and ``-``), and collapses consecutive hyphens/spaces into a single
    hyphen.

    Parameters
    ----------
    smiles : str
        SMILES string to convert.

    Returns
    -------
    str
        Sanitised slug suitable for use as a directory name.
    """
    smiles = (
        unicodedata.normalize("NFKD", smiles).encode("ascii", "ignore").decode("ascii")
    )
    smiles = re.sub(r"[\(\)]", "-", smiles)
    smiles = re.sub(r"[^\w\s@-]", "", smiles)
    return re.sub(r"[-\s]+", "-", smiles).strip("-_")
