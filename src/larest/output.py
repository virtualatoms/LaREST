import logging
import re
import shutil
import unicodedata
from pathlib import Path

logger = logging.getLogger(__name__)


def create_dir(dir_path: Path) -> None:
    logger.debug(f"Creating directory: {dir_path}")
    dir_path.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Directory {dir_path} created")


def remove_dir(dir_path: Path) -> None:
    logger.debug(f"Removing directory: {dir_path}")
    try:
        shutil.rmtree(dir_path, ignore_errors=False)
    except Exception:
        logger.warning(f"Failed to remove directory {dir_path}")
    else:
        logger.debug(f"Directory {dir_path} removed")


# modified from django.utils
def slugify(smiles: str) -> str:
    smiles = (
        unicodedata.normalize("NFKD", smiles).encode("ascii", "ignore").decode("ascii")
    )
    smiles = re.sub(r"[\(\)]", "-", smiles)
    smiles = re.sub(r"[^\w\s@-]", "", smiles)
    return re.sub(r"[-\s]+", "-", smiles).strip("-_")
