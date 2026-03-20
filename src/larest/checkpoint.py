import json
import logging
from collections.abc import Callable
from enum import IntEnum
from pathlib import Path

from larest.constants import PIPELINE_SECTIONS, THERMODYNAMIC_PARAMS
from larest.parsers import parse_best_rdkit_conformer

logger = logging.getLogger(__name__)


class PipelineStage(IntEnum):
    RDKIT = 1
    CREST_CONFGEN = 2
    CENSO = 3
    CREST_ENTROPY = 4
    FINISH = 5


def restore_results(
    dir_path: Path,
) -> tuple[dict[str, dict[str, float | None]], PipelineStage]:
    results: dict[str, dict[str, float | None]] = {
        section: dict.fromkeys(THERMODYNAMIC_PARAMS, None)
        for section in PIPELINE_SECTIONS
    }

    if not _load_stage(
        path=dir_path / "xtb" / "rdkit" / "results.csv",
        load_fn=lambda path: _parse_rdkit(results, path),
        label="rdkit",
    ):
        return results, PipelineStage.RDKIT

    if not _load_stage(
        path=dir_path / "xtb" / "crest" / "results.json",
        load_fn=lambda path: results.update({"crest": json.loads(path.read_text())}),
        label="crest_confgen",
    ):
        return results, PipelineStage.CREST_CONFGEN

    if not _load_stage(
        path=dir_path / "censo" / "results.json",
        load_fn=lambda path: results.update(json.loads(path.read_text())),
        label="censo",
    ):
        return results, PipelineStage.CENSO

    if not _load_stage(
        path=dir_path / "crest_entropy" / "results.json",
        load_fn=lambda path: _parse_crest_entropy(results, path),
        label="crest_entropy",
    ):
        return results, PipelineStage.CREST_ENTROPY

    _load_stage(
        path=dir_path / "results.json",
        load_fn=lambda path: results.update(json.loads(path.read_text())),
        label="final",
    )
    return results, PipelineStage.FINISH


def _load_stage(path: Path, load_fn: Callable[[Path], None], label: str) -> bool:
    """Attempt to load a checkpoint stage. Returns True if loaded, False if missing or failed."""
    if not path.exists():
        logger.info(f"No pre-existing {label} results at {path}, running pipeline...")
        return False
    try:
        load_fn(path)
    except Exception:
        logger.exception(f"Failed to load pre-existing {label} results from {path}")
        return False
    logger.info(f"Loaded pre-existing {label} results from {path}")
    return True


def _parse_rdkit(results: dict, path: Path) -> None:
    best = parse_best_rdkit_conformer(path)
    del best["conformer_id"]
    results.update({"rdkit": best})


def _parse_crest_entropy(results: dict, path: Path) -> None:
    crest_entropy_results: dict[str, float | None] = json.loads(path.read_text())
    censo_corrected_results = results["3_REFINEMENT"].copy()
    if censo_corrected_results["S"] is None or crest_entropy_results["S_total"] is None:
        raise ValueError(
            "Failed to apply CREST entropy correction to CENSO results",
        )
    censo_corrected_results["S"] += crest_entropy_results["S_total"]
    results.update({"censo_corrected": censo_corrected_results})
