"""Checkpoint and result-restoration utilities for the LaREST pipeline.

At the start of each molecule's run, :func:`restore_results` walks through the
expected output files in pipeline order and returns both the accumulated results
and the first stage that has not yet been completed.  This allows interrupted
runs to resume from exactly where they left off without re-executing earlier
stages.

:func:`apply_entropy_correction` combines CENSO ``censo_refinement`` enthalpy /
free-energy values with the conformational entropy obtained from a CREST
entropy run, producing the ``censo_corrected`` section stored alongside the
other thermodynamic sections.
"""

import json
import logging
from collections.abc import Callable
from enum import IntEnum
from pathlib import Path

from larest.constants import PIPELINE_SECTIONS, THERMODYNAMIC_PARAMS
from larest.rdkit import parse_best_rdkit_conformer

logger = logging.getLogger(__name__)


class PipelineStage(IntEnum):
    """Ordered enumeration of pipeline stages used for checkpointing.

    Stages are assigned integer values so that simple comparisons (``<=``)
    can determine whether a given stage still needs to be executed.
    """

    RDKIT = 1
    CREST_CONFGEN = 2
    CENSO = 3
    CREST_ENTROPY = 4
    FINISH = 5


def restore_results(
    dir_path: Path,
) -> tuple[dict[str, dict[str, float | None]], PipelineStage]:
    """Load previously completed pipeline results and determine the next stage to run.

    Inspects the expected output files for each pipeline stage in order.  The
    first missing or unreadable file determines the stage that must be (re-)run.
    All successfully loaded results are merged and returned alongside that stage.

    Parameters
    ----------
    dir_path : Path
        Root output directory for the molecule (e.g.
        ``output/Monomer/<slug>``).

    Returns
    -------
    results : dict[str, dict[str, float | None]]
        Mapping of section name to thermodynamic parameters populated from any
        checkpoint files that were found.  Parameters that have not yet been
        computed remain ``None``.
    stage : PipelineStage
        The earliest pipeline stage whose results were not found.  The caller
        should (re-)run this stage and all subsequent ones.
    """
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
    """Attempt to load a single checkpoint stage from disk.

    Parameters
    ----------
    path : Path
        Expected location of the checkpoint file.
    load_fn : Callable[[Path], None]
        Function that reads *path* and mutates the shared results dict in
        the enclosing scope.
    label : str
        Human-readable name of the stage, used only for log messages.

    Returns
    -------
    bool
        ``True`` if the file was found and loaded without error, ``False``
        if the file was absent or ``load_fn`` raised an exception.
    """
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
    """Load the best RDKit/xTB conformer results into *results*.

    Parameters
    ----------
    results : dict
        Shared results dict to update in-place with the ``"rdkit"`` section.
    path : Path
        Path to the xTB results CSV produced during the RDKit stage.
    """
    best = parse_best_rdkit_conformer(path)
    del best["conformer_id"]
    results.update({"rdkit": best})


def apply_entropy_correction(
    refinement_results: dict[str, float | None],
    crest_entropy_results: dict[str, float | None],
) -> dict[str, dict[str, float | None]]:
    """Apply a CREST conformational entropy correction to CENSO refinement results.

    Combines the enthalpy and free energy from the CENSO ``censo_refinement`` stage
    with the total conformational entropy from a CREST entropy calculation to
    produce a corrected thermodynamic section (``"censo_corrected"``).

    Parameters
    ----------
    refinement_results : dict[str, float | None]
        Thermodynamic parameters from the ``censo_refinement`` CENSO stage, keyed
        by ``"H"``, ``"S"``, and ``"G"``.
    crest_entropy_results : dict[str, float | None]
        Entropy parameters from the CREST entropy run, keyed by
        ``"S_conf"``, ``"S_rrho"``, and ``"S_total"``.

    Returns
    -------
    dict[str, dict[str, float | None]]
        A single-entry dict ``{"censo_corrected": {...}}`` where the ``"S"``
        value has been replaced by the CREST total conformational entropy.

    Raises
    ------
    ValueError
        If either ``refinement_results["S"]`` or
        ``crest_entropy_results["S_total"]`` is ``None``.
    """
    if refinement_results["S"] is None or crest_entropy_results["S_total"] is None:
        raise ValueError("Failed to apply CREST entropy correction to CENSO results")
    corrected = refinement_results.copy()
    corrected["S"] += crest_entropy_results["S_total"]
    return {"censo_corrected": corrected}


def _parse_crest_entropy(results: dict, path: Path) -> None:
    """Load CREST entropy results and apply the entropy correction to *results*.

    Parameters
    ----------
    results : dict
        Shared results dict; must already contain a ``"censo_refinement"`` section.
        Updated in-place with a ``"censo_corrected"`` section.
    path : Path
        Path to the CREST entropy ``results.json`` file.
    """
    crest_entropy_results: dict[str, float | None] = json.loads(path.read_text())
    results.update(apply_entropy_correction(results["censo_refinement"], crest_entropy_results))
