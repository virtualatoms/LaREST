from __future__ import annotations

import json
import logging
from functools import cached_property
from pathlib import Path
from typing import Any

from larest.censo import run_censo
from larest.checkpoint import PipelineStage, apply_entropy_correction, restore_results
from larest.crest import run_crest_confgen, run_crest_entropy
from larest.data import Initiator, MolResults, Monomer, Polymer
from larest.output import create_dir, remove_dir, slugify
from larest.rdkit import run_rdkit

logger = logging.getLogger(__name__)


class MolPipeline:
    def __init__(
        self,
        mol: Monomer | Polymer | Initiator,
        output_dir: Path,
        config: dict[str, Any],
    ) -> None:
        self.mol = mol
        self.output_dir = output_dir
        self.config = config

    @cached_property
    def _dir_path(self) -> Path:
        match self.mol:
            case Polymer():
                return (
                    self.output_dir
                    / "Polymer"
                    / f"{slugify(self.mol.monomer_smiles)}_{self.mol.length}"
                )
            case _:
                return (
                    self.output_dir
                    / type(self.mol).__name__
                    / slugify(self.mol.smiles)
                )

    def run(self) -> MolResults:
        """Run the LaREST pipeline for the molecule."""
        try:
            self._setup()
            results, stage = restore_results(self._dir_path)

            if self.config["steps"]["rdkit"] and stage <= PipelineStage.RDKIT:
                results |= run_rdkit(self.mol.smiles, self._dir_path, self.config)
            if self.config["steps"]["crest_confgen"] and stage <= PipelineStage.CREST_CONFGEN:
                results |= run_crest_confgen(self._dir_path, self.config)
            if self.config["steps"]["censo"] and stage <= PipelineStage.CENSO:
                results |= run_censo(self._dir_path, self.output_dir, self.config)
            if self.config["steps"]["crest_entropy"] and stage <= PipelineStage.CREST_ENTROPY:
                crest_entropy_results = run_crest_entropy(self._dir_path, self.config)
                if self.config["steps"]["censo"]:
                    results |= apply_entropy_correction(results["3_REFINEMENT"], crest_entropy_results)

            self._write_final_results(results)
            self._cleanup()
        except Exception:
            logger.exception("Error encountered within pipeline, exiting...")
            raise

        return MolResults(smiles=self.mol.smiles, sections=results)

    def _setup(self) -> None:
        create_dir(self.output_dir / "temp")

    def _cleanup(self) -> None:
        remove_dir(self.output_dir / "temp")

    def _write_final_results(self, results: dict[str, dict[str, float | None]]) -> None:
        results_file = self._dir_path / "results.json"
        logger.debug(f"Writing final results to {results_file}")
        with open(results_file, "w") as fstream:
            json.dump(results, fstream, sort_keys=True, indent=4, allow_nan=True)
