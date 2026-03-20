from __future__ import annotations

import json
import logging
import subprocess
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from rdkit.Chem.AllChem import (
    MMFFGetMoleculeForceField,
    MMFFGetMoleculeProperties,
)
from rdkit.Chem.MolStandardize.rdMolStandardize import StandardizeSmiles
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMoleculeConfs
from rdkit.Chem.rdMolAlign import AlignMolConformers
from rdkit.Chem.rdmolfiles import (
    ForwardSDMolSupplier,
    MolToXYZFile,
    SDWriter,
)
from rdkit.Chem.rdmolops import AddHs

from larest.checkpoint import PipelineStage, apply_entropy_correction, restore_results
from larest.chem import get_mol
from larest.constants import KCALMOL_TO_JMOL, THERMODYNAMIC_PARAMS
from larest.data import Initiator, MolResults, Monomer, Polymer
from larest.output import create_dir, remove_dir, slugify
from larest.parsers import (
    extract_best_conformer_xyz,
    parse_best_censo_conformers,
    parse_best_rdkit_conformer,
    parse_censo_output,
    parse_command_args,
    parse_crest_entropy_output,
    parse_xtb_output,
)
from larest.setup import create_censorc

if TYPE_CHECKING:
    from rdkit.Chem.rdchem import Mol
    from rdkit.ForceField.rdForceField import MMFFMolProperties


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
        self._logger = logging.getLogger(type(mol).__name__)

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
                results |= self._run_rdkit()
            if self.config["steps"]["crest_confgen"] and stage <= PipelineStage.CREST_CONFGEN:
                results |= self._run_crest_confgen()
            if self.config["steps"]["censo"] and stage <= PipelineStage.CENSO:
                results |= self._run_censo()
            if self.config["steps"]["crest_entropy"] and stage <= PipelineStage.CREST_ENTROPY:
                crest_entropy_results = self._run_crest_entropy()
                if self.config["steps"]["censo"]:
                    results |= apply_entropy_correction(results["3_REFINEMENT"], crest_entropy_results)

            self._write_final_results(results)
            self._cleanup()
        except Exception:
            self._logger.exception("Error encountered within pipeline, exiting...")
            raise

        return MolResults(smiles=self.mol.smiles, sections=results)

    def _setup(self) -> None:
        create_dir(self.output_dir / "temp")

    def _cleanup(self) -> None:
        remove_dir(self.output_dir / "temp")

    def _run_rdkit(self) -> dict[str, dict[str, float | None]]:
        dir_path = self._dir_path
        rdkit_dir = dir_path / "rdkit"
        create_dir(rdkit_dir)

        self._logger.debug("Generating conformers and computing energies using RDKit")
        rdkit_mol: Mol = AddHs(get_mol(StandardizeSmiles(self.mol.smiles)))

        n_conformers: int = self.config["rdkit"]["n_conformers"]
        self._logger.debug(f"Generating {n_conformers} conformers")
        conformer_ids: list[int] = EmbedMultipleConfs(
            rdkit_mol,
            n_conformers,
            useRandomCoords=True,
            randomSeed=self.config["rdkit"]["random_seed"],
            boxSizeMult=self.config["rdkit"]["conformer_box_size"],
            numThreads=self.config["rdkit"]["n_cores"],
        )

        self._logger.debug(f"Optimising the {n_conformers} conformers")
        MMFFOptimizeMoleculeConfs(
            rdkit_mol,
            numThreads=self.config["rdkit"]["n_cores"],
            maxIters=self.config["rdkit"]["mmff_iters"],
            mmffVariant=self.config["rdkit"]["mmff"],
        )

        self._logger.debug("Computing molecular properties for MMFF")
        mp: MMFFMolProperties = MMFFGetMoleculeProperties(
            rdkit_mol,
            mmffVariant=self.config["rdkit"]["mmff"],
            mmffVerbosity=0,
        )

        self._logger.debug(f"Computing energies for the {n_conformers} conformers")
        conformer_energies: list[tuple[int, float]] = sorted(
            [
                (
                    conformer_id,
                    MMFFGetMoleculeForceField(rdkit_mol, mp, confId=conformer_id).CalcEnergy()
                    * KCALMOL_TO_JMOL,
                )
                for conformer_id in conformer_ids
            ],
            key=lambda x: x[1],
        )

        self._logger.debug(f"Aligning the {n_conformers} conformers by their geometries")
        AlignMolConformers(rdkit_mol, maxIters=self.config["rdkit"]["align_iters"])

        sdf_file = rdkit_dir / "conformers.sdf"
        self._logger.debug(f"Writing conformers and their energies to {sdf_file}")
        try:
            with open(sdf_file, "w") as fstream:
                writer: SDWriter = SDWriter(fstream)
                for cid, energy in conformer_energies:
                    rdkit_mol.SetIntProp("conformer_id", cid)
                    rdkit_mol.SetDoubleProp("energy", energy)
                    writer.write(rdkit_mol, confId=cid)
                writer.close()
        except Exception:
            self._logger.exception(f"Failed to write RDKit conformers to {sdf_file}")
            raise

        self._logger.debug("Computing thermodynamic parameters of conformers using xTB")
        xtb_default_args: list[str] = parse_command_args(sub_config=["xtb"], config=self.config)

        xtb_base_dir = dir_path / "xtb" / "rdkit"
        xtb_results: dict[str, list[float | None]] = {
            "conformer_id": [],
            **{param: [] for param in THERMODYNAMIC_PARAMS},
        }

        self._logger.debug(f"Getting conformer coordinates from {sdf_file}")
        with open(sdf_file, "rb") as sdfstream:
            mol_supplier: ForwardSDMolSupplier = ForwardSDMolSupplier(
                fileobj=sdfstream,
                sanitize=False,
                removeHs=False,
            )
            for conformer in mol_supplier:
                conformer_id: int = conformer.GetIntProp("conformer_id")
                conformer_xyz_file = rdkit_dir / f"conformer_{conformer_id}.xyz"

                try:
                    MolToXYZFile(
                        mol=rdkit_mol,
                        filename=str(conformer_xyz_file),
                        confId=conformer_id,
                        precision=self.config["rdkit"]["precision"],
                    )
                except Exception:
                    self._logger.exception(
                        f"Failed to write conformer coordinates to {conformer_xyz_file}",
                    )
                    raise

                xtb_dir = xtb_base_dir / f"conformer_{conformer_id}"
                create_dir(xtb_dir)

                xtb_output_file = xtb_dir / f"conformer_{conformer_id}.txt"
                xtb_args: list[str] = [
                    "xtb",
                    str(conformer_xyz_file.absolute()),
                    "--namespace",
                    f"conformer_{conformer_id}",
                ] + xtb_default_args

                try:
                    with open(xtb_output_file, "w") as fstream:
                        subprocess.run(
                            xtb_args,
                            stdout=fstream,
                            stderr=subprocess.STDOUT,
                            cwd=xtb_dir,
                            check=True,
                        )
                except Exception:
                    self._logger.exception(f"Failed to run xTB for conformer {conformer_id}")
                    continue

                try:
                    xtb_output = parse_xtb_output(
                        xtb_output_file=xtb_output_file,
                        temperature=self.config["xtb"]["etemp"],
                    )
                except Exception:
                    self._logger.exception(
                        f"Failed to parse xTB results for conformer {conformer_id}",
                    )
                    continue

                xtb_results["conformer_id"].append(conformer_id)
                for param in THERMODYNAMIC_PARAMS:
                    xtb_results[param].append(xtb_output[param])

        xtb_results_file = xtb_base_dir / "results.csv"
        self._logger.debug(f"Writing results to {xtb_results_file}")

        xtb_results_df = pd.DataFrame(xtb_results, dtype=np.float64).sort_values("G")
        xtb_results_df.to_csv(xtb_results_file, header=True, index=False)

        best_rdkit_conformer_results = xtb_results_df.iloc[0].to_dict()
        del best_rdkit_conformer_results["conformer_id"]
        return {"rdkit": best_rdkit_conformer_results}

    def _run_crest_confgen(self) -> dict[str, dict[str, float | None]]:
        dir_path = self._dir_path
        crest_dir = dir_path / "crest_confgen"
        create_dir(crest_dir)

        xtb_rdkit_results_file = dir_path / "xtb" / "rdkit" / "results.csv"
        best_rdkit_conformer_results = parse_best_rdkit_conformer(xtb_rdkit_results_file)
        if best_rdkit_conformer_results["conformer_id"] is None:
            raise ValueError("Failed to extract best RDKit conformer for crest_confgen")

        best_rdkit_conformer_id = int(best_rdkit_conformer_results["conformer_id"])
        best_rdkit_conformer_xyz_file = (
            dir_path
            / "xtb"
            / "rdkit"
            / f"conformer_{best_rdkit_conformer_id}"
            / f"conformer_{best_rdkit_conformer_id}.xtbopt.xyz"
        )

        crest_output_file = crest_dir / "crest.txt"
        crest_args: list[str] = [
            "crest",
            str(best_rdkit_conformer_xyz_file.absolute()),
        ] + parse_command_args(sub_config=["crest", "confgen"], config=self.config)

        try:
            with open(crest_output_file, "w") as fstream:
                subprocess.run(
                    crest_args,
                    stdout=fstream,
                    stderr=subprocess.STDOUT,
                    cwd=crest_dir,
                    check=True,
                )
        except Exception:
            self._logger.exception("Failed to run CREST confgen")
            raise

        if self.config["steps"]["xtb"]:
            best_crest_conformer_xyz_file = crest_dir / "crest_best.xyz"
            xtb_dir = dir_path / "xtb" / "crest"
            create_dir(xtb_dir)
            xtb_results = self._run_xtb(
                xtb_input_file=best_crest_conformer_xyz_file,
                xtb_dir=xtb_dir,
            )
            return {"crest": xtb_results}

        return {}

    def _run_censo(self) -> dict[str, dict[str, float | None]]:
        dir_path = self._dir_path
        temp_dir = self.output_dir / "temp"
        create_censorc(config=self.config, temp_dir=temp_dir)

        censo_dir = dir_path / "censo"
        create_dir(censo_dir)

        censo_output_file = censo_dir / "censo.txt"
        censo_config_file = temp_dir / ".censo2rc"
        crest_conformers_file = dir_path / "crest_confgen" / "crest_conformers.xyz"

        censo_args: list[str] = [
            "censo",
            "--input",
            str(crest_conformers_file.absolute()),
            "--inprc",
            str(censo_config_file.absolute()),
        ] + parse_command_args(sub_config=["censo", "cli"], config=self.config)

        try:
            with open(censo_output_file, "w") as fstream:
                subprocess.run(
                    censo_args,
                    stdout=fstream,
                    stderr=subprocess.STDOUT,
                    cwd=censo_dir,
                    check=True,
                )
        except Exception:
            self._logger.exception("Failed to run CENSO")
            raise

        censo_results: dict[str, dict[str, float | None]] = parse_censo_output(
            censo_output_file=censo_output_file,
            temperature=self.config["censo"]["general"]["temperature"],
        )

        censo_results_file = censo_dir / "results.json"
        self._logger.debug(f"Writing results to {censo_results_file}")
        with open(censo_results_file, "w") as fstream:
            json.dump(censo_results, fstream, sort_keys=True, allow_nan=True, indent=4)

        return censo_results

    def _run_xtb(self, xtb_input_file: Path, xtb_dir: Path) -> dict[str, float | None]:
        xtb_args: list[str] = [
            "xtb",
            str(xtb_input_file.absolute()),
            "--namespace",
            xtb_input_file.name.split(".")[0],
        ] + parse_command_args(sub_config=["xtb"], config=self.config)

        xtb_output_file = xtb_dir / "xtb.txt"
        try:
            with open(xtb_output_file, "w") as fstream:
                subprocess.run(
                    xtb_args,
                    stdout=fstream,
                    stderr=subprocess.STDOUT,
                    cwd=xtb_dir,
                    check=True,
                )
        except Exception:
            self._logger.exception("Failed to run xTB")
            raise

        try:
            xtb_results: dict[str, float | None] = parse_xtb_output(
                xtb_output_file=xtb_output_file,
                temperature=self.config["xtb"]["etemp"],
            )
        except Exception:
            self._logger.exception(f"Failed to parse xtb results in file {xtb_output_file}")
            raise

        xtb_results_file = xtb_dir / "results.json"
        self._logger.debug(f"Writing results to {xtb_results_file}")
        with open(xtb_results_file, "w") as fstream:
            json.dump(xtb_results, fstream, sort_keys=True, allow_nan=True, indent=4)

        return xtb_results

    def _run_crest_entropy(self) -> dict[str, float | None]:
        dir_path = self._dir_path
        censo_dir = dir_path / "censo"
        best_censo_conformer_xyz_file = censo_dir / "censo_best.xyz"

        try:
            best_censo_conformer = parse_best_censo_conformers(
                censo_output_file=censo_dir / "censo.txt",
            )["3_REFINEMENT"]
            extract_best_conformer_xyz(
                censo_conformers_xyz_file=censo_dir / "3_REFINEMENT.xyz",
                best_conformer_id=best_censo_conformer,
                output_xyz_file=best_censo_conformer_xyz_file,
            )
        except Exception as err:
            raise ValueError(
                "Failed to extract best CENSO conformer for crest_entropy",
            ) from err

        crest_dir = dir_path / "crest_entropy"
        create_dir(crest_dir)

        crest_output_file = crest_dir / "crest.txt"
        crest_args: list[str] = [
            "crest",
            str(best_censo_conformer_xyz_file.absolute()),
        ] + parse_command_args(sub_config=["crest", "entropy"], config=self.config)

        try:
            with open(crest_output_file, "w") as fstream:
                subprocess.run(
                    crest_args,
                    stdout=fstream,
                    stderr=subprocess.STDOUT,
                    cwd=crest_dir,
                    check=True,
                )
        except Exception:
            self._logger.exception("Failed to run CREST entropy")
            raise

        crest_results: dict[str, float | None] = parse_crest_entropy_output(
            crest_output_file=crest_output_file,
        )

        crest_results_file = crest_dir / "results.json"
        self._logger.debug(f"Writing results to {crest_results_file}")
        with open(crest_results_file, "w") as fstream:
            json.dump(crest_results, fstream, sort_keys=True, allow_nan=True, indent=4)

        return crest_results

    def _write_final_results(self, results: dict[str, dict[str, float | None]]) -> None:
        results_file = self._dir_path / "results.json"
        self._logger.debug(f"Writing final results to {results_file}")
        with open(results_file, "w") as fstream:
            json.dump(results, fstream, sort_keys=True, indent=4, allow_nan=True)
