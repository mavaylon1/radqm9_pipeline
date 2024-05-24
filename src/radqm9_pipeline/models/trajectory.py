from typing import List, Optional
from hashlib import blake2b

import numpy as np

from pydantic import Field

from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Molecule

from emmet.core.mpid import MPculeID
from emmet.core.material import PropertyOrigin
from emmet.core.qchem.task import TaskDocument
from emmet.core.molecules.molecule_property import PropertyDoc


__author__ = "Evan Spotte-Smith <ewcspottesmith@lbl.gov>"


class TrajectoryDoc(PropertyDoc):
    property_name: str = "optimization_trajectory"

    num_trajectories: int = Field(
        ...,
        description="Number of separate optimization trajectories extracted from this task"
    )

    species: List[str] = Field(
        ...,
        description="Element symbols for each atom in the molecule"
    )

    geometries: List[List[List[List[float]]]] = Field(
        ...,
        description="XYZ positions of each atom in the molecule for each optimization step for each optimization "
                    "trajectory (units: Angstrom)",
    )

    energies: List[List[float]] = Field(
        ...,
        description="Electronic energies for each optimization step for each optimization trajectory (units: Hartree)"
    )

    forces: List[List[List[List[float]]]] = Field(
        ...,
        description="Forces on each atom for each optimization step for each optimization trajectory (units: Ha/Bohr)"
    )

    pcm_forces: Optional[List[Optional[List[List[List[float]]]]]] = Field(
        None,
        description="Electrostatic atomic forces from polarizable continuum model (PCM) implicit solvation "
                    "for each optimization step for each optimization trajectory (units: Ha/Bohr)."
    )

    cds_forces: Optional[List[Optional[List[List[List[float]]]]]] = Field(
        None,
        description="Atomic force contributions from cavitation, dispersion, and structural rearrangement in the SMx "
                    "family of implicit solvent models, for each optimization step for each optimization trajectory "
                    "(units: Ha/Bohr)"
    )

    mulliken_partial_charges: Optional[List[Optional[List[List[float]]]]] = Field(
        None,
        description="Partial charges of each atom for each optimization step for each optimization trajectory, using "
                    "the Mulliken method"
    )

    mulliken_partial_spins: Optional[List[Optional[List[List[float]]]]] = Field(
        None,
        description="Partial spins of each atom for each optimization step for each optimization trajectory, using "
                    "the Mulliken method"
    )

    resp_partial_charges: Optional[List[Optional[List[List[float]]]]] = Field(
        None,
        description="Partial charges of each atom for each optimization step for each optimization trajectory, using "
                    "the restrained electrostatic potential (RESP) method"
    )

    dipole_moments: Optional[List[Optional[List[List[float]]]]] = Field(
        None,
        description="Molecular dipole moment for each optimization step for each optimization trajectory, "
                    "(units: Debye)"
    )

    resp_dipole_moments: Optional[List[Optional[List[List[float]]]]] = Field(
        None,
        description="Molecular dipole moment for each optimization step for each optimization trajectory, "
                    "using the restrainted electrostatic potential (RESP) method (units: Debye)"
    )


    @classmethod
    def from_task(
        cls,
        task: TaskDocument,
        molecule_id: MPculeID,
        deprecated: bool = False,
        **kwargs,
    ):  # type: ignore[override]
        """
        Construct a trajectory document from a task document

        :param task: document from which force properties can be extracted
        :param molecule_id: MPculeID
        :param deprecated: bool. Is this document deprecated?
        :param kwargs: to pass to PropertyDoc
        :return:
        """

        if task.task_type.value not in [
            "Geometry Optimization",
            "Frequency Flattening Geometry Optimization",
        ]:
            raise ValueError("TrajectoryDoc can only be constructed from geometry optimization calculations,"
                             f"not {task.task_type.value}!")

        if task.output.optimized_molecule is not None:
            mol = task.output.optimized_molecule
        else:
            mol = task.output.initial_molecule

        species = None
        geometries = list()
        energies = list()
        total_gradients = list()
        pcm_gradients = list()
        cds_gradients = list()
        mulliken_partial_charges = list()
        mulliken_partial_spins = list()
        resp_partial_charges = list()
        dipole_moments = list()
        resp_dipole_moments = list()

        for calculation in task.calcs_reversed:
            species = calculation.get("species", species)
            this_geometries = calculation.get("geometries")
            this_energies = calculation.get("energy_trajectory")
            this_total_gradients = calculation.get("gradients")
            this_pcm_gradients = calculation.get("pcm_gradients")
            this_cds_gradients = calculation.get("CDS_gradients")
            
            this_mulliken = calculation.get("Mulliken")
            this_resp = calculation.get("RESP")
            this_dipoles = calculation.get("dipoles")

            valid_trajectory = True
            if this_geometries is None or this_energies is None:
                # No valid geometry optimization found
                valid_trajectory = False
            elif len(this_energies) != len(this_total_gradients):
                # Energies and forces cannot be trivially mapped
                valid_trajectory = False
            elif len(this_geometries) != len(this_energies):
                # Initial geometry not included - common because of how parsing is done
                if len(this_geometries) == len(this_energies) - 1:
                    this_geometries = [calculation["initial_geometry"]] + this_geometries
                # Other issue - no one-to-one mapping of molecule structure and energy
                else:
                    valid_trajectory = False

            if not valid_trajectory:
                continue

            if isinstance(calculation["initial_molecule"], Molecule):
                init_mol = calculation["initial_molecule"]  # type: ignore
            else:
                init_mol = Molecule.from_dict(calculation["initial_molecule"])  # type: ignore

            if species is None:
                species = init_mol.species  # type: ignore

            # Number of steps
            # All data in this (sub)-trajectory must have the same length
            num_steps = len(this_geometries)

            this_dipole_moments = None
            this_resp_dipole_moments = None

            # electric dipoles
            if this_dipoles is not None:
                if this_dipoles.get("dipole") is not None and len(this_dipoles['dipole']) > 0:
                    if (
                        isinstance(this_dipoles["dipole"][0], list)
                        and len(this_dipoles["dipole"]) == num_steps
                    ):
                        this_dipole_moments = this_dipoles["dipole"]
                if (
                    this_dipoles.get("RESP_dipole") is not None
                    and len(this_dipoles["RESP_dipole"]) > 0
                ):
                    if (
                        isinstance(this_dipoles["RESP_dipole"][0], list)
                        and len(this_dipoles["RESP_dipole"]) == num_steps
                    ):
                        this_resp_dipole_moments = this_dipoles["RESP_dipole"]

            this_mulliken_partial_charges = None
            this_mulliken_partial_spins = None
            this_resp_partial_charges = None

            # Partial charges/spins
            if this_mulliken is not None:
                if len(this_mulliken) == num_steps:
                    if int(mol.spin_multiplicity) == 1:
                        this_mulliken_partial_charges = this_mulliken
                    else:
                        # For open-shell molecules, need to split mulliken charges and spins
                        charges = list()
                        spins = list()

                        for step in this_mulliken:
                            step_charges = list()
                            step_spins = list()
                            for atom in step:
                                step_charges.append(atom[0])
                                step_spins.append(atom[1])
                            charges.append(step_charges)
                            spins.append(step_spins)
                        
                        this_mulliken_partial_charges = charges
                        this_mulliken_partial_spins = spins
                elif len(this_mulliken) == num_steps + 1:
                    last = np.asarray(this_mulliken[-1])
                    seclast = np.asarray(this_mulliken[-2])
                    if np.allclose(last, seclast):
                        
                        if int(mol.spin_multiplicity) == 1:
                            this_mulliken_partial_charges = this_mulliken[:-1]
                        
                        else:
                            charges = list()
                            spins = list()

                            for step in this_mulliken[::-1]:
                                step_charges = list()
                                step_spins = list()
                                for atom in step:
                                    step_charges.append(atom[0])
                                    step_spins.append(atom[1])
                                charges.append(step_charges)
                                spins.append(step_spins)
                            
                            this_mulliken_partial_charges = charges
                            this_mulliken_partial_spins = spins

            if this_resp is not None:
                if len(this_resp) == num_steps:
                    this_resp_partial_charges = this_resp
                elif len(this_resp) == num_steps + 1:
                    last = np.asarray(this_resp[-1])
                    seclast = np.asarray(this_resp[-2])
                    if np.allclose(last, seclast):
                        this_resp_partial_charges = this_resp[:-1]

            geometries.append(this_geometries)
            energies.append(this_energies)
            total_gradients.append(this_total_gradients)
            pcm_gradients.append(this_pcm_gradients)
            cds_gradients.append(this_cds_gradients)
            mulliken_partial_charges.append(this_mulliken_partial_charges)
            mulliken_partial_spins.append(this_mulliken_partial_spins)
            resp_partial_charges.append(this_resp_partial_charges)
            dipole_moments.append(this_dipole_moments)
            resp_dipole_moments.append(this_resp_dipole_moments)

        num_trajectories = len(geometries)

        id_string = f"trajectory-{molecule_id}-{task.task_id}-{task.lot_solvent}"
        h = blake2b()
        h.update(id_string.encode("utf-8"))
        property_id = h.hexdigest()

        return super().from_molecule(
            meta_molecule=mol,
            property_id=property_id,
            molecule_id=molecule_id,
            level_of_theory=task.level_of_theory,
            solvent=task.solvent,
            lot_solvent=task.lot_solvent,
            num_trajectories=num_trajectories,
            species=species,
            geometries=geometries,
            energies=energies,
            forces=total_gradients,
            pcm_forces=pcm_gradients,
            cds_forces=cds_gradients,
            mulliken_partial_charges=mulliken_partial_charges,
            mulliken_partial_spins=mulliken_partial_spins,
            resp_partial_charges=resp_partial_charges,
            dipole_moments=dipole_moments,
            resp_dipole_moments=resp_dipole_moments,
            origins=[PropertyOrigin(name="trajectory", task_id=task.task_id)],
            deprecated=deprecated,
            **kwargs,
        )
