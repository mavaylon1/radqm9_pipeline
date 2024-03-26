from typing import List, Optional
from hashlib import blake2b

from pydantic import Field

from pymatgen.core.structure import Molecule

from emmet.core.math import Vector3D
from emmet.core.mpid import MPculeID
from emmet.core.material import PropertyOrigin
from emmet.core.qchem.task import TaskDocument
from emmet.core.molecules.molecule_property import PropertyDoc


__author__ = "Evan Spotte-Smith <ewcspottesmith@lbl.gov>"


class ForcePointDoc(PropertyDoc):
    property_name: str = "force_single_point"

    # Required properties
    species: List[str] = Field(
        ...,
        description="Atomic element symbols",
    )

    coordinates: List[List[float]] = Field(
        ...,
        description="Cartesian coordinates of the atoms at this point (units: Angstrom)"
    )

    energy: float = Field(
        ...,
        description="Electronic energy at this point (units: Ha)"
    )    

    forces: List[List[float]] = Field(
        ...,
        description="Atomic forces at this point (units: Ha/Bohr)"
    )

    mulliken_partial_charges: List[float] = Field(
        ..., description="Atomic partial charges for this point, using the Mulliken method"
    )

    mulliken_partial_spins: List[float] = Field(
        ..., description="Atomic partial spins for this point, using the Mulliken method"
    )

    resp_partial_charges: List[float] = Field(
        ..., description="Atomic partial charges for this point, using the restrained electrostatic potential "
                         "(RESP) method"
    )

    dipole_moment: Vector3D = Field(
        ...,
        description="Molecular dipole moment vector at this point (units: Debye)",
    )

    # Optional properties
    nbo_partial_charges: Optional[List[float]] = Field(
        None, description="Atomic partial charges for this point, using the natural bonding orbital (NBO) method"
    )

    nbo_partial_spins: Optional[List[float]] = Field(
        None, description="Atomic partial spins for this point, using the natural bonding orbital (NBO) method"
    )

    precise_forces: Optional[List[List[float]]] = Field(
        None,
        description="High-precision atomic forces at this point (units: Ha/Bohr)"
    )

    pcm_forces: Optional[List[List[float]]] = Field(
        None,
        description="Electrostatic atomic forces at this point from polarizable continuum model (PCM) implicit "
                    "solvation (units: Ha/Bohr)."
    )

    cds_forces: Optional[List[List[float]]] = Field(
        None,
        description="Atomic force contributions at this point from cavitation, dispersion, and structural "
                    "rearrangement in the SMx family of implicit solvent models (units: Ha/Bohr)"
    )

    resp_dipole_moment: Optional[Vector3D] = Field(
        None,
        description="Molecular dipole moment vector at this point, calculated via restrained electrostatic potential "
                    "(RESP) (units: Debye)",
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
        Construct a force document from a task document

        :param task: document from which force properties can be extracted
        :param molecule_id: MPculeID
        :param deprecated: bool. Is this document deprecated?
        :param kwargs: to pass to PropertyDoc
        :return:
        """

        if task.task_type.value != "Force":
            raise ValueError("ForcePointDoc can only be constructed from force calculations,"
                             f"not {task.task_type.value}!")

        mol = task.output.initial_molecule

        # Basic information
        species = [str(i) for i in mol.species]
        coordinates = mol.cart_coords

        # Energy and forces
        energy = task.output.final_energy
        forces = task.output.gradients
        precise_forces = task.output.precise_gradients
        pcm_forces = task.output.pcm_gradients
        cds_forces = task.output.CDS_gradients

        calc = task.calcs_reversed[0]

        # Precise forces are either in output or don't exist
        # For PCM and CDS forces, can check "calcs_reversed"
        if pcm_forces is None:
            pcm_forces = calc.get("pcm_gradients")
        
        if cds_forces is None:
            cds_forces = calc.get("CDS_gradients")

        # Partial charges and spins
        
        # Mulliken
        if isinstance(task.output.mulliken[0], list):
            # Open-shell
            mulliken_partial_charges = [float(mull[0]) for mull in task.output.mulliken]
            mulliken_partial_spins = [float(mull[1]) for mull in task.output.mulliken]
        else:
            # Closed-shell - no partial spins
            mulliken_partial_charges = [float(i) for i in task.output.mulliken]
            mulliken_partial_spins = [0.0 for i in task.output.mulliken]

        # RESP
        resp_partial_charges = [float(i) for i in task.output.resp]

        # NBO
        if task.output.nbo is None:
            nbo_partial_charges = None
            nbo_partial_spins = None
        else:
            nbo_partial_charges = [
                float(task.output.nbo["natural_populations"][0]["Charge"][str(i)])
                for i in range(len(mol))
            ]
            nbo_partial_spins = [
                float(task.output.nbo["natural_populations"][0]["Density"][str(i)])
                for i in range(len(mol))
            ]

        # Dipoles
        # For SP and force calcs, dipoles stored in output
        dipoles = task.output.dipoles
        if dipoles is None:
            dipole_moment = None
            resp_dipole_moment = None
        else:
            dipole_moment = dipoles.get("dipole")
            resp_dipole_moment = dipoles.get("RESP_dipole")

        # If needed, look for dipoles in calcs_reversed
        calcs_reversed = task.calcs_reversed

        calc = calcs_reversed[0]
        grab_index = -1

        dipoles = calc.get("dipoles", dict())

        if dipole_moment is None:
            dipole_moment = dipoles.get("dipole")
            if isinstance(dipole_moment, list) and len(dipole_moment) == 0:
                dipole_moment = None
            elif isinstance(dipole_moment[0], list):
                dipole_moment = dipole_moment[grab_index]
        if resp_dipole_moment is None:
            resp_dipole_moment = dipoles.get("RESP_dipole")
            if isinstance(resp_dipole_moment, list) and len(resp_dipole_moment) == 0:
                resp_dipole_moment = None
            elif isinstance(resp_dipole_moment[0], list):
                resp_dipole_moment = resp_dipole_moment[grab_index]

        # Important fields missing
        if any(
            [
                x is None
                for x in [
                    energy,
                    forces,
                    mulliken_partial_charges,
                    mulliken_partial_spins,
                    resp_partial_charges,
                    dipole_moment
                ]
            ]
        ):
            return None

        id_string = f"force_single_point-{molecule_id}-{task.task_id}-{task.lot_solvent}"
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
            species=species,
            coordinates=coordinates,
            energy=energy,
            forces=forces,
            mulliken_partial_charges=mulliken_partial_charges,
            mulliken_partial_spins=mulliken_partial_spins,
            resp_partial_charges=resp_partial_charges,
            dipole_moment=dipole_moment,
            nbo_partial_charges=nbo_partial_charges,
            nbo_partial_spins=nbo_partial_spins,
            precise_forces=precise_forces,
            pcm_forces=pcm_forces,
            cds_forces=cds_forces,
            resp_dipole_moment=resp_dipole_moment,
            origins=[PropertyOrigin(name="forces", task_id=task.task_id)],
            deprecated=deprecated,
            **kwargs,
        )
