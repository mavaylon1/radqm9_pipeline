from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, TypeVar
from hashlib import blake2b

from pydantic import Field
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph

from emmet.core.math import Vector3D
from emmet.core.qchem.calc_types import CalcType, LevelOfTheory, TaskType
from emmet.core.molecules.molecule_property import PropertyDoc
from emmet.core.mpid import MPID, MPculeID
from emmet.core.molecules.orbitals import (
    NaturalPopulation,
    LonePair,
    Bond,
    ThreeCenterBond,
    Hyperbond,
    Interaction,
)
from emmet.core.molecules.metal_binding import MetalBindingData


__author__ = "Evan Spotte-Smith <ewcspottesmith@lbl.gov>"


T = TypeVar("T", bound="MoleculeSummaryDoc")


class HasProps(Enum):
    """
    Enum of possible hasprops values.
    """

    molecules = "molecules"
    multipole_moments = "multipole_moments"
    partial_charges = "partial_charges"
    partial_spins = "partial_spins"
    redox = "redox"
    thermo = "thermo"
    vibration = "vibration"


class MoleculeSummaryDoc(PropertyDoc):
    """
    Summary information about molecules and their properties, useful for searching.
    """

    property_name: str = "summary"

    # molecules
    molecules: Dict[str, Molecule] = Field(
        ...,
        description="The lowest energy optimized structures for this molecule for each solvent.",
    )

    molecule_levels_of_theory: Optional[Dict[str, str]] = Field(
        None,
        description="Level of theory used to optimize the best molecular structure for each solvent.",
    )

    species_hash: Optional[str] = Field(
        None,
        description="Weisfeiler Lehman (WL) graph hash using the atom species as the graph "
        "node attribute.",
    )
    coord_hash: Optional[str] = Field(
        None,
        description="Weisfeiler Lehman (WL) graph hash using the atom coordinates as the graph "
        "node attribute.",
    )

    inchi: Optional[str] = Field(
        None, description="International Chemical Identifier (InChI) for this molecule"
    )
    inchi_key: Optional[str] = Field(
        None, description="Standardized hash of the InChI for this molecule"
    )

    task_ids: List[MPID] = Field(
        [],
        title="Calculation IDs",
        description="List of Calculation IDs associated with this molecule.",
    )

    unique_calc_types: Optional[List[CalcType]] = Field(
        None,
        description="Collection of all unique calculation types used for this molecule",
    )

    unique_task_types: Optional[List[TaskType]] = Field(
        None,
        description="Collection of all unique task types used for this molecule",
    )

    unique_levels_of_theory: Optional[List[LevelOfTheory]] = Field(
        None,
        description="Collection of all unique levels of theory used for this molecule",
    )

    unique_solvents: Optional[List[str]] = Field(
        None,
        description="Collection of all unique solvents (solvent parameters) used for this molecule",
    )

    unique_lot_solvents: Optional[List[str]] = Field(
        None,
        description="Collection of all unique combinations of level of theory and solvent used for this molecule",
    )

    # thermo
    thermo_property_ids: Optional[Dict[str, str]] = Field(
        None,
        description="Solvent:property ID map for each MoleculeThermoDoc for this molecule.",
    )

    thermo_levels_of_theory: Optional[Dict[str, str]] = Field(
        None,
        description="Solvent:level of theory map for each MoleculeThermoDoc for this molecule.",
    )

    electronic_energy: Optional[Dict[str, float]] = Field(
        None, description="Electronic energy of the molecule (units: eV)"
    )

    zero_point_energy: Optional[Dict[str, Optional[float]]] = Field(
        None, description="Zero-point energy of the molecule (units: eV)"
    )

    total_enthalpy: Optional[Dict[str, Optional[float]]] = Field(
        None, description="Total enthalpy of the molecule at 298.15K (units: eV)"
    )
    total_entropy: Optional[Dict[str, Optional[float]]] = Field(
        None, description="Total entropy of the molecule at 298.15K (units: eV/K)"
    )

    translational_enthalpy: Optional[Dict[str, Optional[float]]] = Field(
        None,
        description="Translational enthalpy of the molecule at 298.15K (units: eV)",
    )
    translational_entropy: Optional[Dict[str, Optional[float]]] = Field(
        None,
        description="Translational entropy of the molecule at 298.15K (units: eV/K)",
    )
    rotational_enthalpy: Optional[Dict[str, Optional[float]]] = Field(
        None, description="Rotational enthalpy of the molecule at 298.15K (units: eV)"
    )
    rotational_entropy: Optional[Dict[str, Optional[float]]] = Field(
        None, description="Rotational entropy of the molecule at 298.15K (units: eV/K)"
    )
    vibrational_enthalpy: Optional[Dict[str, Optional[float]]] = Field(
        None, description="Vibrational enthalpy of the molecule at 298.15K (units: eV)"
    )
    vibrational_entropy: Optional[Dict[str, Optional[float]]] = Field(
        None, description="Vibrational entropy of the molecule at 298.15K (units: eV/K)"
    )

    free_energy: Optional[Dict[str, Optional[float]]] = Field(
        None, description="Gibbs free energy of the molecule at 298.15K (units: eV)"
    )

    # vibrational properties
    vibration_property_ids: Optional[Dict[str, str]] = Field(
        None,
        description="Solvent:property ID map for each VibrationDoc for this molecule.",
    )

    vibration_levels_of_theory: Optional[Dict[str, str]] = Field(
        None,
        description="Solvent:level of theory map for each VibrationDoc for this molecule.",
    )

    frequencies: Optional[Dict[str, List[float]]] = Field(
        None, description="List of molecular vibrational frequencies"
    )

    frequency_modes: Optional[Dict[str, List[List[List[float]]]]] = Field(
        None,
        description="Vibrational frequency modes of the molecule (units: Angstrom)",
    )

    ir_intensities: Optional[Dict[str, List[float]]] = Field(
        None,
        title="IR intensities",
        description="Intensities for infrared vibrational spectrum peaks",
    )

    ir_activities: Optional[Dict[str, List]] = Field(
        None,
        title="IR activities",
        description="List indicating if frequency-modes are IR-active",
    )

    raman_intensities: Optional[Dict[str, List[float]]] = Field(
        None,
        title="Raman intensities",
        description="Intensities for Raman spectrum peaks",
    )

    raman_activities: Optional[Dict[str, List]] = Field(
        None,
        title="Raman activities",
        description="List indicating if frequency-modes are Raman-active",
    )

    # partial charges
    partial_charges_property_ids: Optional[Dict[str, Dict[str, str]]] = Field(
        None,
        description="Solvent:method:property ID map for each PartialChargesDoc for this molecule.",
    )

    partial_charges_levels_of_theory: Optional[Dict[str, Dict[str, str]]] = Field(
        None,
        description="Solvent:method:level of theory map for each PartialChargesDoc for this molecule.",
    )

    partial_charges: Optional[Dict[str, Dict[str, List[float]]]] = Field(
        None,
        description="Atomic partial charges for the molecule using different partitioning schemes "
        "(Mulliken, Restrained Electrostatic Potential, Natural Bonding Orbitals, etc.)",
    )

    # partial spins
    partial_spins_property_ids: Optional[Dict[str, Dict[str, str]]] = Field(
        None,
        description="Solvent:method:property ID map for each PartialSpinsDoc for this molecule.",
    )

    partial_spins_levels_of_theory: Optional[Dict[str, Dict[str, str]]] = Field(
        None,
        description="Solvent:method:level of theory map for each PartialSpinsDoc for this molecule.",
    )

    partial_spins: Optional[Dict[str, Dict[str, List[float]]]] = Field(
        None,
        description="Atomic partial spins for the molecule using different partitioning schemes "
        "(Mulliken, Natural Bonding Orbitals, etc.)",
    )

    # electric multipoles
    multipole_moments_property_ids: Optional[Dict[str, str]] = Field(
        None,
        description="Solvent:method:property ID map for each ElectricMultipoleDoc for this molecule.",
    )

    multipole_moments_levels_of_theory: Optional[Dict[str, str]] = Field(
        None,
        description="Solvent:method:level of theory map for each ElectricMultipoleDoc for this molecule.",
    )

    total_dipole: Optional[Dict[str, float]] = Field(
        None,
        description="Total molecular dipole moment (Debye)",
    )

    dipole_moment: Optional[Dict[str, Vector3D]] = Field(
        None,
        description="Molecular dipole moment vector (Debye)",
    )

    resp_total_dipole: Optional[Dict[str, float]] = Field(
        None,
        description="Total dipole moment, calculated via restrained electrostatic potential (RESP) (Debye)",
    )

    resp_dipole_moment: Optional[Dict[str, Vector3D]] = Field(
        None,
        description="Molecular dipole moment vector, calculated via RESP (Debye)",
    )

    quadrupole_moment: Optional[Dict[str, Dict[str, float]]] = Field(
        None,
        description="Quadrupole moment components (Debye Ang)",
    )

    octopole_moment: Optional[Dict[str, Dict[str, float]]] = Field(
        None,
        description="Octopole moment components (Debye Ang^2)",
    )

    hexadecapole_moment: Optional[Dict[str, Dict[str, float]]] = Field(
        None,
        description="Hexadecapole moment tensor components (Debye Ang^2)",
    )

    # redox properties
    redox_property_ids: Optional[Dict[str, str]] = Field(
        None, description="Solvent:property ID map for each RedoxDoc for this molecule."
    )

    redox_levels_of_theory: Optional[Dict[str, str]] = Field(
        None,
        description="Solvent:level of theory map for each RedoxDoc for this molecule.",
    )

    electron_affinity: Optional[Dict[str, float]] = Field(
        None, description="Vertical electron affinity in eV"
    )

    ea_task_id: Optional[Dict[str, MPID]] = Field(
        None, description="Molecule ID for electron affinity"
    )

    ionization_energy: Optional[Dict[str, float]] = Field(
        None, description="Vertical ionization energy in eV"
    )

    ie_task_id: Optional[Dict[str, MPID]] = Field(
        None, description="Molecule ID for ionization energy"
    )

    reduction_free_energy: Optional[Dict[str, float]] = Field(
        None, description="Adiabatic free energy of reduction"
    )

    red_molecule_id: Optional[Dict[str, MPculeID]] = Field(
        None, description="Molecule ID for adiabatic reduction"
    )

    oxidation_free_energy: Optional[Dict[str, float]] = Field(
        None, description="Adiabatic free energy of oxidation"
    )

    ox_molecule_id: Optional[Dict[str, MPculeID]] = Field(
        None, description="Molecule ID for adiabatic oxidation"
    )

    reduction_potential: Optional[Dict[str, float]] = Field(
        None,
        description="Reduction potential referenced to the standard hydrogen electrode (SHE) (units: V)",
    )

    oxidation_potential: Optional[Dict[str, float]] = Field(
        None,
        description="Oxidation potential referenced to the standard hydrogen electrode (SHE) (units: V)",
    )

    # has props
    has_props: Optional[Dict[str, bool]] = Field(
        None,
        description="Properties available for this molecule",
    )

    @classmethod
    def from_docs(cls, molecule_id: MPculeID, docs: Dict[str, Any]):
        """Converts a bunch of property docs into a SummaryDoc"""

        doc = _copy_from_doc(docs)

        if len(doc["has_props"]) == 0:
            raise ValueError("Missing minimal properties!")

        id_string = f"summary-{molecule_id}"
        h = blake2b()
        h.update(id_string.encode("utf-8"))
        property_id = h.hexdigest()
        doc["property_id"] = property_id

        return MoleculeSummaryDoc(molecule_id=molecule_id, **doc)


# Key mapping
summary_fields: Dict[str, list] = {
    HasProps.molecules.value: [
        "charge",
        "spin_multiplicity",
        "natoms",
        "elements",
        "nelements",
        "composition",
        "composition_reduced",
        "formula_alphabetical",
        "chemsys",
        "symmetry",
        "molecules",
        "deprecated",
        "task_ids",
        "species_hash",
        "coord_hash",
        "inchi",
        "inchi_key",
        "unique_calc_types",
        "unique_task_types",
        "unique_levels_of_theory",
        "unique_solvents",
        "unique_lot_solvents",
        "molecule_levels_of_theory",
    ],
    HasProps.thermo.value: [
        "electronic_energy",
        "zero_point_energy",
        "total_enthalpy",
        "total_entropy",
        "translational_enthalpy",
        "translational_entropy",
        "rotational_enthalpy",
        "rotational_entropy",
        "vibrational_enthalpy",
        "vibrational_entropy",
        "free_energy",
    ],
    HasProps.vibration.value: [
        "frequencies",
        "frequency_modes",
        "ir_intensities",
        "ir_activities",
        "raman_intensities",
        "raman_activities",
    ],
    HasProps.partial_charges.value: ["partial_charges"],
    HasProps.partial_spins.value: ["partial_spins"],
    HasProps.multipole_moments.value: [
        "total_dipole",
        "dipole_moment",
        "resp_total_dipole",
        "resp_dipole_moment",
        "quadrupole_moment",
        "octopole_moment",
        "hexadecapole_moment",
    ],
    HasProps.redox.value: [
        "electron_affinity",
        "ea_task_id",
        "ionization_energy",
        "ie_task_id",
        "reduction_free_energy",
        "red_molecule_id",
        "oxidation_free_energy",
        "ox_molecule_id",
        "reduction_potential",
        "oxidation_potential",
    ],
}


def _copy_from_doc(doc: Dict[str, Any]):
    """Helper function to copy the list of keys over from amalgamated document"""

    # Doc format:
    # {property0: {...},
    #  property1: {solvent1: {...}, solvent2: {...}},
    #  property2: {solvent1: [{...}, {...}], solvent2: [{...}, {...}]}
    # }

    has_props = {str(val.value): False for val in HasProps}
    d = {"has_props": has_props, "origins": []}

    # Function to grab the keys and put them in the root doc
    for doc_key in summary_fields:
        sub_doc = doc.get(doc_key, None)

        if doc_key == "molecules":
            # Molecules is special because there should only ever be one
            # MoleculeDoc for a given molecule
            # There are not multiple MoleculeDocs for different solvents
            if sub_doc is None:
                break

            d["has_props"][doc_key] = True
            for copy_key in summary_fields[doc_key]:
                d[copy_key] = sub_doc[copy_key]
        else:
            # No information for this particular set of properties
            # Shouldn't happen, but can
            if sub_doc is None:
                continue

            d["has_props"][doc_key] = True
            sd, by_method = sub_doc

            if isinstance(sd, dict) and len(sd) > 0:
                for copy_key in summary_fields[doc_key]:
                    d[copy_key] = dict()

                    if by_method:
                        for solvent, solv_entries in sd.items():
                            d[copy_key][solvent] = dict()
                            for method, entry in solv_entries.items():
                                if entry.get(copy_key) is not None:
                                    d[copy_key][solvent][method] = entry[copy_key]
                            if len(d[copy_key][solvent]) == 0:
                                # If this key was not populated at all for this solvent, get rid of it
                                del d[copy_key][solvent]
                    else:
                        for solvent, entry in sd.items():
                            if entry.get(copy_key) is not None:
                                d[copy_key][solvent] = entry[copy_key]

                    if len(d[copy_key]) == 0:
                        # If this key was not populated at all, set it to None
                        d[copy_key] = None

                # Populate property id and level of theory values
                d[doc_key + "_property_ids"] = dict()
                d[doc_key + "_levels_of_theory"] = dict()
                if by_method:
                    for solvent, solv_entries in sd.items():
                        d[doc_key + "_property_ids"][solvent] = dict()
                        d[doc_key + "_levels_of_theory"][solvent] = dict()
                        for method, entry in solv_entries.items():
                            d[doc_key + "_property_ids"][solvent][method] = entry[
                                "property_id"
                            ]
                            d[doc_key + "_levels_of_theory"][solvent][method] = entry[
                                "level_of_theory"
                            ]
                        if len(d[doc_key + "_property_ids"][solvent]) == 0:
                            del d[doc_key + "_property_ids"][solvent]
                        if len(d[doc_key + "_levels_of_theory"][solvent]) == 0:
                            del d[doc_key + "_levels_of_theory"][solvent]

                else:
                    for solvent, entry in sd.items():
                        d[doc_key + "_property_ids"][solvent] = entry["property_id"]
                        d[doc_key + "_levels_of_theory"][solvent] = entry[
                            "level_of_theory"
                        ]

                if len(d[doc_key + "_property_ids"]) == 0:
                    d[doc_key + "_property_ids"] = None
                if len(d[doc_key + "_levels_of_theory"]) == 0:
                    d[doc_key + "_levels_of_theory"] = None

    return d
