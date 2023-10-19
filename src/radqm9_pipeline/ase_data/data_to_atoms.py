import ase
import copy
import numpy as np

def build_atoms(data: dict,
                energy: str = None,
                forces: str = None,
                charge:str = None,
                spin:str = None) -> ase.Atoms:
    """ 
    Populate Atoms class with atoms in molecule.
        atoms.info : global variables
        atoms.array : variables for individual atoms
        
    Both "energy" and "forces" are the dict strings in data.
    """
    atom_list = []
    for i in range(3):
        atoms = ase.atoms.Atoms(
            symbols=data['species'],
            positions=data['geometries'][i]
        )
        if energy is not None:
            atoms.info['energy'] = data[energy][i]
        if forces is not None:
            atoms.arrays['forces'] = np.array(data[forces][i])
        if charge is not None:
             atoms.info['charge'] = data[charge]
        if spin is not None:
            atoms.info['spin'] = data[spin]
        if i == 0:
            atoms.info['position_type'] = 'start'
        if i == 1:
            atoms.info['position_type'] = 'energetic_middle'
        if i == 2:
            atoms.info['position_type'] = 'minimum'
        atom_list.append(atoms)
    return atom_list

def build_atoms_iterator(data: list):
    """
    This method assumes the data has been validated. This will create ASE atoms to be written.
    
    The input needs to be a list of lists that contain the event dictionaries. Each inner list needs to represent all the events for a single
    mol_id.
    """
    data_set=[]
    for mol_id_list in tqdm(data):
        for pair in mol_id_list:
            atoms=build_atoms(pair, energy='energies', forces='gradients', charge='charge', spin='spin')
            data_set+=atoms
    return data_set