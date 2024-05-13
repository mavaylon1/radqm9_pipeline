from monty.serialization import loadfn
from tqdm import tqdm
import collections
from glob import glob
import numpy as np
import ase
import ase.io

from radqm9_pipeline.elements import read_elements

import os
CUR_DIR = os.path.dirname(os.path.realpath(__file__))

elements_dict = read_elements(os.path.join(CUR_DIR, 'elements.pkl'))

def merge_data(folder: str):
    """
    Load and merge the data into a single list from a folder of json files.
    """
    files = glob(folder+'/*')
    merged_data = []
    
    for file in tqdm(files):
        data = loadfn(file)
        merged_data+=data
    
    return merged_data

def bucket_mol_id(data: list):
    """
    Return: 
    - Good_data--> This means data that either had no duplicates or data that has been able to be mended 
    - bad_ids --> are events that are not able to be mended
    - length is the total number of events both mended and unmended
    
    TODO: fix the merge duplication in geometries, for now this does not matter since we downsample.
  
    """
    bucket = {}
    
    # Step 1    
    """
    Bucket into mol ids
    """
    for event in tqdm(data):
        try:
            bucket[event['mol_id']].append(event) 
        except KeyError:
            bucket[event['mol_id']] = [event]
    
    # Step 2
    """
    Find duplicate pairs in each mol id.
    What are duplicate pairs? These are training sessions that most likely continued as a separate job. 
    They need to be attached if possible. 
    """
    length=0
    bad_ids = []
    for mol_id in tqdm(bucket):
        pairs = [event['charge_spin'] for event in bucket[mol_id]]
        # Get a list of all charge_spins that have duplicates.
        duplicate_pairs = [item for item, count in collections.Counter(pairs).items() if count > 1] 
        if len(duplicate_pairs)!=0:
            """
            Handle the duplicate pairs to see if they can be merged
            """
            len_p = len(bucket[mol_id])
            for dup in duplicate_pairs:
                bad_data = []
                # Order events
                case_events = [event for event in bucket[mol_id] if event['charge_spin']==dup]
                bad_events = []
                bad_events += case_events

                for event in case_events: # remove events to be fixed. If fixable, add back at the end.
                    bucket[mol_id].remove(event) 
        
                ordered = [case_events[0]]
                del case_events[0]
                
                counter = 0
                threshold = 30
                while len(case_events)!=0:
                    if len(bad_data)==0:
                        for event in case_events:
                            beg = event['geometries'][0]
                            end = event['geometries'][len(event['geometries'])-1]

                            ordered_beg = ordered[0]['geometries'][0]
                            ordered_end = ordered[len(ordered)-1]['geometries'][len(ordered[len(ordered)-1]['geometries'])-1]

                            if beg==ordered_end:
                                ordered.append(event)
                                case_events.remove(event)
                            elif end==ordered_beg:
                                ordered.insert(0, event)
                                case_events.remove(event)
                            else:
                                counter+=1
                                if counter>threshold:
                                    bad_data.append(mol_id)
                                else:
                                    continue
                    else:
                        break
                    
                if len(bad_data)==0:                
                    # Merge the ordered events: forces, geometries
                    merged_event = {}
                    merged_event['task_id'] = ordered[0]['task_id']
                    merged_event['mol_id'] = mol_id
                    merged_event['name'] = ordered[0]['name']
                    merged_event['charge'] = ordered[0]['charge']
                    merged_event['spin'] = ordered[0]['spin']
                    merged_event['charge_spin'] = ordered[0]['charge_spin']
                    merged_event['species'] = ordered[0]['species']

                    geometries = []
                    energies = []
                    grads = []
                    mulliken = []
                    resp = []
                    dipole_moments = []
                    dipole_moments_resp = []
                    for event in ordered:
                        geometries += event['geometries']
                        energies += event['energies']
                        grads += event['gradients']
                        mulliken += event['mulliken']
                        resp += event['resp']
                        dipole_moments += event['dipole_moments']
                        dipole_moments_resp += event['dipole_moments_resp']

                    merged_event['geometries'] = geometries
                    merged_event['energies'] = energies
                    merged_event['gradients'] = grads
                    merged_event['mulliken'] = mulliken
                    merged_event['resp'] = resp
                    merged_event['dipole_moments'] = dipole_moments
                    merged_event['dipole_moments_resp'] = dipole_moments_resp

                    bucket[mol_id].append(merged_event)
                else:
                    bad_ids += bad_events
            len_r = len(bucket[mol_id])
            length += len_p-len_r
    good_data = flatten_filter(bucket)
        
    return good_data, bad_ids, length


def flatten_filter(data: dict):
    """
    Flatten bucket (dict) to list
    """
    data_to_be_parsed = []
    for mol_id in tqdm(data):
        for pair in data[mol_id]:
            data_to_be_parsed.append(pair)
    return data_to_be_parsed    

def add_unique_id(data: list):
    for item in data:
        item['mol_cs'] = str(item['mol_id']) + str(item['charge_spin'])

def average_force_trajectory(pair):
    """
    This method will take a specfic spin charge pair. At each point in the optimization trajectory, the 
    """
    forces = {}
    for i in range(len(pair['gradients'])):
        temp = []
        for atom in pair['gradients'][i]:
            res = np.sqrt(sum([j**2 for j in atom]))
            temp.append(res)
        forces[i] = np.mean(temp)
    del forces[0]
    return forces

def sparse_trajectory(bucket: list):
    """
    This takes the cleaned data and will sparsifiy the optimization trajectories. How this is done will depend on the
    charge_spin pair:
    - Neutral Singlet (0,1): First and Last
    - Other: First, Last, and structure with the highest molecular force other than the First.
    
    Note: Molecular Force is just the average of the force magnitudes of each atom in the molecule:
    """
    
    for pair in tqdm(bucket):
        if pair['charge_spin'] == '0,1':
            geometries = [pair['geometries'][0], pair['geometries'][-1]]
            energies = [pair['energies'][0], pair['energies'][-1]]
            grads = [pair['gradients'][0], pair['gradients'][-1]]
            mulliken = [pair['mulliken'][0], pair['mulliken'][-1]]
            resp = [pair['resp'][0], pair['resp'][-1]]
            dipole_moments = [pair['dipole_moments'][0], pair['dipole_moments'][-1]]
            dipole_moments_resp = [pair['dipole_moments_resp'][0], pair['dipole_moments_resp'][-1]]

            pair['geometries'] = geometries
            pair['energies'] = energies
            pair['gradients'] = grads
            pair['mulliken'] = mulliken
            pair['resp'] = resp
            pair['dipole_moments'] = dipole_moments
            pair['dipole_moments_resp'] = dipole_moments_resp
        else:
            force_dict = average_force_trajectory(pair)
            max_index = max(force_dict, key=force_dict.get)

            geometries = [pair['geometries'][0], pair['geometries'][max_index], pair['geometries'][-1]]
            energies = [pair['energies'][0], pair['energies'][max_index], pair['energies'][-1]]
            grads = [pair['gradients'][0], pair['gradients'][max_index], pair['gradients'][-1]]
            mulliken = [pair['mulliken'][0], pair['mulliken'][max_index], pair['mulliken'][-1]]
            resp = [pair['resp'][0], pair['resp'][max_index], pair['resp'][-1]]
            dipole_moments = [pair['dipole_moments'][0], pair['dipole_moments'][max_index], pair['dipole_moments'][-1]]
            dipole_moments_resp = [pair['dipole_moments_resp'][0], pair['dipole_moments_resp'][max_index], pair['dipole_moments_resp'][-1]]

            pair['geometries'] = geometries
            pair['energies'] = energies
            pair['gradients'] = grads
            pair['mulliken'] = mulliken
            pair['resp'] = resp
            pair['dipole_moments'] = dipole_moments
            pair['dipole_moments_resp'] = dipole_moments_resp
    
def force_magnitude_filter(cutoff: float,
                           data: list):
    """
    This method returns both data that meets the cuttoff value and data that is equal to or above the cuttoff value.
    If this is run before downsampling, it removes the entire data point trajectory.
    
    Returns: lists
    """
    good = []
    bad = []
    
    for item in tqdm(data):
        forces = item['gradients']
        for path_point in forces:
            next_item = False
            for atom in path_point:
                res = np.sqrt(sum([i**2 for i in atom]))
                if res >= cutoff:
                    bad.append(item)
                    next_item = True
                    break
            if next_item:
                break
        if not next_item:
            good.append(item)
                            
    return good, bad

def build_graph(species, position):
    atoms = ase.atoms.Atoms(symbols=species,
                            positions=position)
    mol = AseAtomsAdaptor.get_molecule(atoms)
    graph = MoleculeGraph.with_local_env_strategy(mol, OpenBabelNN())
    return graph

def filter_broken_graphs(data: list):
    broken = []
    good = []
    
    for item in tqdm(data):
        for traj_point in item['geometries'][1:]: # Ignore the first is in the filter_broken_graphs
            graph = build_graph(item['species'], traj_point)
            connected = nx.is_connected(graph.graph.to_undirected())
            if not connected:
                broken.append(item)
                break
        if connected:
            good.append(item)

    return good, broken



def charge_filter(charges: list, data):
    """
    Takes both a list of charges to filter by and a dataset that is a list of data points.
    """
    filtered_data = []
    for point in tqdm(data):
        if point['charge'] in charges:
            filtered_data.append(point)
    
    return filtered_data

def chunk_train_multiple(data: list, percentage: list):
    """
    Percentage: list of percentages e.g., [.05, .1, .25, .5, .75]
    """
    mol_id_bucket = {}
    for point in tqdm(data):
        try:
            mol_id_bucket[point['mol_id']].append(point)
        except KeyError:
            mol_id_bucket[point['mol_id']] = [point]
    
    
    elements_dict = read_elements('/pscratch/sd/m/mavaylon/sam_ldrd/radqm9_pipeline/src/radqm9_pipeline/modules/elements.pkl')
    
    ##################
    # Create a weight dictionary such that the keys are the unique weights using atomic mass of the molecule
    # as a float and the values are list of mol_ids that correspond to said weight. 
    ##################
    weight_dict = {}
    for point in tqdm(data):
        species = point['species']
        species_num = 0
        for element in species:
            species_num+=elements_dict[element]
        try:
            weight_dict[str(species_num)].append(point['mol_id'])
            weight_dict[str(species_num)] = list(set(weight_dict[str(species_num)]))
        except KeyError:
            weight_dict[str(species_num)] = [point['mol_id']]
    
    ##################
    # Calculate total data points
    ##################
    total=0
    for pair in tqdm(data):
        total+=len(pair['geometries'])
    print(total)
    
    ##################
    # Calculate the size for the data chunk for each percentage
    ##################
    sizes = []
    for item in percentage:
        temp_size = round(total*item)
        sizes.append(temp_size)
    
    ##################
    # Get the chunked mol_ids for each size
    ##################
    chunks = []
    count = 0
    
    for size in tqdm(sizes):
        print("size:", size)
        chunked_mol_id_data = []
        # while count<size:
        for i in tqdm(range(total)):
            if count<size:
                for key in weight_dict:
                    if len(weight_dict[key])!=0:
                        _id = weight_dict[key][0]
                        for point in mol_id_bucket[_id]:
                            count += len(point['geometries'])
                        # print(count)

                        chunked_mol_id_data.append(_id)
                        weight_dict[key] = weight_dict[key][1:]

                    else:
                        pass
                weight_dict = {k: v for (k,v) in weight_dict.items() if len(v)!=0}

            else:
                break
        print('count:', count)
        # print(len(chunked_mol_id_data))
        chunks.append(chunked_mol_id_data)
        
        for item in chunks:
            print(len(item))
        
    chunked_data = []
    for chunk_set in tqdm(chunks):
        chunk = []
        for item in tqdm(chunk_set):
            chunk+=mol_id_bucket[item]
        chunked_data.append(chunk)
    
    return chunked_data, weight_dict

def get_molecule_weight(data: list):
    """
    The method takes in a list of data (either trajectories or single points) and sorts into a distribution
    dictionary. The keys are the species/formula and the value of each key is the weight. appearing number of times the species
    appears in the dataset.
    """
    dict_dist = {}
    for item in tqdm(data):
        species_num = []
        species=''.join((sorted(item['species'])))
        
        for element in item['species']:
            species_num.append(elements_dict[element])

        species_sum = sum(species_num)
        try:
            dict_dist[species].append(species_sum)
            dict_dist[species] = [dict_dist[species][0]]*len(dict_dist[species])
        except KeyError:
            dict_dist[species] = [species_sum]
        
    return dict_dist

def molecule_weight(data: list, weight_dict):
    """
    This method takes in data and assigns the mass.
    Python does a weird thing floats e.g., {126.15499999999993, 126.15499999999994}, having this and
    get_molecule_weight gurantees that species that are the same are not being assigned different weights.
    """
    for item in tqdm(data):
        weight = weight_dict[''.join((sorted(item['species'])))][0]
        item['molecule_mass'] = weight
        
def weight_to_data(data: list):
    """
    This method buckets the data by the mass such that the dict key is the mass and the values are the data
    points.
    """
    dict_data = {}
    for item in tqdm(data):
        try:
            dict_data[item['molecule_mass']].append(item)
        except KeyError:
            dict_data[item['molecule_mass']] = [item]
    return dict_data

def length_dict(data: dict):
    """
    This method takes in the output of weight_to_data and returns a dictionary that is sorted from largest
    to smallest mass. The keys are the mass and the values are the number of appearances.
    """
    length_dict = {key: len(value) for key, value in data.items()}
    sorted_length_dict = {k: length_dict[k] for k in sorted(length_dict, reverse=True)}
    
    return sorted_length_dict

def build_minimal_atoms(data: dict,
                energy: str = None,
                forces: str = None,
                charge:str = None,
                spin:str = None,
                train = False) -> ase.Atoms:
    """ 
    Populate Atoms class with atoms in molecule.
        atoms.info : global variables
        atoms.array : variables for individual atoms
        
    Both "energy" and "forces" are the dict strings in data.
    """
    atom_list = []
    for i in range(len(data['geometries'])):
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
            if data['charge_spin'] == '0,1':
                atoms.info['position_type'] = 'end'
            else:
                atoms.info['position_type'] = 'middle'
        if i == 2:
            atoms.info['position_type'] = 'end'
        atom_list.append(atoms)
    return atom_list


def build_minimal_atoms_iterator(data: list,
                         train=False):
    """
    This method assumes the data has been validated. This will create ASE atoms to be written.
    
    The input needs to be a list of lists that contain the event dictionaries. Each inner list needs to represent all the events for a single
    mol_id.
    """
    data_set=[]
    for point in tqdm(data):
        atoms=build_minimal_atoms(point, energy='energies', forces='gradients', charge='charge', spin='spin', train=train)
        data_set+=atoms
    return data_set


def create_dataset(data: dict,
                   file_name:str,
                   path:str):
    """
    This method will handle the I/O for writing the data to xyz files to the path provided.
    """
    train_data = data['train']
    val_data = data['val']
    test_data = data['test']
    
    train_file = os.path.join(path,file_name+'_train.xyz')
    ase.io.write(train_file, train_data,format="extxyz")
     
    val_file = os.path.join(path,file_name+'_val.xyz')
    ase.io.write(val_file, val_data,format="extxyz")
    
    test_file = os.path.join(path,file_name+'_test.xyz')
    ase.io.write(test_file, test_data,format="extxyz")