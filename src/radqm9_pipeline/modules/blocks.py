from monty.serialization import loadfn
from tqdm import tqdm
import collections
from PyAstronomy import pyasl
from glob import glob
import numpy as np

from radqm9_pipeline.ase_data.data_to_atoms import build_atoms
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
    Method requires a list of events such that each item in the list is a single event.
    
    Steps:
    - Bucket the events in a dictionary such that the keys are the mol_ids and the values are all events with that mol_id. Create
      a new field for each event called "charge_spin" for easy down_stream checks.
    - For each mol_id, check the list of events for duplicate "charge_spin". 
      If there are duplicates, for each duplicate pair:
                1. Get all events that have the duplicate charge_spin
                2. Remove these events from the values of the mol_id in the bucket. This is to process them and re-add them later or move to
                   bad data.
                3. Search for the order in which these duplicate events for the specific charge_spin should be attached. If there is no match,
                   this is considered bad data and will be moved out of the bucket.
                4. If an order can be established, create a new "merged" event. 
                5. Append the event to the values of the mol_id of the bucket.
                6. Repeat for every duplicate pair.
    """
    bucket = {}
    
    # Step 1    
    for event in tqdm(data):
        event['charge_spin'] = str(event['charge'])+','+str(event['spin'])
        try:
            bucket[event['mol_id']].append(event) 
        except KeyError:
            bucket[event['mol_id']] = [event]

    # Step 2
    bad_ids = []
    for mol_id in tqdm(bucket):
        pairs = [event['charge_spin'] for event in bucket[mol_id]]
        # Get a list of all charge_spins that have duplicates.
        duplicate_pairs = [item for item, count in collections.Counter(pairs).items() if count > 1] 
        non_dup_pairs = []
        if len(duplicate_pairs)!=0:
            non_dup_pairs = list(set(pairs))
            for i in duplicate_pairs:
                non_dup_pairs.remove(i)
        
        if len(duplicate_pairs)==0:
            # replace each charge_spin event with the beginning, middle, end
            for pair in bucket[mol_id]:
                middle_index = int(len(pair['geometries'])/2)
                reference_index = pair['energies'].index(sorted(pair['energies'])[middle_index])
                
                geometries = [pair['geometries'][0], pair['geometries'][reference_index], pair['geometries'][-1]]
                energies = [pair['energies'][0], pair['energies'][reference_index], pair['energies'][-1]]
                grads = [pair['gradients'][0], pair['gradients'][reference_index], pair['gradients'][-1]]
                mulliken = [pair['mulliken'][0], pair['mulliken'][reference_index], pair['mulliken'][-1]]
                resp = [pair['resp'][0], pair['resp'][reference_index], pair['resp'][-1]]
                
                # add weight tag
                species = pair['species']
                species_num = []
                species_sorted = ''.join(sorted(set(species)))
                for element in species:
                    species_num.append(elements_dict[element])
                    
                pair['weight_tag'] = round(sum(species_num))
                
                pair['geometries'] = geometries
                pair['energies'] = energies
                pair['gradients'] = grads
                pair['mulliken'] = mulliken
                pair['resp'] = resp
                
        else:
            if len(non_dup_pairs)!=0:
                for pair in non_dup_pairs:
                    case_events = [event for event in bucket[mol_id] if event['charge_spin']==pair]
                    for event in case_events:
                        middle_index = int(len(event['geometries'])/2)
                        reference_index = event['energies'].index(sorted(event['energies'])[middle_index])

                        geometries = [event['geometries'][0], event['geometries'][reference_index], event['geometries'][-1]]
                        energies = [event['energies'][0], event['energies'][reference_index], event['energies'][-1]]
                        grads = [event['gradients'][0], event['gradients'][reference_index], event['gradients'][-1]]
                        mulliken = [event['mulliken'][0], event['mulliken'][reference_index], event['mulliken'][-1]]
                        resp = [event['resp'][0], event['resp'][reference_index], event['resp'][-1]]

                        event['geometries'] = geometries
                        event['energies'] = energies
                        event['gradients'] = grads
                        event['mulliken'] = mulliken
                        event['resp'] = resp
                        
                        # add weight tag
                        species = event['species']
                        species_num = []
                        species_sorted = ''.join(sorted(set(species)))
                        for element in species:
                            species_num.append(elements_dict[element])

                        event['weight_tag'] = round(sum(species_num))
                

            bad_data = []
            for dup in duplicate_pairs:
                # Order events
                case_events = [event for event in bucket[mol_id] if event['charge_spin']==dup]
                for event in case_events:
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
                    pass
                else:
                    bad_ids += bad_data
                    continue

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
                for event in ordered:
                    geometries += event['geometries']
                    energies += event['energies']
                    grads += event['gradients']
                    mulliken += event['mulliken']
                    resp += event['resp']
                
                #sort the beginning, middle, end
                middle_index = int(len(geometries)/2)
                reference_index = energies.index(sorted(energies)[middle_index])
                
                geometries = [geometries[0], geometries[reference_index], geometries[-1]]
                energies = [energies[0], energies[reference_index], energies[-1]]
                grads = [grads[0], grads[reference_index], grads[-1]]
                mulliken = [mulliken[0], mulliken[reference_index], mulliken[-1]]
                resp = [resp[0], resp[reference_index], resp[-1]]
               
                merged_event['geometries'] = geometries
                merged_event['energies'] = energies
                merged_event['gradients'] = grads
                merged_event['mulliken'] = mulliken
                merged_event['resp'] = resp
                
                species = merged_event['species']
                species_num = []
                species_sorted = ''.join(sorted(set(species)))
                for element in species:
                    species_num.append(elements_dict[element])
                merged_event['weight_tag'] = round(sum(species_num))
                
                bucket[mol_id].append(merged_event)

    if len(bad_ids)!=0:
        ids = list(set(bad_ids))
        for _id in ids:
            bucket.pop(_id)
        
        return ids, bucket
    else:
        for mol_id in bucket:
            if len(bucket[mol_id])>4:
                print(mol_id)
        return bucket

    
def force_magnitude_filter(cutoff: float,
                           data: dict):
    """
    This method returns both dict of mol_id, force_values, and general info of data that is equal to or above the cuttoff value.
    Cutoff: <float> ev per angstrom 
    
    Returns: Dict
    """
    force_dict = {}
    
    for mol_id in tqdm(data):
        for config in data[mol_id]:
            forces = config['gradients']
            for path_point in forces:
                for atom in path_point:
                    res = np.sqrt(sum([i**2 for i in atom]))
                    if res >= cutoff:
                        try:
                            force_dict['removed_forces'].append(res)
                            force_dict['removed_mol_ids'].append(mol_id)
                            force_dict['info'].append([mol_id, config['charge_spin'], res])
                        except KeyError:
                            force_dict['removed_forces'] = [res]
                            force_dict['removed_mol_ids'] = [mol_id]
                            force_dict['info'] = [mol_id, config['charge_spin'], res]
    return force_dict
    

def __mol_id_weight_bins(bucket: dict):
    """
    This method takes in the output from bucket_mol_id.
    For each mol_id, calculate the molecule weight based off the atoms. Combine that with that type of atoms used.
    This weight+type serves as a key for a dict. The values are then the mol_ids that match the key.
    
    The intent is create a dict such that we can sample from evenly based on weight. We also want to ensure even
    representation of atom type across train/val/test, hence why we include the atoms used in the key.
    """
    
    weight_dist = {}
    weight_dict = {}
    for mol_id in tqdm(bucket):
        species = bucket[mol_id][0]['species']
        species_num = []
        species_sorted = ''.join(sorted(set(species)))
        for element in species:
            species_num.append(elements_dict[element])

        try:
            weight_dist[round(sum(species_num))]+=1
        except KeyError:
            weight_dist[round(sum(species_num))]=1
            
        try:
            weight_dict[str(sum(species_num))+'_'+species_sorted].append(mol_id)

        except KeyError:
            weight_dict[str(sum(species_num))+'_'+species_sorted] = [mol_id]
    
    return weight_dict, weight_dist

def train_val_test_split(bucket: dict,
                         train_size: float,
                         val_size: float):
    """
    This method takes in the output from mol_id_weight_bins.
    This method will sample from each key-value pair from the input dict based on the train_size, val_size.
    The method requires a validation set, but the user can combine it with test if they so choose.
    """
    
    weight_dict, weight_dist = __mol_id_weight_bins(bucket)
    
    train_marker = train_size
    val_marker = train_size + val_size
    
    split={}

    import random
    random.seed(10)
    for strata in tqdm(weight_dict):
        random.shuffle(weight_dict[strata])
        train_index = round(len(weight_dict[strata])*train_marker)
        val_index = round(len(weight_dict[strata])*val_marker)
        
        try:
            split['train']+=(weight_dict[strata][:train_index])
            split['val']+=(weight_dict[strata][train_index:val_index])
            split['test']+=(weight_dict[strata][val_index:])
        except KeyError:
            split['train'] = weight_dict[strata][:train_index]
            split['val'] = weight_dict[strata][train_index:val_index]
            split['test'] = weight_dict[strata][val_index:]
    
    train_data = [bucket[i] for i in split['train']]
    val_data = [bucket[i] for i in split['val']]
    test_data = [bucket[i] for i in split['test']]
                 
    split['train']=train_data
    split['val']=val_data
    split['test']=test_data
    
    return split, weight_dist

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

def flatten(data):
    return [x for y in data for x in y]

def tagger(train_data, dist):
    
    # Flatten train_data into a pool of data points
    flat = [x for y in train_data for x in y] 
    
    # Find all the weights that appear less than 20 occurances
    cutoff_bin = [x for x in dist if .05*dist[x]<1]
    
    bins = {}
    cutoff_data = {}
    
    for point in flat:
        try:
            bins[point['weight_tag']].append(point['mol_id'])
        except KeyError:
            bins[point['weight_tag']] = [point['mol_id']]
        
    
    mol_id_tag = {}
    
    for weight_group in tqdm(bins):
        tag_groups = list(split(list(set(bins[weight_group])),20))
        chunk = 0
        for group in tag_groups:
            for mol_id in group:
                mol_id_tag[mol_id] = chunk
            chunk += 1

    for point in tqdm(flat):
        if point['mol_id'] in mol_id_tag:
            point['chunk'] = mol_id_tag[point['mol_id']]      
    
    return flat

def build_atoms(data: dict,
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
        if train:
            atoms.info['chunk'] = data['chunk']
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

def build_manager(data: dict, weight_dist, train):
    """
    Manage building atoms for train/val/test splits
    """
    
    data['train'] = tagger(data['train'], weight_dist)
    data['val'] = flatten(data['val'])
    data['test'] = flatten(data['test'])
    
    build = {}
    for split in data:
        if split == 'train':
            build[split] = build_atoms_iterator(data[split], train=train)
        else:
            build[split] = build_atoms_iterator(data[split])
    return build

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