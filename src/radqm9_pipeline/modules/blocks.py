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

def bucket_mol_id_all_pts(data: list):
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
            for pair in bucket[mol_id]:
                # add weight tag
                species = pair['species']
                species_num = []
                species_sorted = ''.join(sorted(set(species)))
                for element in species:
                    species_num.append(elements_dict[element])
                    
                pair['weight_tag'] = round(sum(species_num))
        else:
            if len(non_dup_pairs)!=0:
                for pair in non_dup_pairs:
                    case_events = [event for event in bucket[mol_id] if event['charge_spin']==pair]
                    for event in case_events:
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

def sparse_trajectory(bucket: dict):
    """
    This takes the cleaned data and will sparsifiy the optimization trajectories. How this is done will depend on the
    charge_spin pair:
    - Neutral Singlet (0,1): First and Last
    - Other: First, Last, and structure with the highest molecular force other than the First.
    
    Note: Molecular Force is just the average of the force magnitudes of each atom in the molecule:
    """
    
    for mol_id in tqdm(bucket):
        for pair in bucket[mol_id]:
            if pair['charge_spin'] == '0,1':
                geometries = [pair['geometries'][0], pair['geometries'][-1]]
                energies = [pair['energies'][0], pair['energies'][-1]]
                grads = [pair['gradients'][0], pair['gradients'][-1]]
                mulliken = [pair['mulliken'][0], pair['mulliken'][-1]]
                resp = [pair['resp'][0], pair['resp'][-1]]
                
                pair['geometries'] = geometries
                pair['energies'] = energies
                pair['gradients'] = grads
                pair['mulliken'] = mulliken
                pair['resp'] = resp
            else:
                force_dict = average_force_trajectory(pair)
                max_index = max(force_dict, key=force_dict.get)
                
                geometries = [pair['geometries'][0], pair['geometries'][max_index], pair['geometries'][-1]]
                energies = [pair['energies'][0], pair['energies'][max_index], pair['energies'][-1]]
                grads = [pair['gradients'][0], pair['gradients'][max_index], pair['gradients'][-1]]
                mulliken = [pair['mulliken'][0], pair['mulliken'][max_index], pair['mulliken'][-1]]
                resp = [pair['resp'][0], pair['resp'][max_index], pair['resp'][-1]]
                
                pair['geometries'] = geometries
                pair['energies'] = energies
                pair['gradients'] = grads
                pair['mulliken'] = mulliken
                pair['resp'] = resp
    
def force_magnitude_filter(cutoff: float,
                           data: dict):
    """
    This method returns both data that meets the cuttoff value and data that is equal to or above the cuttoff value.
    
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
                        # print(mol_id, config)
                        try:
                            data[mol_id].remove(config)
                        except ValueError:
                            pass
                    else:
                        try:
                            force_dict['filtered_force_magnitude'].append(res)
                        except KeyError:
                            force_dict['filtered_force_magnitude'] = [res]
                            
    return force_dict, data

def prepare_graph_filter(data: dict):
    """
    Go through the data and retrieve the data needed for graph separation check.
    1. Ignore charge_spin pair = 0,1
    2. Ignore all starting points in the configurations.
    """
    data_to_be_parsed = []
    for mol_id in tqdm(data):
        for pair in data[mol_id]:
            if pair['charge_spin'] != '0,1':
                data_to_be_parsed.append(pair)
    
    # Ignore the first is in the filter_broken_graphs
    return data_to_be_parsed

def filter_broken_graphs(data: list):
    broken = []
    good = []
    
    for item in tqdm(data):
        for traj_point in item['geometries'][1:]:
                atoms = ase.atoms.Atoms(symbols=item['species'],
                                        positions=traj_point)
                mol = AseAtomsAdaptor.get_molecule(atoms)
                graph = MoleculeGraph.with_local_env_strategy(mol, OpenBabelNN())

                connected = nx.is_connected(graph.graph.to_undirected())
                if not connected:
                    broken.append(item)
                else:
                    good.append(item)

    return good, broken

def removed_broken_graph_data(data: dict,
                              broken: list):
    """
    
    """
    data_to_be_parsed = {}
    for mol_id in tqdm(data):
        for pair in data[mol_id]:
            # data_to_be_parsed.append(pair)
            data_to_be_parsed[str(pair['task_id']) +'_'+str(pair['mol_id'])+'_'+pair['charge_spin']] = pair
    print(len(data_to_be_parsed))
    
    print(len(broken))
    for item in tqdm(broken):
        data_to_be_parsed.pop(item)
    
    print(len(data_to_be_parsed))
    
    return data_to_be_parsed

def __mol_id_weight_bins(data: dict):
    """
    This method takes in the output from removing the broken graphs.
    
    1. Bin the data by mol_ids in a dict.

    For each mol_id, calculate the molecule weight based off the atoms. Combine that with that type of atoms used.
    This weight+type serves as a key for a dict. The values are then the mol_ids that match the key.
    
    The intent is create a dict such that we can sample from evenly based on weight. We also want to ensure even
    representation of atom type across train/val/test, hence why we include the atoms used in the key.
    """
    bucket={}
    for item in tqdm(data):
        try:
            bucket[data[item]['mol_id']].append(data[item])
        except KeyError:
            bucket[data[item]['mol_id']] = [data[item]]
    
    
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
    
    return weight_dict, weight_dist, bucket


def train_val_test_split(bucket: dict,
                         train_size: float,
                         val_size: float):
    """
    This method takes in the output from mol_id_weight_bins.
    This method will sample from each key-value pair from the input dict based on the train_size, val_size.
    The method requires a validation set, but the user can combine it with test if they so choose.
    """
    
    weight_dict, weight_dist, bucket = __mol_id_weight_bins(bucket)
    
    train_marker = train_size
    val_marker = train_size + val_size
    
    split={}

    import random
    random.seed(10)
    for strata in tqdm(weight_dict):
        random.shuffle(weight_dict[strata])
        train_index = round(len(weight_dict[strata])*train_marker)
        val_index = round(len(weight_dict[strata])*val_marker)
        # print(len(weight_dict[strata]))
        # print(train_index)
        # print(val_index)
        # break
              
        
        try:
            train_split = (weight_dict[strata][:train_index])
            val_split = (weight_dict[strata][train_index:val_index+1])
            test_split = (weight_dict[strata][val_index+1:])
            
            if len(test_split)> len(val_split):
                print('bleh')
                return [weight_dict[strata], train_split, val_split, test_split, train_index, val_index+1]
            
            split['train']+=train_split
            split['val']+=val_split
            split['test']+=test_split
            
            
        except KeyError:
            split['train'] = weight_dict[strata][:train_index]
            split['val'] = weight_dict[strata][train_index:val_index+1]
            split['test'] = weight_dict[strata][val_index+1:]
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
    # cutoff_bin = [x for x in dist if .05*dist[x]<1]
    
    bins = {}
    chunked_training = {}
    
    for point in flat:
        try:
            bins[point['weight_tag']].append(point['mol_id'])
        except KeyError:
            bins[point['weight_tag']] = [point['mol_id']]
        
    
    mol_id_tag = {}
    chunks = []
    for weight_group in tqdm(bins):
        tag_groups = list(split(list(set(bins[weight_group])),20))
        chunk = 0
        for group in tag_groups:
            for mol_id in group:
                chunks.append(chunk)
                mol_id_tag[mol_id] = chunk
            chunk += 5
    print(set(chunks))
    
    for point in tqdm(flat):
        if point['mol_id'] in mol_id_tag:
            point['chunk'] = mol_id_tag[point['mol_id']] 
            try:
                chunked_training[mol_id_tag[point['mol_id']]].append(point)
            except KeyError:
                chunked_training[mol_id_tag[point['mol_id']]] = [point]
    return flat, chunked_training

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
    try:
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
            if train:
                atoms.info['chunk'] = data['chunk']
            if i == 0:
                atoms.info['position_type'] = 'start'
            atoms.info['mol_id'] = data['mol_id']
            if i == 1:
                if data['charge_spin'] == '0,1':
                    atoms.info['position_type'] = 'end'
                else:
                    atoms.info['position_type'] = 'middle'
            if i == 2:
                atoms.info['position_type'] = 'end'
            atom_list.append(atoms)
    except IndexError:
        print(i)
        print(data['mol_id'])

def build_atoms_iterator(data: list,
                         train=False):
    """
    This method assumes the data has been validated. This will create ASE atoms to be written.
    
    The input needs to be a list of lists that contain the event dictionaries. Each inner list needs to represent all the events for a single
    mol_id.
    """
    data_set=[]
    for point in tqdm(data):
        atoms=build_atoms(point, energy='energies', forces='gradients', charge='charge', spin='spin', train=train)
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