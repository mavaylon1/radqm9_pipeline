import collections
from radqm9_pipeline.elements import read_elements

def charge_spin_validator(bucket: dict):
    """
    The method takes in data that has been sorted into a python dictionary. They keys are the mol_ids and the values are the list of events.
    This method will check for duplicate charge_spin pairs.
    
    return: None (if valid) | bad_data-->list of invalid mol_ids (if invalid)
    """

    for mol_id in bucket:
        pairs = [event['charge_spin'] for event in bucket[mol_id]]
        
        # Get a list of all charge_spins that have duplicates.
        duplicate_pairs = [item for item, count in collections.Counter(pairs).items() if count > 1]
        
        if len(duplicate_pairs)==0:
            pass
        else:
            bad_data.append(mol_id)
    
    if len(bad_data)==0:
        return None
    else:
        return bad_data

def composition_validator(train, val, test):
    """
    Each are lists. Validate atomic weight distribution and element distribution.
    
    - check that for each mol_id there are four instances enclosed in a split
    
    - Weight count
    - Number of element occurances in each split
    """
    elements_dict = read_elements('./src/radqm9_pipeline/elements.pkl')
    # weight count
    sets = [train, val, test]
    an = pyasl.AtomicNo()
    
    total = [item for sublist in sets for item in sublist]
    total_dict = {}
    total_element_dict = {}
    for mol_id in tqdm(total):
        # check for four charge_spin pairs
        if len(mol_id)!=4:
            #add
            pass

        # weight distribution
        species = mol_id[0]['species']
        species_num = []
        species_sorted = ''.join(sorted(set(species)))
        for element in species:
            species_num.append(elements_dict[element])
        try:
            total_dict[str(sum(species_num))]+= 1
        except KeyError:
            total_dict[str(sum(species_num))]=1
        
        for element in list(set(species)):
            try:
                total_element_dict[element] += 1
            except KeyError:
                total_element_dict[element] =1 
    # return total_dict
    
    set_weight = []
    element_list = []
    for split in sets:
        weight_dict = {}
        ratio_dict = {}
        split_element_dict = {}
        ratio_element_dict = {}
        for mol_id in tqdm(split):
            # check for four charge_spin pairs
            if len(mol_id)!=4:
                #add
                pass
            
            # weight distribution
            species = mol_id[0]['species']
            species_num = []
            species_sorted = ''.join(sorted(set(species)))
            for element in species:
                species_num.append(elements_dict[element])
            try:
                weight_dict[str(sum(species_num))]+= 1
            except KeyError:
                weight_dict[str(sum(species_num))]=1
            
            for element in list(set(species)):
                try:
                    split_element_dict[element] += 1
                except KeyError:
                    split_element_dict[element] =1 
        
        for key in weight_dict:
            try:
                ratio = weight_dict[key]/total_dict[key]
                ratio_dict[key] = ratio
            except ZeroDivisionError:
                return key
            
        set_weight.append(ratio_dict)
        
        for element in split_element_dict:
            element_ratio = split_element_dict[element]/total_element_dict[element]
            ratio_element_dict[element] = element_ratio
        
        element_list.append(ratio_element_dict)
    
    return set_weight, total_dict , total_element_dict, element_list

def format_validator(data: list):
    """
    This is the final check to make sure the format of each potential data point is correct. 
    This will check to see there is only a beginning, middle, and end for the positions, energies, and forces.
    
    The input needs to be a list of lists that contain the event dictionaries. Each inner list needs to represent all the events for a single
    mol_id.
    """
    for mol_id_list in tqdm(data):
        for event in mol_id_list:
            if len(event['geometries'])!=3: #TODO: make 'geometries', 'energies', 'gradients' function parameters to generalize
                msg = 'The positions are not length 3'
                raise ValueError(msg)
            if len(event['energies'])!=3:
                msg = 'The energies are not length 3'
                raise ValueError(msg)
            if len(event['gradients'])!=3:
                msg = 'The forces are not length 3'
                raise ValueError(msg)