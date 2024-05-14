import logging
import ast
import numpy as np
import json
import random
import tqdm
from glob import glob
import h5py
from ase.io import read
import torch
import multiprocessing as mp
import os
from typing import List, Tuple


from mace.tools import to_numpy
from mace import tools, data
from mace.data.utils import (
    save_AtomicData_to_HDF5,
    save_configurations_as_HDF5,
)
from mace.tools.scripts_utils import get_dataset_from_xyz, get_atomic_energies
from mace.tools.utils import AtomicNumberTable
from mace.tools import torch_geometric
from mace.modules import compute_statistics    

from data_arg_parser import build_default_arg_parser


args = build_default_arg_parser().parse_args()

try:
    config_type_weights = ast.literal_eval(args.config_type_weights)
    assert isinstance(config_type_weights, dict)
except Exception as e:  # pylint: disable=W0703
    # Very 
    config_type_weights = {"Default": 1.0}

atomic_energies_dict, all_data_configs = data.load_from_xyz(
        file_path=args.file,
        config_type_weights=config_type_weights,
        energy_key=args.energy_key,
        forces_key=args.forces_key,
        stress_key=args.stress_key,
        virials_key=args.virials_key,
        dipole_key=args.dipole_key,
        charges_key=args.charges_key,
        extract_atomic_energies=True,
    )

#########################
# H5 dataset multiprocess
#########################

# # split collections.train into batches and save them to hdf5
# split_train = np.array_split(all_data_configs, args.num_process)
# drop_last = False
# if len(all_data_configs) % 2 == 1:
#     drop_last = True

# # Define Task for Multiprocessiing
# def multi_train_hdf5(process):
#     with h5py.File(args.directory + '/' + args.prefix + str(process)+".h5", "w") as f:
#         f.attrs["drop_last"] = drop_last
#         save_configurations_as_HDF5(split_train[process], process, f)

# processes = []
# for i in range(args.num_process):
#     p = mp.Process(target=multi_train_hdf5, args=[i])
#     p.start()
#     processes.append(p)

# for i in processes:
#     i.join()
    
# ##############################
# # Optional: Compute Statistics
# ##############################

# if args.statistics:
#     # Atomic number table
#     # yapf: disable
#     if args.atomic_numbers is None:
#         z_table = tools.get_atomic_number_table_from_zs(
#             z
#             for z in all_data_configs.atomic_numbers
#         )
#     else:
#         # logging.info("Using atomic numbers from command line argument")
#         zs_list = ast.literal_eval(args.atomic_numbers)
#         assert isinstance(zs_list, list)
#         z_table = tools.get_atomic_number_table_from_zs(zs_list)
    
#     logging.info("Computing statistics")
#     if len(atomic_energies_dict) == 0:
#         atomic_energies_dict = get_atomic_energies(args.E0s, all_data_configs, z_table)
#     atomic_energies: np.ndarray = np.array(
#         [atomic_energies_dict[z] for z in z_table.zs]
#     )
#     logging.info(f"Atomic energies: {atomic_energies.tolist()}")
#     _inputs = [args.h5_prefix+'train', z_table, args.r_max, atomic_energies, args.batch_size, args.num_process]
#     avg_num_neighbors, mean, std=pool_compute_stats(_inputs)
#     # logging.info(f"Average number of neighbors: {avg_num_neighbors}")
#     # logging.info(f"Mean: {mean}")
#     # logging.info(f"Standard deviation: {std}")

#     # save the statistics as a json
#     statistics = {
#         "atomic_energies": str(atomic_energies_dict),
#         "avg_num_neighbors": avg_num_neighbors,
#         "mean": mean,
#         "std": std,
#         "atomic_numbers": str(z_table.zs),
#         "r_max": args.r_max,
#     }
    
#     with open(args.directory + '/' + args.prefix + "statistics.json", "w") as f:
#         json.dump(statistics, f)
    