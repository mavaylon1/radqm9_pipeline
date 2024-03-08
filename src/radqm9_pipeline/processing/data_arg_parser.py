import argparse
import os

def build_default_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "--file",
        help="XYZ file.",
        type=str,
        required = True
    )
    parser.add_argument(
        "--prefix",
        help="The prefix for the h5 file.",
        type=str,
        required = True
    )
    parser.add_argument(
        "--atomic_numbers",
        help="List of atomic numbers",
        type=str,
        default=None,
        required=False,
    )
    parser.add_argument(
        "--num_process",
        help="The user defined number of processes to use(the number of files created).", 
        type=int, 
        default=int(os.cpu_count()/4)
    )
    parser.add_argument(
        "--E0s",
        help="Dictionary of isolated atom energies",
        type=str,
        default=None,
        required=False,
    )
    parser.add_argument(
        "--r_max",
        help="distance cutoff (in Ang)", 
        type=float, 
        default=5.0
    )
    parser.add_argument(
        "--batch_size", 
        help="batch size", 
        type=int, 
        default=10)
    parser.add_argument(
        "--config_type_weights",
        help="String of dictionary containing the weights for each config type",
        type=str,
        default='{"Default":1.0}',
    )
    parser.add_argument(
        "--energy_key",
        help="Key of reference energies in training xyz",
        type=str,
        default="energy",
    )
    parser.add_argument(
        "--forces_key",
        help="Key of reference forces in training xyz",
        type=str,
        default="forces",
    )
    parser.add_argument(
        "--virials_key",
        help="Key of reference virials in training xyz",
        type=str,
        default="virials",
    )
    parser.add_argument(
        "--stress_key",
        help="Key of reference stress in training xyz",
        type=str,
        default="stress",
    )
    parser.add_argument(
        "--dipole_key",
        help="Key of reference dipoles in training xyz",
        type=str,
        default="dipole",
    )
    parser.add_argument(
        "--charges_key",
        help="Key of atomic charges in training xyz",
        type=str,
        default="charges",
    )
    parser.add_argument(
        "--directory",
        help="The directory to store the h5 files.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--statistics",
        help="Bool to trigger statistics.json",
        type=bool,
        default=False
    )
   
    return parser

        
        