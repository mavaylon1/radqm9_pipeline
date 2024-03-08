#!/bin/bash
#SBATCH -A m4298
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 5:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 128

export SLURM_CPU_BIND="cores"
srun python ../src/radqm9_pipeline/processing/process_dataset.py \
#     --file="path to dataset " \
#     --forces_key="forces" \
#     --prefix="<add prefix>" \
#     --directory='path to directory to store hdf5 files. Do not end with "/"' \

