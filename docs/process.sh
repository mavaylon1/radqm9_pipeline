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
    --file="/pscratch/sd/m/mavaylon/sam_ldrd/all_clean_charge_full_chunked/charge_spin_test/charge_spin_test_0_3" \
    --forces_key="forces" \
    --prefix="charge_spin_test_0_3_" \
    --directory="/pscratch/sd/m/mavaylon/chem_directory/charge_spin_test_0_3" \

