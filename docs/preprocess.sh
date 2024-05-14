#!/bin/bash
#SBATCH -A m4298
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 10:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 128

export SLURM_CPU_BIND="cores"
srun python /pscratch/sd/m/mavaylon/kareem_mace/mace/scripts/preprocess_data.py \
    --train_file="/pscratch/sd/m/mavaylon/chem_final_data/Traj/Charge_data/rad_qm9_65_10_25_charge_data_train.xyz" \
    --valid_file="/pscratch/sd/m/mavaylon/chem_final_data/Traj/Charge_data/rad_qm9_65_10_25_charge_data_val.xyz" \
    --test_file="/pscratch/sd/m/mavaylon/chem_final_data/Traj/Charge_data/rad_qm9_65_10_25_charge_data_test.xyz" \
    --num_process=64 \
    --atomic_numbers="[1, 6, 7, 8, 9]" \
    --total_charges="[-1, 0, 1,]" \
    --spins="[1, 2, 3]" \
    --r_max=5.0 \
    --h5_prefix="/pscratch/sd/m/mavaylon/chem_final_data/Traj/Charge_data/h5/charge_traj" \
    --seed=123 \
    --E0s="average" 