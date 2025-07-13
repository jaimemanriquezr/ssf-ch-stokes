#!/bin/sh
#SBATCH -n 40
#SBATCH -t 48:00:00
#SBATCH -J SSF_2D
#SBATCH -o output/slurm-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaime.manriquez@math.lth.se
#cat $0

# Activate FEniCS environment
module load buildtool-easybuild/4.8.0-hpce082752a2
module load Miniforge/24.7.1-2-hpc1
conda activate fenicsproject

ssf="-n 1 -N 1 --exclusive python main.py"
output_dir="./results/"

# Submit jobs
N_of_steps=10
dt=1E-07
t_snap=1E-03

name=blobby
c_1=1.0
c_2=1.0

absolute_tolerance=1E-9
relative_tolerance=1E-8

inflow_velocity=1.0
inflow_concentration=2.0

length=1.0
preset=blobs
preset_param=0.05

N_list=(35 40)
for N in "${N_list[@]}"; do
    _name=${name}_N_${N}
    srun ${ssf} --name=${_name} --output_dir=${output_dir} \
        --preset_name=${preset} --preset_parameter=${preset_param} \
        --inflow_velocity=${inflow_velocity} --inflow_profile=parabolic \
        --inflow_concentration=${inflow_concentration} \
        --time_step_snap=${t_snap} \
        --abs_tol=${absolute_tolerance} --rel_tol=${relative_tolerance} \
        --c_1=${c_1} --c_2=${c_2} \
        --n_mesh=${N} --n_steps=${N_of_steps} --time_step_size=${dt} \
        &
    sleep 1
done

wait
