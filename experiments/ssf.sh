#!/bin/sh
#SBATCH -n 20
#SBATCH -t 48:00:00
#SBATCH -J SSF_2D
#SBATCH -o output/slurm-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaime.manriquez@math.lth.se
#cat $0

NAME=filter
OUTPUT_DIR="./results/"
PRESET=filter
PRESET_PARAMETER=0.3 # 0.6 * 0.5

FILTER_LENGTH=0.5
FILTER_AREA=0.5

INLET_AREA=0.1
INLET_CENTER=0.1
INLET_LENGTH=0.025

N_ELEMENTS=30
TIME_STEP=1E-07
N_STEPS=10
TIME_SNAP=1E-03

ABS_TOL=1E-9
REL_TOL=1E-8

Q_MAX_INFLOW=300.0
S_INFLOW=988.02 # 0.99 * 998

# Activate FEniCS environment
module load buildtool-easybuild/4.8.0-hpce082752a2
module load Miniforge/24.7.1-2-hpc1
conda activate fenicsproject
srun --ntasks=1 --nodes=1 --exclusive python main.py --name=zero_rx_${NAME} \
--output_dir=${OUTPUT_DIR} \
--filter_length=${FILTER_LENGTH} --filter_area=${FILTER_AREA} \
--inlet_area=${INLET_AREA} --inlet_length=${INLET_LENGTH} --inlet_center=${INLET_CENTER} \
--preset_name=${PRESET} --preset_parameter=${PRESET_PARAMETER} \
--inflow_velocity=${Q_MAX_INFLOW} --inflow_profile=parabolic --inflow_concentration=${S_INFLOW} \
--time_step_snap=${TIME_SNAP} \
--abs_tol=${ABS_TOL} --rel_tol=${REL_TOL} \
--n_mesh=${N_ELEMENTS} --n_steps=${N_STEPS} --time_step_size=${TIME_STEP} \
--c_1=0.0 --c_2=0.0 & sleep 1

srun --ntasks=1 --nodes=1 --exclusive python main.py --name=norm_rx_${NAME} \
--output_dir=${OUTPUT_DIR} \
--filter_length=${FILTER_LENGTH} --filter_area=${FILTER_AREA} \
--inlet_area=${INLET_AREA} --inlet_length=${INLET_LENGTH} --inlet_center=${INLET_CENTER} \
--preset_name=${PRESET} --preset_parameter=${PRESET_PARAMETER} \
--inflow_velocity=${Q_MAX_INFLOW} --inflow_profile=parabolic --inflow_concentration=${S_INFLOW} \
--time_step_snap=${TIME_SNAP} \
--abs_tol=${ABS_TOL} --rel_tol=${REL_TOL} \
--n_mesh=${N_ELEMENTS} --n_steps=${N_STEPS} --time_step_size=${TIME_STEP} \
--c_1=2E1 --c_2=1E3 & sleep 1

srun --ntasks=1 --nodes=1 --exclusive python main.py --name=high_rx_${NAME} \
--output_dir=${OUTPUT_DIR} \
--filter_length=${FILTER_LENGTH} --filter_area=${FILTER_AREA} \
--inlet_area=${INLET_AREA} --inlet_length=${INLET_LENGTH} --inlet_center=${INLET_CENTER} \
--preset_name=${PRESET} --preset_parameter=${PRESET_PARAMETER} \
--inflow_velocity=${Q_MAX_INFLOW} --inflow_profile=parabolic --inflow_concentration=${S_INFLOW} \
--time_step_snap=${TIME_SNAP} \
--abs_tol=${ABS_TOL} --rel_tol=${REL_TOL} \
--n_mesh=${N_ELEMENTS} --n_steps=${N_STEPS} --time_step_size=${TIME_STEP} \
--c_1=1E3 --c_2=5E3 & sleep 1
wait
