#!/bin/sh
#SBATCH -n 6
#SBATCH -t 24:00:00
#SBATCH -J BLOB_CLOSE
#SBATCH -o output/slurm-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaime.manriquez@math.lth.se
#cat $0

NAME=blob_close
OUTPUT_DIR="./results/"
PRESET=blobs
PRESET_PARAMETER=0.05

FILTER_LENGTH=1.0
FILTER_AREA=1.0

INLET_AREA=0.0
INLET_CENTER=0.0
INLET_LENGTH=0.0

N_ELEMENTS=30
TIME_STEP=1E-07
N_STEPS=10000
TIME_SNAP=1E-05

ABS_TOL=1E-9
REL_TOL=1E-8

# Activate FEniCS environment
module load buildtool-easybuild/4.8.0-hpce082752a2
module load Miniforge/24.7.1-2-hpc1
conda activate fenicsproject

srun --ntasks=1 --nodes=1 --exclusive python main.py --name=zero_rx_${NAME} \
--output_dir=${OUTPUT_DIR} \
--filter_length=${FILTER_LENGTH} --filter_area=${FILTER_AREA} \
--inlet_area=${INLET_AREA} --inlet_length=${INLET_LENGTH} --inlet_center=${INLET_CENTER} \
--preset_name=${PRESET} --preset_parameter=${PRESET_PARAMETER} \
--time_step_snap=${TIME_SNAP} \
--abs_tol=${ABS_TOL} --rel_tol=${REL_TOL} \
--n_mesh=${N_ELEMENTS} --n_steps=${N_STEPS} --time_step_size=${TIME_STEP} \
--c_1=0.0 --c_2=0.0 & sleep 1

srun --ntasks=1 --nodes=1 --exclusive python main.py --name=norm_rx_${NAME} \
--output_dir=${OUTPUT_DIR} \
--filter_length=${FILTER_LENGTH} --filter_area=${FILTER_AREA} \
--inlet_area=${INLET_AREA} --inlet_length=${INLET_LENGTH} --inlet_center=${INLET_CENTER} \
--preset_name=${PRESET} --preset_parameter=${PRESET_PARAMETER} \
--time_step_snap=${TIME_SNAP} \
--abs_tol=${ABS_TOL} --rel_tol=${REL_TOL} \
--n_mesh=${N_ELEMENTS} --n_steps=${N_STEPS} --time_step_size=${TIME_STEP} \
--c_1=1E2 --c_2=1E3 & sleep 1

srun --ntasks=1 --nodes=1 --exclusive python main.py --name=high_rx_${NAME} \
--output_dir=${OUTPUT_DIR} \
--filter_length=${FILTER_LENGTH} --filter_area=${FILTER_AREA} \
--inlet_area=${INLET_AREA} --inlet_length=${INLET_LENGTH} --inlet_center=${INLET_CENTER} \
--preset_name=${PRESET} --preset_parameter=${PRESET_PARAMETER} \
--time_step_snap=${TIME_SNAP} \
--abs_tol=${ABS_TOL} --rel_tol=${REL_TOL} \
--n_mesh=${N_ELEMENTS} --n_steps=${N_STEPS} --time_step_size=${TIME_STEP} \
--c_1=1E3 --c_2=1E4 & sleep 1

srun --ntasks=1 --nodes=1 --exclusive python main.py --name=zero_rx_${NAME}_gamma1 \
--output_dir=${OUTPUT_DIR} \
--filter_length=${FILTER_LENGTH} --filter_area=${FILTER_AREA} \
--inlet_area=${INLET_AREA} --inlet_length=${INLET_LENGTH} --inlet_center=${INLET_CENTER} \
--preset_name=${PRESET} --preset_parameter=${PRESET_PARAMETER} \
--time_step_snap=${TIME_SNAP} \
--abs_tol=${ABS_TOL} --rel_tol=${REL_TOL} \
--n_mesh=${N_ELEMENTS} --n_steps=${N_STEPS} --time_step_size=${TIME_STEP} \
-G=1 --c_1=0.0 --c_2=0.0 & sleep 1

srun --ntasks=1 --nodes=1 --exclusive python main.py --name=norm_rx_${NAME}_gamma1 \
--output_dir=${OUTPUT_DIR} \
--filter_length=${FILTER_LENGTH} --filter_area=${FILTER_AREA} \
--inlet_area=${INLET_AREA} --inlet_length=${INLET_LENGTH} --inlet_center=${INLET_CENTER} \
--preset_name=${PRESET} --preset_parameter=${PRESET_PARAMETER} \
--time_step_snap=${TIME_SNAP} \
--abs_tol=${ABS_TOL} --rel_tol=${REL_TOL} \
--n_mesh=${N_ELEMENTS} --n_steps=${N_STEPS} --time_step_size=${TIME_STEP} \
-G=1 --c_1=1E2 --c_2=1E3 & sleep 1

srun --ntasks=1 --nodes=1 --exclusive python main.py --name=high_rx_${NAME}_gamma1 \
--output_dir=${OUTPUT_DIR} \
--filter_length=${FILTER_LENGTH} --filter_area=${FILTER_AREA} \
--inlet_area=${INLET_AREA} --inlet_length=${INLET_LENGTH} --inlet_center=${INLET_CENTER} \
--preset_name=${PRESET} --preset_parameter=${PRESET_PARAMETER} \
--time_step_snap=${TIME_SNAP} \
--abs_tol=${ABS_TOL} --rel_tol=${REL_TOL} \
--n_mesh=${N_ELEMENTS} --n_steps=${N_STEPS} --time_step_size=${TIME_STEP} \
-G=1 --c_1=1E3 --c_2=1E4 & sleep 1
wait
