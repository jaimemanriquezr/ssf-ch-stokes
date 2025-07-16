#!/bin/bash
#SBATCH -n 24
#SBATCH -t 48:00:00
#SBATCH -J BLOB_FAR
#SBATCH -o output/slurm-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaime.manriquez@math.lth.se
#cat $0
module load buildtool-easybuild/4.8.0-hpce082752a2
module load Miniforge/24.7.1-2-hpc1
conda activate fenicsproject
SIM_NAME=blob_far
N_ELEMENTS=30
TIME_STEP=1E-07

N_STEPS=${1:-20}
TIME_SNAP=2E-04

ABS_TOL=1E-8
REL_TOL=1E-7

OUTPUT_DIR="./results/"
PRESET=blobs
PRESET_PARAMETER=0.15
FILTER_LENGTH=1.0
FILTER_AREA=1.0
INLET_AREA=0.0
INLET_CENTER=0.0
INLET_LENGTH=0.0

function run_experiment () {
	PREFIX="${1:-}"
	N_ELEMENTS="${2:-}"
	ABS_TOL="${3:-}"
	REL_TOL="${4:-}"
	C1="${5:-0.0}"
	C2="${6:-0.0}"
	ADDITIONAL_ARGS="${7:-"--is_mass_lumped"}"
	NAME="${SIM_NAME}_N_${N_ELEMENTS}_dt_${TIME_STEP}"
	srun --ntasks=1 --nodes=1 --exclusive python main.py --name="${PREFIX}"_"${NAME}" \
	--output_dir=${OUTPUT_DIR} \
	--filter_length=${FILTER_LENGTH} --filter_area=${FILTER_AREA} \
	--inlet_area=${INLET_AREA} --inlet_length=${INLET_LENGTH} --inlet_center=${INLET_CENTER} \
	--preset_name=${PRESET} --preset_parameter=${PRESET_PARAMETER} \
	--time_step_snap=${TIME_SNAP} \
	--abs_tol="${ABS_TOL}" --rel_tol="${REL_TOL}" \
	--n_mesh="${N_ELEMENTS}" --n_steps="${N_STEPS}" --time_step_size=${TIME_STEP} \
	--c_1="$C1" --c_2="$C2" "${ADDITIONAL_ARGS}" & sleep 1
}

run_experiment "zero_rx_TOL9" 50 1E-9 1E-8 0.0 0.0
run_experiment "norm_rx_TOL9" 50 1E-9 1E-8 1E2 1E3
run_experiment "high_rx_TOL9" 50 1E-9 1E-8 1E3 1E4

run_experiment "zero_rx_TOL9" 40 1E-9 1E-8 0.0 0.0
run_experiment "norm_rx_TOL9" 40 1E-9 1E-8 1E2 1E3
run_experiment "high_rx_TOL9" 40 1E-9 1E-8 1E3 1E4

run_experiment "zero_rx_TOL9" 35 1E-9 1E-8 0.0 0.0
run_experiment "norm_rx_TOL9" 35 1E-9 1E-8 1E2 1E3
run_experiment "high_rx_TOL9" 35 1E-9 1E-8 1E3 1E4

run_experiment "zero_rx_TOL9" 30 1E-9 1E-8 0.0 0.0
run_experiment "norm_rx_TOL9" 30 1E-9 1E-8 1E2 1E3
run_experiment "high_rx_TOL9" 30 1E-9 1E-8 1E3 1E4

run_experiment "zero_rx_TOL8" 30 1E-8 1E-7 0.0 0.0
run_experiment "norm_rx_TOL8" 30 1E-8 1E-7 1E2 1E3
run_experiment "high_rx_TOL8" 30 1E-8 1E-7 1E3 1E4

run_experiment "zero_rx_TOL7" 30 1E-7 1E-6 0.0 0.0
run_experiment "norm_rx_TOL7" 30 1E-7 1E-6 1E2 1E3
run_experiment "high_rx_TOL7" 30 1E-7 1E-6 1E3 1E4

run_experiment "zero_rx_TOL9" 25 1E-9 1E-8 0.0 0.0
run_experiment "norm_rx_TOL9" 25 1E-9 1E-8 1E2 1E3
run_experiment "high_rx_TOL9" 25 1E-9 1E-8 1E3 1E4

run_experiment "zero_rx_TOL9_gamma1" 30 1E-9 1E-8 0.0 0.0 "-G=1"
run_experiment "norm_rx_TOL9_gamma1" 30 1E-9 1E-8 1E2 1E3 "-G=1"
run_experiment "high_rx_TOL9_gamma1" 30 1E-9 1E-8 1E3 1E4 "-G=1"
wait
