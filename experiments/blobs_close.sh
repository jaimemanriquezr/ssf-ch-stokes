#!/bin/bash
#SBATCH -n 3
#SBATCH -t 48:00:00
#SBATCH -J BLOB_CLOSE
#SBATCH -o output/slurm-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaime.manriquez@math.lth.se
#cat $0
module load buildtool-easybuild/4.8.0-hpce082752a2
module load Miniforge/24.7.1-2-hpc1
conda activate fenicsproject
SIM_NAME=blob_close
N_ELEMENTS=30
TIME_STEP=1E-07

N_STEPS=${1:-20}
TIME_SNAP=${2:-"2E-04"}

ABS_TOL=1E-8
REL_TOL=1E-7

OUTPUT_DIR="./results/"
PRESET=blobs
PRESET_PARAMETER=0.05
FILTER_LENGTH=1.0
FILTER_AREA=1.0
INLET_AREA=0.0
INLET_CENTER=0.0
INLET_LENGTH=0.0

NAME="${SIM_NAME}_N_${N_ELEMENTS}_dt_${TIME_STEP}"
function run_experiment () {
	PREFIX="${1:-}"
	C1="${2:-0.0}"
	C2="${3:-0.0}"
	ADDITIONAL_ARGS="${4:-"--is_mass_lumped"}"
	srun --ntasks=1 --nodes=1 --exclusive python main.py --name="${PREFIX}"_${NAME} \
	--output_dir=${OUTPUT_DIR} \
	--filter_length=${FILTER_LENGTH} --filter_area=${FILTER_AREA} \
	--inlet_area=${INLET_AREA} --inlet_length=${INLET_LENGTH} --inlet_center=${INLET_CENTER} \
	--preset_name=${PRESET} --preset_parameter=${PRESET_PARAMETER} \
	--time_step_snap=${TIME_SNAP} \
	--abs_tol=${ABS_TOL} --rel_tol=${REL_TOL} \
	--n_mesh=${N_ELEMENTS} --n_steps="${N_STEPS}" --time_step_size=${TIME_STEP} \
	--c_1="$C1" --c_2="$C2" "${ADDITIONAL_ARGS}" & sleep 1
}

run_experiment "zero_rx_quads" 0.0 0.0
run_experiment "norm_rx_quads" 1E2 1E3
run_experiment "high_rx_quads" 1E3 1E4

# run_experiment "zero_rx_gamma1_quads" 0.0 0.0 "-G=1"
# run_experiment "norm_rx_gamma1_quads" 1E2 1E3 "-G=1"
# run_experiment "high_rx_gamma1_quads" 1E3 1E4 "-G=1"
wait
