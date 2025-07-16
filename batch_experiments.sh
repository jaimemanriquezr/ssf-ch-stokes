#!/bin/bash
if [[ "$(whoami)" == *"tetralith"* ]]; then
    module load buildtool-easybuild/4.8.0-hpce082752a2
    module load ParaView/5.13.1-hpc1-bdist
fi
N_STEPS=${1:-10}
TIME_SNAP=${2:-"2E-04"}
DIR=${3:-./experiments/}

for filename in "$DIR"*.sh; do
    echo "Batching ${filename}"
    sbatch "${filename}" "${N_STEPS}" "${TIME_SNAP}"
done
