#!/bin/bash
if [[ "$(whoami)" == *"tetralith"* ]]; then
    module load buildtool-easybuild/4.8.0-hpce082752a2
    module load ParaView/5.13.1-hpc1-bdist
fi
DIR=${1:-./experiments/}

for filename in "$DIR"*.sh; do
    echo "Batching ${filename}"
    sbatch "${filename}"
done
