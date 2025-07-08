#!/bin/bash
if [ $(whoami) == *"tetralith"* ]; then
    module load buildtool-easybuild/4.8.0-hpce082752a2
    module load ParaView/5.13.1-hpc1-bdist
fi
DIR=${1:-.}

MAXIMUM_PARTICLES=25.0
MAXIMUM_LIQUIDS=998.0
FONT_SIZE=36
TICKS=3
FINAL_TIME=1E-2
N_STEPS=10

for filename in $DIR/*.xdmf; do
    name=${filename#"$DIR/"}
    name=${name%".xdmf"}
    pvpython ssf_figure_export.py ${name} \
        --maximum_particle_concentration=${MAXIMUM_PARTICLES} \
        --maximum_fluid_concentration=${MAXIMUM_LIQUIDS} \
        --axes_font=${FONT_SIZE} \
        --axes_ticks=${TICKS} \
        --final_time=${FINAL_TIME} \
        --time_steps=${N_STEPS} \
        --output=pdf \
        --all \
        # --compress
done
