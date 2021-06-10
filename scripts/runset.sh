#!/bin/bash

# Run a set of sims for one model

MDL=configs/a1/m_el.json

for EXPT in e_wt.json e_tecswap.json e_tecrot90.json e_tecrot180.json e_tecablate.json e_single.json e_retablate.json e_reber.json e_mismatch.json e_knockout1.json e_knockin1.json e_compound.json; do
    echo "Model $MDL, Experiment $EXPT..."
    ./build/sim/agent/agent1c $MDL configs/a1/$EXPT
done
