#!/bin/sh

# Rotational wildtype, but with different selected axons:
./build/sim/agent/agent1 ./configs/a1/m_eE_GJ.json ./configs/a1/e_retablate.json
./build/sim/agent/agent1 ./configs/a1/m_eE_GJ.json ./configs/a1/e_tecablate.json
./build/sim/agent/agent1 ./configs/a1/m_eE_GJ.json ./configs/a1/e_mismatch.json
./build/sim/agent/agent1 ./configs/a1/m_eE_GJ.json ./configs/a1/e_compound.json
