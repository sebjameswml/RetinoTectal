#!/bin/sh

# Rotational wildtype, but with different selected axons:
./build/sim/agent/agent1 ./configs/a1/m_eE_GC.json ./configs/a1/e_wt_rotcmp.json
./build/sim/agent/agent1 ./configs/a1/m_eE_GC.json ./configs/a1/e_tecrot90.json
./build/sim/agent/agent1 ./configs/a1/m_eE_GC.json ./configs/a1/e_tecrot180.json
./build/sim/agent/agent1 ./configs/a1/m_eE_GC.json ./configs/a1/e_tecswap.json
