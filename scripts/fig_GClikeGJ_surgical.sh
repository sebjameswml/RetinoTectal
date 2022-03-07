#!/bin/sh

./build/sim/agent/agent1 configs/simpler/m_ee_GClikeGJ.json  configs/simpler/e_tecrot90_figcomp.json  -co:steps=1000
./build/sim/agent/agent1 configs/simpler/m_ee_GClikeGJ.json  configs/simpler/e_tecrot180_figcomp.json -co:steps=1000
./build/sim/agent/agent1 configs/simpler/m_ee_GClikeGJ.json  configs/simpler/e_tecswap_figcomp.json   -co:steps=1000
./build/sim/agent/agent1 configs/simpler/m_ee_GClikeGJ.json  configs/simpler/e_retablat_figcomp.json  -co:steps=1000
./build/sim/agent/agent1 configs/simpler/m_ee_GClikeGJ.json  configs/simpler/e_tecablat_figcomp.json  -co:steps=1000
./build/sim/agent/agent1 configs/simpler/m_ee_GClikeGJ.json  configs/simpler/e_mismatch_figcomp.json  -co:steps=1000
