#!/bin/sh

./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json  configs/simpler/e_tecrot90_figcomp.json  -co:steps=1000
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json  configs/simpler/e_tecrot180_figcomp.json -co:steps=1000
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json  configs/simpler/e_tecswap_figcomp.json   -co:steps=1000
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json  configs/simpler/e_retablat_figcomp.json  -co:steps=1000
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json  configs/simpler/e_tecablat_figcomp.json  -co:steps=1000
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json  configs/simpler/e_mismatch_figcomp.json  -co:steps=1000

#./build/sim/agent/tissuevis configs/simpler/m_ee_GJ_best_1.json  configs/simpler/e_tecrot90_fig4row1.json  1000
#./build/sim/agent/tissuevis configs/simpler/m_ee_GJ_best_1.json  configs/simpler/e_tecrot180_fig4row2.json  1000
#./build/sim/agent/tissuevis configs/simpler/m_ee_GJ_best_1.json  configs/simpler/e_tecswap_fig4row3.json  1000
#./build/sim/agent/tissuevis configs/simpler/m_ee_GJ_best_1.json  configs/simpler/e_retablat_fig4row4.json  1000
#./build/sim/agent/tissuevis configs/simpler/m_ee_GJ_best_1.json  configs/simpler/e_tecablat_fig4row5.json  1000
#./build/sim/agent/tissuevis configs/simpler/m_ee_GJ_best_1.json  configs/simpler/e_mismatch_fig4row6.json  1000
