#!/bin/sh

./build/sim/agent/agent1 configs/simpler/m_eE_G.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:noise_gain=0.0 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_eE_G.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:noise_gain=0.2 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_eE_G.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:noise_gain=0.4 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_eE_G.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:noise_gain=0.6 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_eE_G.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:noise_gain=0.8 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_eE_G.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:noise_gain=1.0 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_eE_G.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:noise_gain=1.5 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_eE_G.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:noise_gain=2.0 -co:exit=true