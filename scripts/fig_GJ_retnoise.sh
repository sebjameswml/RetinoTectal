#!/bin/sh

./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:ret_rcpt_noise_gain=0.0 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:ret_rcpt_noise_gain=0.001 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:ret_rcpt_noise_gain=0.002 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:ret_rcpt_noise_gain=0.005 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:ret_rcpt_noise_gain=0.01 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:ret_rcpt_noise_gain=0.02 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:ret_rcpt_noise_gain=0.05 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_figcomp2.json -co:steps=1000 -co:ret_rcpt_noise_gain=0.1 -co:exit=true
