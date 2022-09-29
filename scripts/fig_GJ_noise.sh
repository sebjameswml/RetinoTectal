#!/bin/sh

./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fignoise.json -co:steps=2000 -co:noise_gain=0.19 -co:ret_rcpt_noise_gain=0.0  -co:tec_lgnd_noise_gain=0.0 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fignoise.json -co:steps=2000 -co:noise_gain=0.0  -co:ret_rcpt_noise_gain=0.08 -co:tec_lgnd_noise_gain=0.0 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fignoise.json -co:steps=2000 -co:noise_gain=0.0  -co:ret_rcpt_noise_gain=0.0  -co:tec_lgnd_noise_gain=0.08 -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fignoise.json -co:steps=2000 -co:noise_gain=0.19 -co:ret_rcpt_noise_gain=0.08 -co:tec_lgnd_noise_gain=0.08 -co:exit=true

./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fignoise.json -co:steps=2000 -co:noise_gain=1 -co:ret_rcpt_noise_gain=0.0  -co:tec_lgnd_noise_gain=0.0 -co:exit=true

#./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fignoise.json -co:steps=1000 -co:noise_gain=0.2 -co:exit=true
#./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fignoise.json -co:steps=1000 -co:noise_gain=0.4 -co:exit=true
#./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fignoise.json -co:steps=1000 -co:noise_gain=0.6 -co:exit=true
#./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fignoise.json -co:steps=1000 -co:noise_gain=0.8 -co:exit=true
#./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fignoise.json -co:steps=1000 -co:noise_gain=1.0 -co:exit=true
#./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fignoise.json -co:steps=1000 -co:noise_gain=1.5 -co:exit=true
#./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fignoise.json -co:steps=1000 -co:noise_gain=2.0 -co:exit=true
