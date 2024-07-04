#!/bin/bash

#
# Generate all the paper figure graph elements, one after the other.
#

# First clear out some of rcnt
rm -f ./rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_ki-kd_exit_true_steps_1500_*.h5
rm -f ./rcnt/j4_ee_GJ_best_1_eph_ki-kd_exit_true_steps_1500_*.h5

# GJ model. EphA3 and EphA4 manipulations
# 1. EphA3 ki/+ and EphA4 +/+ (i.e. EphA4 is wildtype/unmodified)
# 2. EphA3 ki/ki and EphA4 +/+
# 3. EphA3 ki/+ and EphA4 +/-

#./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_eph_ki-wt.json -co:exit=true -co:steps=1500
#./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_eph_kiki-wt.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_eph_ki-kd.json -co:exit=true -co:steps=1500

# Extended GJ model (EphA cluster size) plus the r2 collapse condition (EphA3 and EphA4 manipulations)
#./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4_r2collapse.json configs/simpler/e_eph_ki-wt.json -co:exit=true -co:steps=1500
#./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4_r2collapse.json configs/simpler/e_eph_kiki-wt.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4_r2collapse.json configs/simpler/e_eph_ki-kd.json -co:exit=true -co:steps=1500
#./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4_r2collapse.json configs/simpler/e_eph_wt-kd.json -co:exit=true -co:steps=1500
