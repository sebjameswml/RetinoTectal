#!/bin/bash

# GJ model. EphA3 and EphA4 manipulations
# EphA3 ki/+ and EphA4 +/+ (i.e. EphA4 is wildtype/unmodified)
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_eph_ki-wt.json -co:exit=true -co:steps=1500
# EphA3 ki/ki and EphA4 +/+
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_eph_kiki-wt.json -co:exit=true -co:steps=1500
# EphA3 ki/+ and EphA4 +/-
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_eph_ki-kd.json -co:exit=true -co:steps=1500

# Extended GJ model (EphA cluster size). EphA3 and EphA4 manipulations
# NB: Can't use agent1c here, as I hacked the rcnt/ save out into the visualisation code.
for i in 1 2 3 4
do
    echo "Run $i"
    ./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4.json configs/simpler/e_eph_ki-wt.json -co:exit=true -co:steps=1500
    ./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4.json configs/simpler/e_eph_kiki-wt.json -co:exit=true -co:steps=1500
    ./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4.json configs/simpler/e_eph_ki-kd.json -co:exit=true -co:steps=1500
    ./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4.json configs/simpler/e_eph_wt-kd.json -co:exit=true -co:steps=1500
done
