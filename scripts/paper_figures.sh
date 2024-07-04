#!/bin/bash

#
# Generate all the paper figure graph elements, one after the other.
#

# First clear out rcnt
rm -rf ./rcnt/*.h5

# Tissue visualisation (no longer generated for Fig 1)
# ./scripts/fig_expression.sh

# Chemoaffinity model, gradient expression from literature. WT.
./scripts/fig_G0_wildtype.sh # possibly included as supplementary figure

# GJ model, wildtype
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_fig2.json -co:exit=true -co:steps=1500
# GJ model, Surgical manipulations
./scripts/fig_GJ_surgical.sh
# GJ model tectal rotation visualisation for supplementary figure
./build/sim/agent/tissuevis configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_tecrot90_suppfig.json

# GJ model. EphA3 and EphA4 manipulations
# 1. EphA3 ki/+ and EphA4 +/+ (i.e. EphA4 is wildtype/unmodified)
# 2. EphA3 ki/ki and EphA4 +/+
# 3. EphA3 ki/+ and EphA4 +/-
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_eph_ki-wt.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_eph_kiki-wt.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_eph_ki-kd.json -co:exit=true -co:steps=1500

# Extended GJ model (EphA cluster size) plus the r2 collapse condition (EphA3 and EphA4 manipulations)
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4_r2collapse.json configs/simpler/e_eph_ki-wt.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4_r2collapse.json configs/simpler/e_eph_kiki-wt.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4_r2collapse.json configs/simpler/e_eph_ki-kd.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4_r2collapse.json configs/simpler/e_eph_wt-kd.json -co:exit=true -co:steps=1500

# Single axons
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_single-centre.json -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_single-caudal.json -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_single-rostral.json -co:exit=true
