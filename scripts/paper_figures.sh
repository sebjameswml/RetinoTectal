#!/bin/bash

#
# Generate all the paper figures, one after the other.
#

# Tissue visualisation
./scripts/fig_expression.sh

# Chemoaffinity model, gradient expression from literature. WT.
./scripts/fig_G0_wildtype.sh
# Chemoaffinity model, hand-tuned gradient expression. WT.
./scripts/fig_G_wildtype.sh
# Chemoaffinity model, hand-tuned gradient expression. Surgical manipulations
./scripts/fig_G_surgical.sh

# GC model, wildtype
./build/sim/agent/agent1 configs/simpler/m_ee_GC_best_1.json configs/simpler/e_wt_figcomp1.json -co:exit=true
# GC model, Surgical manipulations
./scripts/fig_GC_surgical.sh

# GI model, Surgical manipulations (not included in paper)
#./build/sim/agent/agent1 configs/simpler/m_ee_GI_best_1.json configs/simpler/e_wt_figcomp2.json -co:exit=true
#./scripts/fig_GI_surgical.sh

# GJ model, wildtype
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_wt_figcomp3.json -co:exit=true -co:steps=1500
# GJ model, Surgical manipulations
./scripts/fig_GJ_surgical.sh

# GJ model. EphA3 and EphA4 manipulations
# EphA3 ki/+ and EphA4 +/+ (i.e. EphA4 is wildtype/unmodified)
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_eph_ki-wt.json -co:exit=true -co:steps=1500
# EphA3 ki/ki and EphA4 +/+
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_eph_kiki-wt.json -co:exit=true -co:steps=1500
# EphA3 ki/+ and EphA4 +/-
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_eph_ki-kd.json -co:exit=true -co:steps=1500

# Extended GJ model (EphA cluster size). EphA3 and EphA4 manipulations
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4.json configs/simpler/e_eph_ki-wt.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4.json configs/simpler/e_eph_kiki-wt.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_EphA4.json configs/simpler/e_eph_ki-kd.json -co:exit=true -co:steps=1500

# Compare GC model
./build/sim/agent/agent1 configs/simpler/m_ee_GC_best_1.json configs/simpler/e_eph_ki-wt.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GC_best_1.json configs/simpler/e_eph_kiki-wt.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GC_best_1.json configs/simpler/e_eph_ki-kd.json -co:exit=true -co:steps=1500

# Compare GC model with same paramters as GJ model:
./build/sim/agent/agent1 configs/simpler/m_ee_GClikeGJ.json configs/simpler/e_eph_ki-wt.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GClikeGJ.json configs/simpler/e_eph_kiki-wt.json -co:exit=true -co:steps=1500
./build/sim/agent/agent1 configs/simpler/m_ee_GClikeGJ.json configs/simpler/e_eph_ki-kd.json -co:exit=true -co:steps=1500

# Single axons
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_single-centre.json -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_single-caudal.json -co:exit=true
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1.json configs/simpler/e_single-rostral.json -co:exit=true
