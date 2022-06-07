#!/bin/bash

# GJ model. EphA3 and EphA4 manipulations
# EphA3 ki/+ and EphA4 +/+ (i.e. EphA4 is wildtype/unmodified)
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_tec_ligand_exp3.json configs/simpler/e_eph_ki-wt.json -co:exit=false -co:steps=1500

# EphA3 ki/ki and EphA4 +/+
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_tec_ligand_exp3.json configs/simpler/e_eph_kiki-wt.json -co:exit=false -co:steps=1500

# EphA3 ki/+ and EphA4 +/-
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_tec_ligand_exp3.json configs/simpler/e_eph_ki-kd.json -co:exit=false -co:steps=1500


# EphA3 ki/+ and EphA4 +/+ (i.e. EphA4 is wildtype/unmodified)
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_tec_ligand_exp4.json configs/simpler/e_eph_ki-wt.json -co:exit=false -co:steps=5000
# EphA3 ki/ki and EphA4 +/+
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_tec_ligand_exp4.json configs/simpler/e_eph_kiki-wt.json -co:exit=false -co:steps=5000
# EphA3 ki/+ and EphA4 +/-
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_tec_ligand_exp4.json configs/simpler/e_eph_ki-kd.json -co:exit=false -co:steps=5000



# EphA3 ki/+ and EphA4 +/+ (i.e. EphA4 is wildtype/unmodified)
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_tec_ligand_sigmoid.json configs/simpler/e_eph_ki-wt.json -co:exit=false -co:steps=5000
# EphA3 ki/ki and EphA4 +/+
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_tec_ligand_sigmoid.json configs/simpler/e_eph_kiki-wt.json -co:exit=false -co:steps=5000
# EphA3 ki/+ and EphA4 +/-
./build/sim/agent/agent1 configs/simpler/m_ee_GJ_best_1_tec_ligand_sigmoid.json configs/simpler/e_eph_ki-kd.json -co:exit=false -co:steps=5000
