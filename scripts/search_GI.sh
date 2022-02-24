#!/bin/sh

#
# Use anneal1 to search for the best parameters for the GI model.
#
# m_ee_GI.json contains the starting model paramters. s_GI.json
# contains the (simulated annealing) search parameters
#
./build/sim/agent/anneal1c ./optimization/m_ee_GI.json ./optimization/s_GI.json
