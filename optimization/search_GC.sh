#!/bin/sh

#
# Use anneal1 to search for the best parameters for the GC model.
#
# m_ee_GC.json contains the starting model paramters. s_GC.json
# contains the (simulated annealing) search parameters
#
./build/sim/agent/anneal1c ./optimization/m_ee_GC.json ./optimization/s_GC.json
