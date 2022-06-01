#!/bin/sh

#
# Use anneal1 to search for the best parameters for the GJ model with
# the intermediate exponential function, exp3.
#
# m_e3_GJ.json contains the starting model paramters. s_GJ.json
# contains the (simulated annealing) search parameters
#
./build/sim/agent/anneal1c ./optimization/m_e3_GJ.json ./optimization/s_GJ.json
