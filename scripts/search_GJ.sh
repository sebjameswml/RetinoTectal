#!/bin/sh

#
# Use anneal1 to search for the best parameters for the GJ model.
#
# GJ.json contains the starting model paramteres. s_GJ.json contains
# the (simulated annealing) search parameters
#
./build/sim/agent/anneal1c ./optimization/m_GJ.json ./optimization/s_GJ.json
