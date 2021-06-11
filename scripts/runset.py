#!/usr/bin/python

import subprocess

# Run a set of sims
models = ['el', 'lin', 'ee']
expts = ['wt', 'tecswap', 'tecrot90', 'tecrot180',
         'tecablate', 'single', 'retablate',
         'reber', 'mismatch', 'knockout1',
         'knockin1', 'compound']

for m in models:
    for e in expts:
        print ("Running Model {0}, Expt {1}".format(m, e))
        subprocess.call(["./build/sim/agent/agent1c",
                         "configs/a1/m_{0}.json".format(m),
                         "configs/a1/e_{0}.json".format(e)])
