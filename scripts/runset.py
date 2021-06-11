#!/usr/bin/python

import json
import subprocess

# Run a set of sims with models and expts specified in json
jf = open('runset.json')
data = json.load(jf)

for m in data['models']:
    for e in data['expts']:
        print ("Running Model {0}, Expt {1}".format(m, e))
        subprocess.call(["../build/sim/agent/agent1c",
                         "../configs/a1/m_{0}.json".format(m),
                         "../configs/a1/e_{0}.json".format(e)])
