#!/usr/bin/python

import sys
import json
import subprocess

# Get a json file from cmd, if possible, falling back to runset.json if not
jsonfile = 'runset.json'
if len(sys.argv) > 1:
    jsonfile = sys.argv[1]

# Run a set of sims with models and expts specified in json
jfh = open(jsonfile)
data = json.load(jfh)

for m in data['models']:
    for e in data['expts']:
        print ("Running Model {0}, Expt {1}".format(m, e))
        subprocess.call(["../build/sim/agent/agent1c",
                         "../configs/a1/m_{0}.json".format(m),
                         "../configs/a1/e_{0}.json".format(e)])
