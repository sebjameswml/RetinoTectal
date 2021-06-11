import numpy as np

import sebcolour
C = sebcolour.Colour

# Set plotting font defaults
import matplotlib
fs = 18
fnt = {'family' : 'Arial',
       'weight' : 'regular',
       'size'   : fs}
matplotlib.rc('font', **fnt)

# Important for svg output of text as 'things that can be edited in inkscape'
import pylab as pl
pl.rcParams['svg.fonttype'] = 'none'

# Control overall figure size here
fig = pl.figure(figsize=(12,5))

import h5py
import json

models = ['explog', 'linear', 'expexp']
expts = ['wildtype', 'tecswap', 'tecrot90', 'tecrot180',
         'tecablate', 'single', 'retablate',
         'reber', 'mismatch', 'knockout1',
         'knockin1', 'compound']

jf = open('runset.json')
jconf = json.load(jf)
models = jconf['models']
expts = jconf['expts']

# a models (rows) by expts (cols) array for results
r = np.zeros((len(models),len(expts)))

for mi in range(0,len(models)):
    for ei in range(0,len(expts)):
        filename = '../log/agent/'+models[mi]+'_'+expts[ei]+'.h5'
        with h5py.File (filename, 'r') as f:
            r[mi][ei] = list(f['/sos'])[0]

print ('r: {0}'.format(r))

# Normalise each col of r to show best model for each condition?
r_n = r / r.max(axis=0)

ax = fig.add_subplot (1,1,1);
im = ax.imshow (r, cmap='plasma')
ax.set_xticks(range(0, len(expts)))
ax.set_yticks(range(0, len(models)))
ax.set_xticklabels (expts)
ax.set_yticklabels (models)
pl.xticks(rotation=90)
pl.colorbar(im)
ax.set_title('SOS pattern match. Purple is best. Raw.')


#pl.tight_layout()
#pl.savefig('somefile.svg', transparent=True)
pl.show()
