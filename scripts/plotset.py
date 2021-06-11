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
fig = pl.figure(figsize=(8,10))

import h5py

models = ['explog', 'linear', 'expexp']
expts = ['wildtype', 'tecswap', 'tecrot90', 'tecrot180',
         'tecablate', 'single', 'retablate',
         'reber', 'mismatch', 'knockout1',
         'knockin1', 'compound']

# a models (rows) by expts (cols) array for results
r = np.zeros((len(models),len(expts)))

for mi in range(0,len(models)):
    for ei in range(0,len(expts)):
        filename = '../log/agent/'+models[mi]+'_'+expts[ei]+'.h5'
        with h5py.File (filename, 'r') as f:
            r[mi][ei] = list(f['/sos'])[0]

print ('r: {0}'.format(r))

ax = fig.add_subplot (1,1,1);
im = ax.imshow (r)

#pl.tight_layout()
#pl.savefig('somefile.svg', transparent=True)
pl.show()
