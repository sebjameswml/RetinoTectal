import numpy as np
import sebcolour
C = sebcolour.Colour

import matplotlib
fs = 18
fnt = {'family' : 'Arial',
       'weight' : 'regular',
       'size'   : fs}
matplotlib.rc('font', **fnt) # set default font

import pylab as pl
pl.rcParams['svg.fonttype'] = 'none' # for good svg output
fig = pl.figure(figsize=(12,5))

import json
jf = open('runset.json')
jconf = json.load(jf)
models = jconf['models']
expts = jconf['expts']

# a models (rows) by expts (cols) array for results
r = np.zeros((len(models),len(expts)))

import h5py
for mi in range(0,len(models)):
    for ei in range(0,len(expts)):
        filename = './log/agent/'+models[mi]+'_'+expts[ei]+'.h5'
        with h5py.File (filename, 'r') as f:
            r[mi][ei] = list(f['/rms'])[0]

print ('r: {0}'.format(r))

# Normalise each col of r to show best model for each condition?
r_n = r / r.max(axis=0)

ax = fig.add_subplot (1,1,1);
im = ax.imshow (r, cmap='jet')
ax.set_xticks(range(0, len(expts)))
ax.set_yticks(range(0, len(models)))
ax.set_xticklabels (expts)
ax.set_yticklabels (models)
pl.xticks(rotation=90)
pl.colorbar(im)
ax.set_title('RMS error. Low is good. Raw.')
for mi in range(0,len(models)):
    for ei in range(0,len(expts)):
        ax.text (ei-0.25,mi,'{0:.3f}'.format(r[mi][ei]), fontsize=12)

rs = np.sum(r, axis=1)
print ('rs: {0}'.format (rs))

#pl.tight_layout()
#pl.savefig('somefile.svg', transparent=True)
pl.show()
