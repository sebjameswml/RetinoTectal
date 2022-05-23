import numpy as np
from ephcolour import ephcol as C

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
fig = pl.figure(figsize=(10,10))

# Graphing code goes here.
import h5py
filename = 'Schmitt_Wnt-Ryk_2006_Fig2a.h5'

# File 1 has panel a - chick Wnt
with h5py.File (filename, 'r') as f:
    nf = list(f['/nframes'])[0]
    # Frame001 has the curve from panel B
    frame = 'Frame{0:03}'.format(1)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    ephrinB1_means = np.array(f[key])

    ax1 = fig.add_subplot (1,1,1)
    ax1.plot (ephrinB1_means, color=C.Wnt, label='Chick Ryk', linewidth=3)

    ax1.set_xticks([0,200])
    ax1.set_xticklabels(['Dors.','Vent.'])
    ax1.set_xlabel ('Retinal D-V axis')
    ax1.set_ylabel ('Normalised expression')
    ax1.legend()

pl.savefig('Schmitt_Wnt-Ryk_2006_Fig2a_Ryk.svg', transparent=True)
pl.show()
