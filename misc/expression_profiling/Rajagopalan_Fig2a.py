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
fig = pl.figure(figsize=(16,10))

# Graphing code goes here.
import h5py
filename = 'Rajagopalan_Fig2a.h5'

# File contains two frames, with two separate curves for different expression (one deep, one superficial)
with h5py.File (filename, 'r') as f:
    nf = list(f['/nframes'])[0]
    # Frame001 has the all layers curve
    frame = 'Frame{0:03}'.format(1)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    neogen_sup_means = np.array(f[key])
    # Frame002 has the deep expression only
    frame = 'Frame{0:03}'.format(2)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    neogen_deep_means = np.array(f[key])

    ax1 = fig.add_subplot (1,2,1)
    ax1.plot (neogen_sup_means, color=C.neogen, label='neogen all layers.', linewidth=3)

    ax2 = fig.add_subplot (1,2,2)
    ax2.plot (neogen_deep_means, color=C.neogen, label='neogen deep', linewidth=3)

    ax1.set_xticks([0,200])
    ax1.set_xticklabels(['Temp.','Nas.'])
    ax1.set_xlabel ('Retinal T-N axis')
    ax1.set_ylabel ('Normalised expression')
    ax1.legend()

    ax2.set_xticks([0,200])
    ax2.set_xticklabels(['Temp.','Nas.'])
    ax2.set_xlabel ('Retinal T-N axis')
    ax2.set_ylabel ('Normalised expression')
    ax2.legend()

pl.savefig('Rajagopalan_Fig2a.svg', transparent=True)
pl.show()
