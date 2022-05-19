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
filename = 'Monnier_Fig2_RGM_expression_cutout_grey.h5'

# File contains two frames, with two separate curves for different expression (one deep, one superficial)
with h5py.File (filename, 'r') as f:
    nf = list(f['/nframes'])[0]
    # Frame001 has the superficial curve
    frame = 'Frame{0:03}'.format(1)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    RGM_sup_means = np.array(f[key])
    # Frame002 has the deeper curve
    frame = 'Frame{0:03}'.format(2)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    RGM_deep_means = np.array(f[key])

    ax1 = fig.add_subplot (1,2,1)
    ax1.plot (RGM_sup_means, color=C.RGM, label='RGM sup.', linewidth=3)

    ax2 = fig.add_subplot (1,2,2)
    ax2.plot (RGM_deep_means, color=C.RGM, label='RGM deep', linewidth=3)

    ax1.set_xticks([0,200])
    ax1.set_xticklabels(['Post.','Ant.'])
    ax1.set_xlabel ('Tectal P-A axis')
    ax1.set_ylabel ('Normalised expression')
    ax1.legend()

    ax2.set_xticks([0,200])
    ax2.set_xticklabels(['Post.','Ant.'])
    ax2.set_xlabel ('Tectal P-A axis')
    ax2.set_ylabel ('Normalised expression')
    ax2.legend()

pl.savefig('Monnier_Fig2_RGM_expression.svg', transparent=True)
pl.show()
