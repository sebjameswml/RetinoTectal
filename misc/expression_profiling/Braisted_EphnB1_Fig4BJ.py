import numpy as np
from sebcolour import Colour as C

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

# Graphing code goes here.
import h5py
filename = 'Braisted_EphnB1_Fig4B.h5'
filename2 = 'Braisted_EphnB1_Fig4J.h5'


# File 1 has panel B; transverse E6 expression
with h5py.File (filename, 'r') as f:
    nf = list(f['/nframes'])[0]
    # Frame001 has the curve from panel B
    frame = 'Frame{0:03}'.format(1)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    ephrinB1_means = np.array(f[key])

    ax1 = fig.add_subplot (2,1,1)
    ax1.plot (ephrinB1_means, color=C.cornflowerblue, label='ephrinB1', linewidth=3)

    ax1.set_xticks([0,200])
    ax1.set_xticklabels(['Dors.','Vent.'])
    ax1.set_xlabel ('Tectal D-V axis')
    ax1.set_ylabel ('Normalised expression (E6)')
    ax1.legend()

# File 2 has panel J; sagg. E8 expression
with h5py.File (filename2, 'r') as f:
    nf = list(f['/nframes'])[0]
    # Frame001 has the curve from panel J
    frame = 'Frame{0:03}'.format(1)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    ephrinB1_means = np.array(f[key])

    ax2 = fig.add_subplot (2,1,2)
    ax2.plot (ephrinB1_means, color=C.cornflowerblue, label='ephrinB1', linewidth=3)

    ax2.set_xticks([0,200])
    ax2.set_xticklabels(['Ant.','Post.'])
    ax2.set_xlabel ('Tectal A-P axis')
    ax2.set_ylabel ('Normalised expression (E8)')
    ax2.legend()

pl.savefig('Braisted_EphnB1_Fig4BJ.svg', transparent=True)
pl.show()
