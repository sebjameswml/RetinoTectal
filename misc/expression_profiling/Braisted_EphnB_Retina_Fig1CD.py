import numpy as np
from ephcolour import ephcol as C

# Set plotting font defaults
import matplotlib
fs = 28
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
filename = 'Braisted_EphnB_Retina_Fig1CD.h5'

with h5py.File (filename, 'r') as f:
    print("Keys: {0}".format(list(f.keys())))

    nf = list(f['/nframes'])[0]
    print ('Number of frames: {0}'.format(nf))

    # Frame001 has the curve from panel C, Ephrin B1, Retina
    frame = 'Frame{0:03}'.format(1)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    ephrinB1_means = np.array(f[key])

    # Frame002 has the curve from panel D, Ephrin B2, Retina
    frame = 'Frame{0:03}'.format(2)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    ephrinB2_means = np.array(f[key])

    ax1 = fig.add_subplot (1,1,1)

    ax1.plot (ephrinB1_means, color=C.ephrinB1, label='ephrinB1', linewidth=3)
    ax1.plot (ephrinB2_means, color=C.ephrinB2, label='ephrinB2', linewidth=3)

    ax1.set_xticks([0,200])
    ax1.set_xticklabels(['Dors.','Vent.'])
    ax1.set_xlabel ('Retinal D-V axis')
    ax1.set_ylabel ('Normalised expression')
    ax1.legend()

    pl.tight_layout()
    pl.savefig('Braisted_EphnB_Retina_Fig1CD.svg', transparent=True)
    pl.show()
