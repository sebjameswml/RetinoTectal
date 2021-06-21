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
fig = pl.figure(figsize=(8,10))

# Graphing code goes here.
import h5py
filename = 'feldheim_genetic_2000_Fig4.h5'

with h5py.File (filename, 'r') as f:
    print("Keys: {0}".format(list(f.keys())))

    nf = list(f['/nframes'])[0]
    print ('Number of frames: {0}'.format(nf))

    # Frame001 has the curve from panel A, ephrin A2.
    frame = 'Frame{0:03}'.format(1)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    ephrinA2_means = np.array(f[key])

    # Frame002 has the curve from panel B, ephrin A5.
    frame = 'Frame{0:03}'.format(2)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    ephrinA5_means = np.array(f[key])

    ax1 = fig.add_subplot (1,1,1)

    ax1.plot (ephrinA2_means, color=C.ephrinA2, label='ephrinA2', linewidth=3)
    ax1.plot (ephrinA5_means, color=C.ephrinA5, label='ephrinA5', linewidth=3)
    #ax1.plot ((ephrinA5_means+ephrinA2_means)/2.0, color=C.mediumpurple2, label='(A5+A2)/2', linewidth=2, linestyle='dotted')

    ax1.set_xticks([0,100])
    ax1.set_xticklabels(['Ant.','Post.'])
    ax1.set_xlabel ('Collicular A-P axis')
    ax1.set_ylabel ('Normalised expression')
    ax1.legend()

    pl.savefig('feldheim_genetic_2000_Fig4.svg', transparent=True)
    pl.show()
