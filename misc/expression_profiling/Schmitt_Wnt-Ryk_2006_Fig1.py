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
filename = 'Schmitt_Wnt-Ryk_2006_Fig1a.h5'
filename2 = 'Schmitt_Wnt-Ryk_2006_Fig1bc.h5'

# File 1 has panel a - chick Wnt
with h5py.File (filename, 'r') as f:
    nf = list(f['/nframes'])[0]
    # Frame001 has the curve from panel B
    frame = 'Frame{0:03}'.format(1)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    ephrinB1_means = np.array(f[key])

    ax1 = fig.add_subplot (1,2,1)
    ax1.plot (ephrinB1_means, color=C.Wnt, label='Chick Wnt', linewidth=3)

    ax1.set_xticks([0,200])
    ax1.set_xticklabels(['Med.','Lat.'])
    ax1.set_xlabel ('Tectal M-L axis')
    ax1.set_ylabel ('Normalised expression')
    ax1.legend()

# File 2 has panel bc - Mouse Wnt and ephrinB1
with h5py.File (filename2, 'r') as f:
    nf = list(f['/nframes'])[0]
    # Frame001 has the curve from panel J
    frame = 'Frame{0:03}'.format(1)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    wnt_means = np.array(f[key])

    frame = 'Frame{0:03}'.format(2)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    ephrinB1_means = np.array(f[key])

    ax2 = fig.add_subplot (1,2,2)
    ax2.plot (wnt_means, color=C.Wnt, label='Mouse Wnt', linewidth=3)
    ax2.plot (ephrinB1_means, color=C.ephrinB1, label='Mouse ephrin-B1', linewidth=3)

    ax2.set_xticks([0,50])
    ax2.set_xticklabels(['Med.','Lat.'])
    ax2.set_xlabel ('Tectal M-L axis')
    ax2.set_ylabel ('Normalised expression')
    ax2.legend()

pl.savefig('Schmitt_Wnt-Ryk_2006.svg', transparent=True)
pl.show()
