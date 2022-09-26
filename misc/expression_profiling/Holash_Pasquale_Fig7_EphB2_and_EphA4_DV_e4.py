import numpy as np
from ephcolour import ephcol as C
from numpy.polynomial import Polynomial as Poly

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
fig = pl.figure(figsize=(18,18))

# Graphing code goes here.
import h5py
filename = 'Holash_Pasquale_Fig7_EphB2_and_EphA4_DV_e4.h5'
filename2 = 'Holash_Pasquale_Fig7_EphB2_and_EphA4_DV_e6.h5'
filename3 = 'Holash_Pasquale_Fig7_EphB2_and_EphA4_DV_e8.h5'

# File contains two frames, with two separate curves for different expression (one deep, one superficial)
with h5py.File (filename, 'r') as f, h5py.File (filename2, 'r') as f2, h5py.File (filename3, 'r') as f3 :
    nf = list(f['/nframes'])[0]

    frame = 'Frame{0:03}'.format(1)
    key = '{0}/signal/postproc/boxes/means_autoscaled'.format(frame)
    # Note np.flip, as curve was drawn from ventral to dorsal:
    EphB2_means = np.flip(np.array(f[key]))
    EphB2_means_part = EphB2_means[60:-60]
    x_ = range(0,len(EphB2_means_part))

    EphB2_e6_means = np.array(f2[key])
    EphB2_e8_means = np.array(f3[key])


    p = Poly.fit(x_, EphB2_means_part, deg=1)
    plin = p.linspace(80) # obtains values suitable for plotting
    print ('EphB2_means_part size: {0} and plin[0] size: {1}'.format (len(EphB2_means_part), len(plin[0])))

    d = plin[1] - EphB2_means_part
    dsq = d * d
    drms = np.sqrt(dsq)

    ax1 = fig.add_subplot (2,2,1)
    ax1.plot (EphB2_means, color=C.EphB2, label='EphB2 (E4)', linewidth=3)

    ax2 = fig.add_subplot (2,2,2)
    ax2.plot (EphB2_means_part, color=C.EphB2, label='EphB2 portion (E4)', linewidth=3)
    ax2.plot (plin[0], plin[1], color=C.EphB2, label='linear fit', linewidth=3, linestyle=':')

    ax2.plot (plin[0], drms, color=C.EphB2, label='RMS err', linewidth=1, linestyle='-.')
    me = np.mean (drms)
    ax2.plot ([x_[0],x_[-1]], [me, me], color=C.EphB2, label='Mean RMS err', linewidth=1, linestyle='--')
    print ('The mean RMS error is {0}'.format(me))

    ax3 = fig.add_subplot (2,2,3)
    ax3.plot (EphB2_e6_means, color=C.EphB2, label='EphB2 (E6)', linewidth=3, linestyle='-')

    ax4 = fig.add_subplot (2,2,4)
    ax4.plot (EphB2_e8_means, color=C.EphB2, label='EphB2 (E8)', linewidth=3, linestyle='-')

    ax1.set_xticks([0,200])
    # Ant. to right
    ax1.set_xticklabels(['Dor.','Vent.'])
    ax1.set_xlabel ('Retinal axis')
    ax1.set_ylabel ('Normalised expression')
    ax1.legend()

    ax2.set_xticks([0,len(EphB2_means_part)])
    ax2.set_xticklabels(['Dor.','Vent.'])
    ax2.set_xlabel ('Retinal axis')
    ax2.set_ylabel ('Normalised expression')
    ax2.legend()

    ax3.set_xticks([0,200])
    # Ant. to right
    ax3.set_xticklabels(['Dor.','Vent.'])
    ax3.set_xlabel ('Retinal axis')
    ax3.set_ylabel ('Normalised expression')
    ax3.legend()

    ax4.set_xticks([0,200])
    # Ant. to right
    ax4.set_xticklabels(['Dor.','Vent.'])
    ax4.set_xlabel ('Retinal axis')
    ax4.set_ylabel ('Normalised expression')
    ax4.legend()

pl.savefig('Holash_Pasquale_Fig7_EphB2_and_EphA4_DV_e468.svg', transparent=True)
pl.show()
