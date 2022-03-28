# Load rcnt.h5 and then compare, for each N-T position, whether the manipulated
# R-C position is statistically significantly different from the wildtype R-C
# position.
import numpy as np
import matplotlib.pyplot as plt
import csv
from sebcolour import Colour as C
clr1 = C.royalblue
clr2 = C.crimson
# Set plotting font defaults
import matplotlib
fs = 18
fnt = {'family' : 'Arial',
       'weight' : 'regular',
       'size'   : fs}
matplotlib.rc('font', **fnt)
# Important for svg output of text as 'things that can be edited in inkscape'
import pylab
pylab.rcParams['svg.fonttype'] = 'none'

import h5py
with h5py.File ("../rcnt.h5", 'r') as f:
    rc   = np.array(list(f['/rc']))
    rc_m = np.array(list(f['/rc_m']))
    nt   = np.array(list(f['/nt']))
    nt_m = np.array(list(f['/nt_m']))

    plt.plot (nt, rc)
    plt.plot (nt_m, rc_m)

    # find unique values in nt
    print ("nt: {0}".format(np.unique(nt)))
    print ("nt_m: {0}".format(np.unique(nt_m)))

    nt_vals = np.unique(nt)
    # for each in nt_vals, find rc and rc_m's that match
    for ntv in nt_vals:
        idx1 = np.asarray(nt==ntv).nonzero()
        idx1 = np.asarray(nt_m==ntv).nonzero()
        print ("for ntv {0}, idx: {1}".format(ntv, idx1))

    plt.show()
