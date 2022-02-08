import numpy as np
import matplotlib.pyplot as plt
import csv
from sebcolour import Colour as C
clr1 = C.royalblue
clr2 = C.crimson
blk = C.black
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

# Important for svg output of text as 'things that can be edited in inkscape'
import pylab as pl
pl.rcParams['svg.fonttype'] = 'none'


# Will need to extract data from .h5 files:
'''
files=['../log/agent/eE_G_wt_figcomp2_exit_true_noise_gain_0.0_steps_1000.h5',
'../log/agent/eE_G_wt_figcomp2_exit_true_noise_gain_0.2_steps_1000.h5',
'../log/agent/eE_G_wt_figcomp2_exit_true_noise_gain_0.4_steps_1000.h5',
'../log/agent/eE_G_wt_figcomp2_exit_true_noise_gain_0.6_steps_1000.h5',
'../log/agent/eE_G_wt_figcomp2_exit_true_noise_gain_0.8_steps_1000.h5',
'../log/agent/eE_G_wt_figcomp2_exit_true_noise_gain_1.0_steps_1000.h5',
'../log/agent/eE_G_wt_figcomp2_exit_true_noise_gain_1.5_steps_1000.h5',
'../log/agent/eE_G_wt_figcomp2_exit_true_noise_gain_2.0_steps_1000.h5',]
'''

files=['../log/agent/eE_G_wt_figcomp2_exit_true_noise_gain_0.0_steps_1000.h5']

# Containers for means
import h5py
for fn in files:
    with h5py.File (fn, 'r') as f:
        print("Keys: {0}".format(list(f.keys())))
        # It's a bit ugly, the way you get the value of the nframes
        # object; you have to convert the hdf5 object into a list, then
        # get the 0th value of the list (which always has one element)
        crosscount = list(f['/crosscount'])
        rms = list(f['/crosscount'])
        t = list(f['/t'])
        # Matlab: c = crosscount(crosscount~=-1)
        print ('crosscount: {0}'.format(crosscount))
        cc_where = np.where(crosscount != -1)
        print ('t where {0}'.format(cc_where))


'''
with open ('GI_vs_noise.txt') as csvfile:
    rdr = csv.reader (csvfile, delimiter=',')
    first = True
    noise_gain = []
    mean_eta = []
    mean_eps = []
    for row in rdr:
        if first:
            # Don't process header
            first = False
        else:
            print ('noise_gain: {0} mean_eta {1}, mean_eps {2}'.format (row[0], row[1], row[2]))
            noise_gain.append(float(row[0]))
            mean_eta.append(float(row[1]))
            mean_eps.append(float(row[2]))

    fig, ax1 = plt.subplots()

    color = clr2
    ax1.set_xlabel('Noise gain')
    # The bar reproduces well in matplotlib, but not in svg/inkscape
    #ax1.set_ylabel(r'$\bar{\eta}$', rotation=0,  color=blk)
    ax1.set_ylabel(r'$\eta$', rotation=0,  color=blk)
    ax1.plot(noise_gain, mean_eta, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = clr1
    ax2.set_ylabel(r'$\epsilon$', rotation=0, color=blk)  # we already handled the x-label with ax1
    ax2.plot(noise_gain, mean_eps, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig('GI_vs_noise.svg', transparent=True)
    plt.show()
'''
