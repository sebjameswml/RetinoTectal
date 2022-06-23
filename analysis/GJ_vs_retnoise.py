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


# Will need to extract data from .h5 files
files=['../log/agent/ee_GJ_best_1_wt_figcomp2_exit_true_ret_rcpt_noise_gain_0.0_steps_1000.h5',
'../log/agent/ee_GJ_best_1_wt_figcomp2_exit_true_ret_rcpt_noise_gain_0.001_steps_1000.h5',
'../log/agent/ee_GJ_best_1_wt_figcomp2_exit_true_ret_rcpt_noise_gain_0.002_steps_1000.h5',
'../log/agent/ee_GJ_best_1_wt_figcomp2_exit_true_ret_rcpt_noise_gain_0.005_steps_1000.h5',
'../log/agent/ee_GJ_best_1_wt_figcomp2_exit_true_ret_rcpt_noise_gain_0.01_steps_1000.h5',
'../log/agent/ee_GJ_best_1_wt_figcomp2_exit_true_ret_rcpt_noise_gain_0.02_steps_1000.h5',
'../log/agent/ee_GJ_best_1_wt_figcomp2_exit_true_ret_rcpt_noise_gain_0.05_steps_1000.h5',
'../log/agent/ee_GJ_best_1_wt_figcomp2_exit_true_ret_rcpt_noise_gain_0.1_steps_1000.h5']

noise_gain = [0,.001,.002,.005,.01,.02,.05,.1]

#files=['../log/agent/ee_GJ_best_1_wt_figcomp2_exit_true_noise_gain_0.0_steps_1000.h5']

mean_eta = []
mean_eps = []
import h5py
for fn in files:
    with h5py.File (fn, 'r') as f:
        t = np.array(list(f['/t']))
        crosscount = np.array(list(f['/crosscount']))

        cc_idx = np.nonzero(np.where (crosscount > -1, 1, 0)) # incantation to make an index array
        cc_red = crosscount[cc_idx]
        t_red = t[cc_idx]
        rms = list(f['/rms'])
        # 500-1000 time slice is [41:] 750-1000 slice is [66:]
        mean_eta.append(np.mean(cc_red[66:]))
        # 500-1000 time slice is [100:] 750-1000 slice is [150:]
        mean_eps.append(np.mean(rms[150:]))

fig, ax1 = plt.subplots()

color = clr2
ax1.set_xlabel(r'$\nu$')
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

#ax2.set_ylim([0, 0.1])
#ax1.set_ylim([0,400])

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('GJ_vs_retnoise.svg', transparent=True)
plt.show()