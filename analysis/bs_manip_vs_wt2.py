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
fnt = {'family' : 'DejaVu Sans',
       'weight' : 'regular',
       'size'   : fs}
matplotlib.rc('font', **fnt)
# Important for svg output of text as 'things that can be edited in inkscape'
import pylab
pylab.rcParams['svg.fonttype'] = 'none'


rng = np.random.default_rng(8947584) # number is seed
def resample_with_replacement (data, B):
    data_n = len(data)
    resampled = np.zeros ((B, data_n), dtype=data.dtype) # match data type
    for i in range(0,B):
        rints = rng.integers(low=0, high=data_n, size=data_n)
        #print ('data = {0}'.format (data))
        #print ('rints = {0}'.format (rints))
        resampled[i,:] = data[rints]
        #print ('resampled {0}: {1}'.format(i, resampled[i,:]))
    return resampled

# Compute a bootstrapped two sample t statistic as per algorithm 16.2
# in Efron & Tibshirani.
# zdata is treatment; ydata is control. B is the number of bootstrap samples to make
#
# This tests for equality of means, not whether the populations are identical
# and may be used on distributions where the variances may not be equal.
#
# Re-coded from R (See https://github.com/ABRG-Models/GPR-BSB/blob/master/labbook/bootstrap_functions.r)
def bootstrap_ttest_equalityofmeans (zdata_, ydata_, B):
    # Ensure that the group which we name zdata is the larger one.
    if (np.mean(zdata_) > np.mean(ydata_)):
        zdata = zdata_;
        ydata = ydata_;
    else:
        zdata = ydata_;
        ydata = zdata_;

    n = len(zdata)
    m = len(ydata)
    print ('zdata: {0}'.format (zdata))
    print ('ydata: {0}'.format (ydata))

    # combine the data as if they were drawn from a single distribution
    x = np.concatenate ((zdata, ydata))
    xmean = np.mean(x)

    ymean = np.mean(ydata)
    zmean = np.mean(zdata)
    print ('means: xmean = {0} ymean = {1} ({3} dpoints) zmean = {2} ({4} dpoints)'.format(xmean,ymean,zmean,n,m))

    # Compute variances for the observed values:
    obsvarz = np.sum( np.power((zdata-zmean), 2) ) / (n-1)
    obsvary = np.sum( np.power((ydata-ymean), 2) ) / (m-1)
    print ('obsvarz = {0} and obsvary = {1}'.format(obsvarz, obsvary))

    # Compute the observed value of the studentised statistic (using separate
    # variances, rather than a pooled variance):
    tobs = (zmean - ymean) / np.sqrt(obsvary/m + obsvarz/n)
    print("tobs={0} and xmean={1}".format(tobs,xmean))

    # Create shifted distributions; shifted by group mean and combined mean:
    ztilda = zdata - np.mean(zdata) + xmean
    ytilda = ydata - np.mean(ydata) + xmean

    # Resample from the shifted (tilda) distributions:
    zstar = resample_with_replacement (ztilda, B)
    ystar = resample_with_replacement (ytilda, B)

    # Create vectors of the means of these resamples:
    zstarmeans = np.mean(zstar, 1)
    ystarmeans = np.mean(ystar, 1)

    # Compute the variances
    #print ('zstar: {0} and zstarmeans: {1}'.format(zstar, np.tile(zstarmeans,(n,1)).T ))
    zvariances = np.sum( np.power((zstar-np.tile(zstarmeans,(n,1)).T), 2), 1 ) / (n-1)
    yvariances = np.sum( np.power((ystar-np.tile(ystarmeans,(m,1)).T), 2), 1 ) / (m-1)

    top = zstarmeans - ystarmeans
    bot = np.sqrt(yvariances/m + zvariances/n)
    txstar = top / bot
    print ('txstar: {0}'.format(txstar))

    numbeyond = np.sum (txstar>=tobs)

    print ("=================")
    print ("{0} out of {1}".format(numbeyond, B))
    print ("=================")
    asl = numbeyond/B
    minasl = 1/B # Smallest possible achieved significance level in case asl==0

    return (asl, minasl)

def bootstrap_error_of_mean (data, B):
    resamples = resample_with_replacement (data, B)
    r_mean = np.mean(resamples, 1)
    std_err = np.sqrt(np.var(r_mean))
    return std_err

import h5py

# These are RC-vs-NT graphs generated by the simulation. Run the three
# simulations for the Genetic manipulation sims to create these.

# 0: ORIGINAL model 1: Clustersize model 2: Clustersize + r2 collapse
model_type = 2
doitall = 1
if doitall:
    filebases = ["../rcnt/j4_ee_GJ_best_1_eph_ki-wt_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_eph_kiki-wt_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_eph_ki-kd_lite_exit_true_steps_1500",
                     "../rcnt/j4_ee_GJ_best_1_EphA4_eph_ki-wt_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_eph_kiki-wt_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_eph_ki-kd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_eph_ki-kdkd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_eph_kiki-kdkd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_eph_kiki-kd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_eph_wt-kd_exit_true_steps_1500",
                     "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_ki-wt_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_kiki-wt_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_ki-kd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_ki-kdkd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_kiki-kdkd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_kiki-kd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_wt-kd_exit_true_steps_1500"]
else:
    if model_type == 0:
        filebases = ["../rcnt/j4_ee_GJ_best_1_eph_ki-wt_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_eph_kiki-wt_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_eph_ki-kd_lite_exit_true_steps_1500"]
    elif model_type == 1:
        # The FIXED model, with clustersize modification
        filebases = ["../rcnt/j4_ee_GJ_best_1_EphA4_eph_ki-wt_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_eph_kiki-wt_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_eph_ki-kd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_eph_ki-kdkd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_eph_kiki-kdkd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_eph_kiki-kd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_eph_wt-kd_exit_true_steps_1500"]
    elif model_type == 2:
        filebases = ["../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_ki-wt_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_kiki-wt_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_ki-kd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_ki-kdkd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_kiki-kdkd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_kiki-kd_exit_true_steps_1500",
                         "../rcnt/j4_ee_GJ_best_1_EphA4_r2collapse_eph_wt-kd_exit_true_steps_1500"]

from numpy import genfromtxt
red = genfromtxt(fname = '../misc/reber_fig3_red.csv', delimiter=',', skip_header=1)
grey = genfromtxt(fname = '../misc/reber_fig3_grey.csv', delimiter=',', skip_header=1)
red = red / 100
grey = grey / 100

#from fnmatch import fnmatch
import os

for fb in filebases:
    print ('File base: {0} (type {1})'.format(fb, type(fb)))

    # For each filebase, collect stats and draw a graph
    rc   = []
    rc_m = []
    nt   = []
    nt_m = []

    # Now find all files that match
    for r, d, fns in os.walk('../rcnt/'):
        for filename in fns:
            #print ('Does filename: {0} match filebase {1}?'.format (os.path.join(r, filename), fb))
            if fb in os.path.join(r, filename):

                filepath = os.path.join(r, filename)
                if '.svg' in filepath:
                    continue
                print('Processing: {0}'.format(filepath))

                with h5py.File (filepath, 'r') as f:
                    rc.append(np.array(list(f['/rc'])))
                    rc_m.append(np.array(list(f['/rc_m'])))
                    nt.append(np.array(list(f['/nt'])))
                    nt_m.append(np.array(list(f['/nt_m'])))

    print ('Sizes of rc/rc_m/nt/nt_m: {0}/{1}/{2}/{3}'.format (len(rc),len(rc_m),len(nt),len(nt_m)))
    print ('Shape of rc[0]: {0}'.format (np.shape (rc[0]))) # rc is a list we can iterate through.

    #fig,axes = plt.subplots (1, len(nt))
    #ii = int(0)
    #for ax in axes:
    #    # Two lines plot the points
    #    ax.plot (nt[ii], rc[ii], marker='x', linestyle="none", color=clr1)
    #    ax.plot (nt_m[ii], rc_m[ii], marker='o', linestyle="none", color=clr2)
    #    ii = ii + int(1)

    # find unique values in nt
    #print ("nt: {0}".format(np.unique(nt)))
    #print ("nt_m: {0}".format(np.unique(nt_m)))

    nt_means_wt_mean = []
    nt_means_wt_stderr = []
    nt_means_wt_sd = []

    nt_means_manip_mean = []
    nt_means_manip_stderr = []
    nt_means_manip_sd = []

    all_rcs_wt = []
    all_rcs_manip = []

    nt_vals = np.unique(nt[0])

    for ntv in nt_vals:
        print ('processing nt value {0}'.format(ntv))

        rcs_wt = np.array([])
        rcs_manip = np.array([])

        for i in range(0,len(rc)):
            #print ('Size of nt_vals is (cf nt): {0}'.format(len(nt_vals)))
            #print ('nt_vals: {0}'.format(nt_vals))
            # for each in nt_vals, find rc and rc_m's that match

            #print ('ntv: {0}'.format(ntv))
            #print ('nt[i]: {0}'.format(nt[i]))
            idx_wt = np.asarray(nt[i]==ntv).nonzero()
            #print ('idx_wt: {0}'.format(idx_wt))
            idx_manip = np.asarray(nt_m[i]==ntv).nonzero()
            rci = rc[i]
            #print ('rci: {0}'.format (rci))
            #print ('rc[i][idx_wt]: {0}'.format (rc[i][idx_wt]))
            rcs_wt = np.append(rcs_wt, rc[i][idx_wt])
            rcs_manip = np.append(rcs_manip, rc_m[i][idx_manip])

        print ('rcs_wt size: {0} shape: {1}'.format (len(rcs_wt), np.shape(rcs_wt)))
        print ('rcs_manip size: {0}'.format (len(rcs_manip)))

        # Now statistically test if rcs_wt and rcs_manip are from the same distribution or not

        print ('rcs_wt: {0}'.format(rcs_wt))
        print ('rcs_manip: {0} NOW BOOTSTRAP'.format(rcs_manip))
        asl, minasl = bootstrap_ttest_equalityofmeans (rcs_wt, rcs_manip, 10000)
        print ("for nt value {0}, asl of bootstrap test = {1} (raw {2})".format(ntv, (asl if asl > 0 else minasl), asl))

        # Now compute the mean and the std err of the mean for each
        nt_means_wt_mean.append (np.mean(rcs_wt))
        nt_means_manip_mean.append(np.mean(rcs_manip))

        nt_means_wt_sd.append (np.std(rcs_wt))
        nt_means_manip_sd.append(np.std(rcs_manip))

        rcs_wtstderr = bootstrap_error_of_mean (rcs_wt, 256)
        nt_means_wt_stderr.append (rcs_wtstderr)

        rcs_manipstderr = bootstrap_error_of_mean (rcs_manip, 256)
        nt_means_manip_stderr.append (rcs_manipstderr)

        #print ('rcs_wt: {0} ntv: {1}'.format (rcs_wt, ntv))
        all_rcs_wt.append (rcs_wt)
        all_rcs_manip.append (rcs_manip)
        #print ("wt mean: {0}({2}), manipulated mean: {1}({3})".format (np.mean(rcs_wt), np.mean(rcs_manip), rcs_wtstderr, rcs_manipstderr))

    fig1,(ax) = plt.subplots (1, 1, figsize=(6,6))

    # 1:1 line for reference
    ax.plot ([1,0],[0,1], color=C.grey45, linestyle='--', dashes=(9, 3))

    # error bar plot
    #plt.errorbar (nt_vals, nt_means_wt_mean, nt_means_wt_sd, capsize=10, color=clr1)
    #plt.errorbar (nt_vals, nt_means_manip_mean, nt_means_manip_sd, capsize=10, color=clr2)

    # Violin plots
    vbits = ax.violinplot (all_rcs_wt, nt_vals, widths=0.05)
    for pc in vbits['bodies']:
        pc.set_facecolor(clr1)
        pc.set_edgecolor(clr1)
    for partname in ('cbars','cmins','cmaxes'):
        vp = vbits[partname]
        vp.set_edgecolor(clr1)
        vp.set_linewidth(1)

    vbits = ax.violinplot (all_rcs_manip, nt_vals, widths=0.05)
    for pc in vbits['bodies']:
        pc.set_facecolor(clr2)
        pc.set_edgecolor(clr2)
    for partname in ('cbars','cmins','cmaxes'):
        vp = vbits[partname]
        vp.set_edgecolor(clr2)
        vp.set_linewidth(1)

    ax.plot (nt_vals, nt_means_wt_mean, color=clr1)
    ax.plot (nt_vals, nt_means_manip_mean, color=clr2)

    # Expt data
    show_expt = 0
    if show_expt:
        if 'eph_ki-kd' in fb:
            ax.plot (red[:,0],red[:,1],'o',color=C.crimson)
        elif 'eph_ki-wt' in fb:
            ax.plot(grey[:,0],grey[:,1],'o',color=C.grey30)

    ax.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax.set_ylabel('R {0} tectum {0} C'.format(u"\u27f6"))

    ax.set_aspect('equal', 'box')

    plt.tight_layout()

    plt.savefig('{0}.svg'.format(fb))
    plt.show()
