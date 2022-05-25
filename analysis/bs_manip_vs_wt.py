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
    #print ('zdata: {0} and ydata: {1}'.format (zdata, ydata))

    # combine the data as if they were drawn from a single distribution
    x = np.concatenate ((zdata, ydata))
    xmean = np.mean(x)

    ymean = np.mean(ydata)
    zmean = np.mean(zdata)

    # Compute variances for the observed values:
    obsvarz = np.sum( np.power((zdata-zmean), 2) ) / (n-1)
    obsvary = np.sum( np.power((ydata-ymean), 2) ) / (m-1)

    # Compute the observed value of the studentised statistic (using separate
    # variances, rather than a pooled variance):
    tobs = (zmean - ymean) / np.sqrt(obsvary/m + obsvarz/n)
    #print("tobs={0} and xmean={1}".format(tobs,xmean))

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

    numbeyond = np.sum (txstar>=tobs)

    print ("{0} out of {1}".format(numbeyond, B))
    asl = numbeyond/B
    minasl = 1/B # Smallest possible achieved significance level in case asl==0

    return (asl, minasl)

def bootstrap_error_of_mean (data, B):
    resamples = resample_with_replacement (data, B)
    r_mean = np.mean(resamples, 1)
    std_err = np.sqrt(np.var(r_mean))
    return std_err

import h5py

files = ["../j4_ee_GIJ_eph_ki-wt_steps_1500_rcnt.h5"]
#files = ["../j4_ee_GJ_best_1_eph_ki-wt_exit_true_steps_1500_rcnt.h5",
#         "../j4_ee_GJ_best_1_eph_kiki-wt_exit_true_steps_1500_rcnt.h5",
#         "../j4_ee_GJ_best_1_eph_ki-kd_exit_true_steps_1500_rcnt.h5"]

for filename in files:
    print ('{0}'.format(filename))

    with h5py.File (filename, 'r') as f:
        rc   = np.array(list(f['/rc']))
        rc_m = np.array(list(f['/rc_m']))
        nt   = np.array(list(f['/nt']))
        nt_m = np.array(list(f['/nt_m']))

        plt.plot (nt, rc, marker='x', linestyle="none", color=clr1)
        plt.plot (nt_m, rc_m, marker='o', linestyle="none", color=clr2)

        # find unique values in nt
        #print ("nt: {0}".format(np.unique(nt)))
        #print ("nt_m: {0}".format(np.unique(nt_m)))

        nt_vals = np.unique(nt)

        nt_means_wt_mean = []
        nt_means_wt_stderr = []
        nt_means_wt_sd = []

        nt_means_manip_mean = []
        nt_means_manip_stderr = []
        nt_means_manip_sd = []

        # for each in nt_vals, find rc and rc_m's that match
        for ntv in nt_vals:
            idx_wt = np.asarray(nt==ntv).nonzero()
            idx_manip = np.asarray(nt_m==ntv).nonzero()

            rcs_wt = rc[idx_wt]
            rcs_manip = rc_m[idx_manip]
            # Now statistically test if rcs_wt and rcs_manip are from the same distribution or not
            asl, minasl = bootstrap_ttest_equalityofmeans (rcs_wt, rcs_manip, 2000)
            print ("for nt value {0}, asl of bootstrap test = {1}".format(ntv, (asl if asl > 0 else minasl)))

            # Now compute the mean and the std err of the mean for each
            nt_means_wt_mean.append (np.mean(rcs_wt))
            nt_means_manip_mean.append(np.mean(rcs_manip))

            nt_means_wt_sd.append (np.std(rcs_wt))
            nt_means_manip_sd.append(np.std(rcs_manip))

            rcs_wtstderr = bootstrap_error_of_mean (rcs_wt, 256)
            nt_means_wt_stderr.append (rcs_wtstderr)

            rcs_manipstderr = bootstrap_error_of_mean (rcs_manip, 256)
            nt_means_manip_stderr.append (rcs_manipstderr)

            #print ("wt mean: {0}({2}), manipulated mean: {1}({3})".format (np.mean(rcs_wt), np.mean(rcs_manip), rcs_wtstderr, rcs_manipstderr))

        # error bar plot
        plt.errorbar (nt_vals, nt_means_wt_mean, nt_means_wt_sd, capsize=10, color=clr1)
        plt.errorbar (nt_vals, nt_means_manip_mean, nt_means_manip_sd, capsize=10, color=clr2)

        plt.show()
