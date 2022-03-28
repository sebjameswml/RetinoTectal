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
def ttest_equalityofmeans (zdata_, ydata_, B):
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

ar = np.array([1.1,2.2,3.3,4.4])
br = np.array([2.2,2.2,3.3,6.6])
arres = resample_with_replacement (ar, 5)
print ('arres: {0}'.format (arres))
print ('arres means: {0}'.format (np.mean(arres, 1)))

ttest_equalityofmeans (ar, br, 10)

import h5py

files = ["../j4_ee_GJ_best_1_eph_ki-wt_exit_true_steps_1500_rcnt.h5",
         "../j4_ee_GJ_best_1_eph_kiki-wt_exit_true_steps_1500_rcnt.h5",
         "../j4_ee_GJ_best_1_eph_ki-kd_exit_true_steps_1500_rcnt.h5"]

for filename in files:
    print ('{0}'.format(filename))

    with h5py.File (filename, 'r') as f:
        rc   = np.array(list(f['/rc']))
        rc_m = np.array(list(f['/rc_m']))
        nt   = np.array(list(f['/nt']))
        nt_m = np.array(list(f['/nt_m']))

        #plt.plot (nt, rc)
        #plt.plot (nt_m, rc_m)

        # find unique values in nt
        #print ("nt: {0}".format(np.unique(nt)))
        #print ("nt_m: {0}".format(np.unique(nt_m)))

        nt_vals = np.unique(nt)
        # for each in nt_vals, find rc and rc_m's that match
        for ntv in nt_vals:
            idx1 = np.asarray(nt==ntv).nonzero()
            idx2 = np.asarray(nt_m==ntv).nonzero()
            #print ("for nt value {0}, idx in nt: {1} and in nt_m: {2}".format(ntv, idx1, idx2))
            #print ("for nt value {0}, idx in nt: {1} and in nt_m: {2}".format(ntv, idx1, idx2))

            #print ("for nt value {0}, rc values are {1} and rc_m values are {2}".format (ntv, rc[idx1],rc[idx2]))
            #plt.plot (nt[idx1], rc[idx1])
            #plt.plot (nt_m[idx2], rc_m[idx2])
            #plt.show()

            rcs1 = rc[idx1]
            rcs2 = rc_m[idx2]
            ## Now statistically test if rcs1 and rcs2 are from the same distribution or not
            asl, minasl = ttest_equalityofmeans (rcs1, rcs2, 2000)
            print ("for nt value {0}, asl of bootstrap test = {1}".format(ntv, (asl if asl > 0 else minasl)))
