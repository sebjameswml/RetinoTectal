# Open CSV of Gosse Fig 2e. and do some fitting.

import numpy as np

wt = np.genfromtxt ('Gosse_Fig2e_WT.csv', delimiter=',', skip_header=1)

print ('wt retpos: {0}'.format(wt[:,0]))
print ('wt tecpos: {0}'.format(wt[:,1]))

# lak is non-competition
lak = np.genfromtxt ('Gosse_Fig2e_lak.csv', delimiter=',', skip_header=1)

import pylab as pl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

F = pl.figure (figsize=(8,4))

bins = 20

# Set x and y for WT
y, x_ = np.histogram(wt[:,1], bins)
x = (x_[1:] + x_[:-1])/2
xos = np.mean(x) # x offset for fitting
x = x - xos

def Gauss(x, A, B):
    y = A*np.exp(-1*B*x**2)
    return y

parameters, covariance = curve_fit(Gauss, x, y)

fit_A = parameters[0]
fit_B = parameters[1]

fit_y = Gauss(x, fit_A, fit_B)

sp = F.add_subplot (1, 2, 1)

sp.plot(x+xos, y, 's', label='WT')
sp.plot(x+xos, fit_y, '-', label='WT fit')
sp.set_xlabel('Retinal posn (percent N-T)')
sp.set_ylabel('Axon density')

# Set x and y for lak (non-competition condition)
y, x_ = np.histogram(lak[:,1], bins)
x = (x_[1:] + x_[:-1])/2
xos = np.mean(x) # x offset for fitting
x = x - xos

parameters, covariance = curve_fit(Gauss, x, y)

fit_A = parameters[0]
fit_B = parameters[1]

fit_y = Gauss(x, fit_A, fit_B)

sp.plot(x+xos, y, 'd', label='lak')
sp.plot(x+xos, fit_y, '-', label='lak fit')
sp.legend()

sp2 = F.add_subplot (1, 2, 2)
sp2.plot (wt[:,0], wt[:,1], 's', label='WT')
sp2.plot (lak[:,0], lak[:,1], 'd', label='lak')
sp2.set_xlabel('Retinal posn (percent N-T)')
sp2.set_ylabel('Tectal posn (percent A-P)')
sp2.legend()

pl.show()
