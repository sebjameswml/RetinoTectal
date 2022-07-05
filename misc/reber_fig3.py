import numpy as np
import matplotlib.pyplot as plt
from sebcolour import Colour as C
# Set plotting font defaults
import matplotlib
fs = 16
fnt = {'family' : 'DejaVu Sans',
       'weight' : 'regular',
       'size'   : fs}
matplotlib.rc('font', **fnt)
# Important for svg output of text as 'things that can be edited in inkscape'
import pylab
pylab.rcParams['svg.fonttype'] = 'none'

from numpy import genfromtxt
red = genfromtxt(fname = 'reber_fig3_red.csv', delimiter=',', skip_header=1)
grey = genfromtxt(fname = 'reber_fig3_grey.csv', delimiter=',', skip_header=1)
red = red / 100
grey = grey / 100
red1 = genfromtxt(fname = 'reber_fig3_red1.csv', delimiter=',', skip_header=1)
grey1 = genfromtxt(fname = 'reber_fig3_grey1.csv', delimiter=',', skip_header=1)
red1 = red1 / 100
grey1 = grey1 / 100

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16,6))

ax1.plot(red[:,0],red[:,1],'o-',color=C.crimson)
ax1.plot(red1[:,0],red1[:,1],'s-',color=C.crimson)
ax1.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
ax1.set_ylabel('R {0} tectum {0} C'.format(u"\u27f6"))
ax1.set_xlim(0,1)
ax1.set_ylim(0,1)


ax2.plot(grey[:,0],grey[:,1],'o-',color=C.grey30)
ax2.plot(grey1[:,0],grey1[:,1],'s-',color=C.grey30)
ax2.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
ax2.set_ylabel('R {0} tectum {0} C'.format(u"\u27f6"))
ax2.set_xlim(0,1)
ax2.set_ylim(0,1)

plt.savefig('reber_fig3.svg')
plt.show()
