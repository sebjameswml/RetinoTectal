import numpy as np
import matplotlib.pyplot as plt
from sebcolour import Colour as C
# Set plotting font defaults
import matplotlib
fs = 10
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

kiki_wt_top = genfromtxt(fname = 'brown_fig5c_kiki-wt-top.csv', delimiter=',', skip_header=1)
kiki_wt_bot = genfromtxt(fname = 'brown_fig5c_kiki-wt-bot.csv', delimiter=',', skip_header=1)

ki_kdkd_top = genfromtxt(fname = 'reber_fig4d_ki-kdkd-top.csv', delimiter=',', skip_header=0)
ki_kdkd_bot = genfromtxt(fname = 'reber_fig4d_ki-kdkd-bot.csv', delimiter=',', skip_header=0)

kiki_kd_top = genfromtxt(fname = 'reber_fig4e_kiki-kd-top.csv', delimiter=',', skip_header=0)
kiki_kd_bot = genfromtxt(fname = 'reber_fig4e_kiki-kd-bot.csv', delimiter=',', skip_header=0)

kiki_kdkd_top = genfromtxt(fname = 'reber_fig4f_kiki-kdkd-top.csv', delimiter=',', skip_header=0)
kiki_kdkd_bot = genfromtxt(fname = 'reber_fig4f_kiki-kdkd-bot.csv', delimiter=',', skip_header=0)

fig, ((ax1, ax2, ax3),(ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(16,6))

# 1:1 line for reference
ax1.plot ([1,0],[0,1], color=C.grey45, linestyle='--', dashes=(9, 3))
ax2.plot ([1,0],[0,1], color=C.grey45, linestyle='--', dashes=(9, 3))
ax3.plot ([1,0],[0,1], color=C.grey45, linestyle='--', dashes=(9, 3))
ax4.plot ([1,0],[0,1], color=C.grey45, linestyle='--', dashes=(9, 3))
ax5.plot ([1,0],[0,1], color=C.grey45, linestyle='--', dashes=(9, 3))
ax6.plot ([1,0],[0,1], color=C.grey45, linestyle='--', dashes=(9, 3))

ax1.plot(red[:,0],red[:,1],'o-',color=C.crimson)
ax1.plot(red1[:,0],red1[:,1],'s-',color=C.crimson)
ax1.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
ax1.set_ylabel('R {0} tectum {0} C'.format(u"\u27f6"))
ax1.set_xlim(0,1)
ax1.set_ylim(0,1)
ax1.set_aspect('equal', 'box')

ax2.plot(grey[:,0],grey[:,1],'o-',color=C.grey30)
ax2.plot(grey1[:,0],grey1[:,1],'s-',color=C.grey30)
ax2.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
ax2.set_ylabel('R {0} tectum {0} C'.format(u"\u27f6"))
ax2.set_xlim(0,1)
ax2.set_ylim(0,1)
ax2.set_aspect('equal', 'box')

ax3.plot(kiki_wt_top[:,0],kiki_wt_top[:,1],'o-',color=C.black)
ax3.plot(kiki_wt_bot[:,0],kiki_wt_bot[:,1],'s-',color=C.black)
ax3.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
ax3.set_ylabel('R {0} tectum {0} C'.format(u"\u27f6"))
ax3.set_xlim(0,1)
ax3.set_ylim(0,1)
ax3.set_aspect('equal', 'box')

ax4.plot(ki_kdkd_top[:,0],ki_kdkd_top[:,1],'o-',color=C.black)
ax4.plot(ki_kdkd_bot[:,0],ki_kdkd_bot[:,1],'s-',color=C.black)
ax4.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
ax4.set_ylabel('R {0} tectum {0} C'.format(u"\u27f6"))
ax4.set_xlim(0,1)
ax4.set_ylim(0,1)
ax4.set_aspect('equal', 'box')

ax5.plot(kiki_kd_top[:,0],kiki_kd_top[:,1],'o-',color=C.black)
ax5.plot(kiki_kd_bot[:,0],kiki_kd_bot[:,1],'s-',color=C.black)
ax5.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
ax5.set_ylabel('R {0} tectum {0} C'.format(u"\u27f6"))
ax5.set_xlim(0,1)
ax5.set_ylim(0,1)
ax5.set_aspect('equal', 'box')

ax6.plot(kiki_kdkd_top[:,0],kiki_kdkd_top[:,1],'o-',color=C.black)
ax6.plot(kiki_kdkd_bot[:,0],kiki_kdkd_bot[:,1],'s-',color=C.black)
ax6.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
ax6.set_ylabel('R {0} tectum {0} C'.format(u"\u27f6"))
ax6.set_xlim(0,1)
ax6.set_ylim(0,1)
ax6.set_aspect('equal', 'box')

show_titles = 0
if show_titles:
    ax1.set_title('ki-kd')
    ax2.set_title('ki-wt')
    ax3.set_title('kiki-wt')
    ax4.set_title('ki-kdkd')
    ax5.set_title('kiki-kd')
    ax6.set_title('kiki-kdkd')

plt.savefig('reber_fig3.svg')
plt.show()
