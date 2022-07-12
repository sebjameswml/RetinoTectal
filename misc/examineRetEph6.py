#
# clustersize = 1/r_A4
#

import numpy as np
import matplotlib.pyplot as plt
from sebcolour import Colour as C
clr_knockin = C.darkorchid1
clr_knockdown = C.cyan2
clr_wt = C.yellowgreen
# Set plotting font defaults
import matplotlib
fs = 8
fnt = {'family' : 'DejaVu Sans',
       'weight' : 'regular',
       'size'   : fs}
matplotlib.rc('font', **fnt)
# Important for svg output of text as 'things that can be edited in inkscape'
import pylab
pylab.rcParams['svg.fonttype'] = 'none'

# Simple 1/EphA4 relationship for clustersize
def clustersz (_EphAx, ki, _EphA4):
    cs = np.ones(len(_EphAx))/_EphA4
    return cs

##
## Epha3/Epha4
##
def examineRetEph():

    x = np.linspace(0,1,21)
    epha4_constant = 3.5
    # params
    axpow = 1
    axmult = 0.01

    # Knockin and knockdown, which should be copied from e_eph_ki-wt.json and e_eph_ki-kd.json (
    kd = 1
    ki = 1.0
    kiki = 3.0

    ## Binding affinity for EphAx. prop. to 1/K_D. see Monschau et al
    ## 1997 for dissociation constants K_D = 6.16e-10 for ephrinA5/EphA5
    ## and 1.44e-10 for ephrinA5/EphA3. Call it 3.8e-10 for EphAx,
    ## ignore the magnitude and let w_EphAx = 0.25.
    w_EphAx = 0.25
    ## Binding affinity parameter for EphA4. K_D = 6.22e-10 for
    ## ephrinA5, so binds slightly less well than EphA3/5. 0.611 = 3.8/6.22.
    w_EphA4 = w_EphAx * 0.611

    _exp = 0.26 * np.exp(2.3*x) + 1.05 # T exponential_expression (const T& x)

    ## Start with a certain amount of ephrinA (ephrinA5 in retina - Frisen et al)
    ephrinA = np.flip(_exp)

    ## Define an even expression of EphA4 (No nasal-temporal gradient)
    _epha4 = np.ones(len(ephrinA)) * epha4_constant
    _epha4_kd = _epha4 - kd

    ## Let EphA4 interact with cis ephrinA (Hornberger et al 1999)
    _p_epha4 = (w_EphA4 * _epha4 * ephrinA)

    ## Remaining epha4 could interact
    EphA4_free = _epha4 - _p_epha4
    EphA4_free_kd = EphA4_free - kd

    ## And some expression of EphA3/x whatever
    EphAx = _exp

    fs1 = (13,2.4)
    fs2 = (15,4)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=fs2)

    # WT
    ax1.plot (x, EphAx, linestyle='-', color=clr_wt, label='$r_0$ (EphA wildtype cells)')
    ax1.plot (x, ephrinA, linestyle=':', color=clr_wt, label='$l_0$ (ephrinA)')
    # EphA3 knockin
    ax1.plot (x, EphAx+ki, linestyle='-', color=clr_knockin, label='$r_0 + ki$ (EphA3 knock-in)')

    ax1.legend()
    ax1.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax1.set_ylabel('Expression')
    ax1.set_ylim(0,5)

    ax2.plot (x, _epha4, linestyle=':', color=clr_wt, label='EphA4')
    ax2.plot (x, EphA4_free, linestyle='-', color=clr_wt, label='$r_{A4}$  (free EphA4, wildtype)')
    ax2.plot (x, EphA4_free_kd, linestyle='-', color=clr_knockdown, label='$r_{A4} - kd$ (EphA4 knock-down)')
    ax2.legend()
    ax2.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax2.set_ylabel('Expression')
    ax2.set_ylim(0,5)

    # Cluster size vs. EphAx expression
    #ax31 = ax3.twiny()
    ax3.plot (EphA4_free, clustersz (EphAx, 0, EphA4_free), label='$c_0 = 1/r_{A4}$', color=clr_wt)
    # Or vs. EphAx/N/T posn
    #ax31.plot (x, clustersz (EphAx, 0, EphA4_free), label='$c_0 = 1/r_{A4}$', color=clr_wt)
    ax3.legend()
    ax3.set_xlabel('EphA4')
    #ax31.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax3.set_ylabel('EphA clustersize ($c_0$)')
    yl = ax3.get_ylim()
    yl = (0, yl[1])
    ax3.set_ylim(yl)

    e4power = 1
    ax4.plot (x, (EphAx) * (clustersz(EphAx, 0, EphA4_free)), linestyle='-', color=clr_wt)
    ax4.plot (x, (EphAx+ki) * (clustersz(EphAx, ki, EphA4_free)), linestyle='-', color=clr_knockin)
    ax4.plot (x, (EphAx) * (clustersz(EphAx, 0, EphA4_free_kd)), linestyle='-', color=clr_knockdown)
    ax4.plot (x, (EphAx+ki) * (clustersz(EphAx, ki, EphA4_free_kd)), linestyle='-', color=clr_knockdown)
    ax4.plot (x, (EphAx+ki) * (clustersz(EphAx, ki, EphA4_free_kd)), linestyle='--', color=clr_knockin, dashes=(5, 5))
    # Have to create a curated legend for the above plots
    filled_line1 = plt.Line2D([], [], linestyle="-", color=clr_wt)
    filled_line2 = plt.Line2D([], [], linestyle="-", color=clr_knockin)
    filled_line3 = plt.Line2D([], [], linestyle="-", color=clr_knockdown)
    dotted_line1 = plt.Line2D([], [], linestyle="-", color=clr_knockdown)
    dotted_line2 = plt.Line2D([], [], linestyle="--", color=clr_knockin, dashes=(5, 5))
    ax4.legend([filled_line1, filled_line2, filled_line3, (dotted_line1, dotted_line2)], ['$r_0/r_{A4}$ (wildtype)','$(r_0 + ki)/r_{A4}$ (EphA3 ki)', '${r_0}/(r_{A4}-kd)$ (EphA4 kd)', '$(r_0+ki)/(r_{A4}-kd)$ (ki + kd)'])

    ax4.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax4.set_ylabel('Signal')
    yl = ax4.get_ylim()
    yl = (0, yl[1])
    ax4.set_ylim(yl)

    plt.tight_layout()
    fn = 'epha_fig_div.svg'
    plt.savefig(fn)
    plt.show()

examineRetEph()
