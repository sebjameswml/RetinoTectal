#
# *Choose* expression functions so that model reproduces data.
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
    kd = 1.2
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

    pAx = 1.0
    pA4 = 0.0

    ## EphA4_expression_function_one
    EphA4_free = 3.5 * (1 - w_EphA4 * (0.26 * np.exp (1.6*(1-x)) + 2.35));

    ## EphA4_knockdown_function w_EphA4 = 0.15275
    EphA4_free_kd = -0.35 + 3.5 * (1 -  w_EphA4 * (0.26 * np.exp (1.0*(1-x)) + 2.35))



    ## And some expression of EphA3/x whatever
    EphAx = _exp

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8,2.4))

    # WT
    ax1.plot (x, EphAx, linestyle='-', color=clr_wt, label='EphAx')
    ax1.plot (x, ephrinA, linestyle=':', color=clr_wt, label='ephrinA')
    # EphA3 knockin
    ax1.plot (x, EphAx+ki, linestyle='-', color=clr_knockin, label='EphAx (knock-in)')

    ax1.legend()
    ax1.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax1.set_ylabel('Expression/Interaction')
    ax1.set_ylim(0,5)

    ax2.plot (x, _epha4, linestyle=':', color=clr_wt, label='EphA4')
    ax2.plot (x, EphA4_free, linestyle='-', color=clr_wt, label='EphA4 free')
    ax2.plot (x, EphA4_free_kd, linestyle='-', color=clr_knockdown, label='EphA4 free (kd)')
    ax2.legend()
    ax2.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax2.set_ylabel('Expression/Interaction')
    ax2.set_ylim(0,5)

    #ax3.plot (x, 1/EphA4_free, linestyle='-.', color=clr_wt, label='Cluster size (1/EphA4 free)')
    #ax3.plot (x, 1/(EphA4_free_kd), linestyle='-.', color=clr_knockdown, label='Cluster size (kd)')
    e4power = 1
    ax3.plot (x, (pAx*EphAx + pA4*EphA4_free)/(np.power(EphA4_free, e4power)), linestyle='-', color=clr_wt, label='Signal (WT)')
    ax3.plot (x, (pAx*EphAx + pA4*EphA4_free_kd)/(np.power(EphA4_free_kd, e4power)), linestyle='-', color=clr_knockdown, label='Signal (kd)')
    ax3.plot (x, (pAx*(EphAx+ki)+pA4*EphA4_free)/(np.power(EphA4_free, e4power)), linestyle='-', color=clr_knockin, label='Signal (ki)')
    ##ax3.plot (x, (EphAx+kiki)/(EphA4_free), linestyle='-', color=C.crimson, label='Signal (kiki)')
    ## How to make this alternating red/green?
    ax3.plot (x, (pAx*(EphAx+ki)+pA4*EphA4_free_kd)/(np.power(EphA4_free_kd, e4power)), linestyle='-', color=clr_knockdown)
    ax3.plot (x, (pAx*(EphAx+ki)+pA4*EphA4_free_kd)/(np.power(EphA4_free_kd, e4power)), linestyle='--', color=clr_knockin, dashes=(5, 5))
    ax3.legend()
    ax3.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax3.set_ylabel('Signal')
    #ax3.set_ylim(0,5)

    plt.tight_layout()
    fn = 'epha_fig.svg'
    plt.savefig(fn)
    plt.show()


examineRetEph()
