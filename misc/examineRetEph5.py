#
# Got idea for the functional form I need - need clustersize (EphAx, EphA4)
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

def clustersz (_EphAx, ki, _EphA4):
    #print ('_EphAx min: {0}'.format (np.min(_EphAx))) # It's 1.31 here.
    cs = 0.2 + 10 * np.exp(0.5*((_EphAx+ki)-np.min(_EphAx)))/(100 * np.power(_EphA4, 3))
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

    ## EphA4_expression_function_one
    ##EphA4_free = 3.5 * (1 - w_EphA4 * (0.26 * np.exp (1.6*(1-x)) + 2.35));

    ## EphA4_knockdown_function w_EphA4 = 0.15275
    ##EphA4_free_kd = -0.35 + 3.5 * (1 -  w_EphA4 * (0.26 * np.exp (1.0*(1-x)) + 2.35))



    ## And some expression of EphA3/x whatever
    EphAx = _exp

    fs1 = (13,2.4)
    fs2 = (18,4)
    fig, (ax1, ax2, ax3, ax31, ax4) = plt.subplots(1, 5, figsize=fs2)

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

    # Cluster size vs. EphAx/N/T posn
    # Or vs. EphAx:
    ax3.plot (EphAx, clustersz (EphAx, 0, 0.5), label='EphA4 = 0.5')
    ax3.plot (EphAx, clustersz (EphAx, 0, 1.5), label='EphA4 = 1.5')
    ax3.legend()
    ax3.set_xlabel('EphAx')
    ax3.set_ylabel('Clustersize')

    ax31.plot (x, clustersz (EphAx, 0, 0.5), label='EphA4 = 0.5')
    ax31.plot (x, clustersz (EphAx, 0, 1.5), label='EphA4 = 1.5')
    print ('EphAx: {0}'.format(EphAx))
    print ('EphAx ki: {0}'.format(EphAx+ki))
    ax31.plot (x, clustersz (EphAx, ki, 0.5), label='ki EphA3, EphA4 = 0.5')
    ax31.plot (x, clustersz (EphAx, ki, 1.5), label='ki EphA3, EphA4 = 1.5')
    ax31.legend()
    ax31.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax31.set_ylabel('Clustersize')

    e4power = 1
    ax4.plot (x, (EphAx) * (clustersz(EphAx, 0, EphA4_free)), linestyle='-', color=clr_wt, label='Signal (WT)')
    ax4.plot (x, (EphAx+ki) * (clustersz(EphAx, ki, EphA4_free)), linestyle='-', color=clr_knockin, label='Signal (ki)')
    ax4.plot (x, (EphAx) * (clustersz(EphAx, 0, EphA4_free_kd)), linestyle='-', color=clr_knockdown, label='Signal (kd)')
    ax4.plot (x, (EphAx+ki) * (clustersz(EphAx, ki, EphA4_free_kd)), linestyle='-', color=clr_knockdown)
    ax4.plot (x, (EphAx+ki) * (clustersz(EphAx, ki, EphA4_free_kd)), linestyle='--', color=clr_knockin, dashes=(5, 5))
    ax4.legend()
    ax4.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax4.set_ylabel('Signal')

    plt.tight_layout()
    fn = 'epha_fig.svg'
    plt.savefig(fn)
    plt.show()


examineRetEph()
