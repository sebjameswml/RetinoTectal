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

##
## Epha3/Epha4
##
def examineRetEph():

    x = np.linspace(0,1,21)
    epha4_constant = 3.5
    # params
    axpow = 1
    axmult = 0.01

    # Knockin and knockdown, whichshould be copied from e_eph_ki-wt.json and e_eph_ki-kd.json (
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

    ## Let EphA4 interact with cis ephrinA (Hornberger et al 1999)
    _p_epha4 = (w_EphA4 * _epha4 * ephrinA)

    ## Remaining epha4 could interact
    EphA4_free = _epha4 - _p_epha4

    ephrinA_free = ephrinA - _p_epha4

    ## And some expression of EphA3/x whatever
    EphAx = _exp

    ## Some ephrinA binds to EphAx. Really, need a first order model that takes into account amount of EphAx and amount of ephrinA
    ephrinA_binds_EphAx = w_EphAx * EphAx * ephrinA_free

    ## Remaining free receptors can interact with tectal ephrins
    #EphAx_free = EphAx - ephrinA_binds_EphAx # unused

    ephrinA_free2 = ephrinA_free - ephrinA_binds_EphAx

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16,6))

    # WT
    ax1.plot (x, EphAx, linestyle='-', color=C.royalblue, label='EphAx')
    ax1.plot (x, ephrinA, linestyle=':', color=C.royalblue, label='ephrinA')
    # EphA3 knockin
    ax1.plot (x, EphAx+ki, linestyle='-', color=C.red, label='EphAx (knock-in)')

    ax1.legend()
    ax1.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax1.set_ylabel('Expression/Interaction')
    ax1.set_ylim(0,5)

    ax2.plot (x, _epha4, linestyle=':', color=C.royalblue, label='EphA4')
    ax2.plot (x, EphA4_free, linestyle='-', color=C.royalblue, label='EphA4 free')
    ax2.plot (x, _p_epha4, linestyle='--', color=C.royalblue, label='EphA4 cis bound')
    ax2.plot (x, EphA4_free-kd, linestyle='-', color=C.springgreen3, label='EphA4 free (kd)')
    ax2.legend()
    ax2.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax2.set_ylabel('Expression/Interaction')
    ax2.set_ylim(0,5)

    ax3.plot (x, 1/EphA4_free, linestyle='-.', color=C.royalblue, label='Cluster size (1/EphA4 free)')
    ax3.plot (x, 1/(EphA4_free-kd), linestyle='-.', color=C.springgreen3, label='Cluster size (kd)')
    ax3.plot (x, EphAx/EphA4_free, linestyle='-', color=C.royalblue, label='Signal (WT)')
    ax3.plot (x, EphAx/(EphA4_free-kd), linestyle='-', color=C.springgreen3, label='Signal (kd)')
    ax3.plot (x, (EphAx+ki)/(EphA4_free), linestyle='-', color=C.red, label='Signal (ki)')
    ax3.plot (x, (EphAx+kiki)/(EphA4_free), linestyle='-', color=C.crimson, label='Signal (kiki)')
    ax3.plot (x, (EphAx+ki)/(EphA4_free-kd), linestyle='-', color=C.goldenrod1, label='Signal (ki+kd)')
    ax3.legend()
    ax3.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax3.set_ylabel('Cluster size')
    ax3.set_ylim(0,5)

    plt.tight_layout()
    fn = 'epha_fig.svg'
    plt.savefig(fn)
    plt.show()


examineRetEph()
