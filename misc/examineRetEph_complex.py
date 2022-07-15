#
# This script explores clustersize as a function of both EphAx, and EphA4
#

import numpy as np
import matplotlib.pyplot as plt
from sebcolour import Colour as C
clr_knockin = C.darkorchid1
clr_knockdown = C.cyan2
clr_wt = C.yellowgreen
# Set plotting font defaults
import matplotlib
fs = 12
fnt = {'family' : 'DejaVu Sans',
       'weight' : 'regular',
       'size'   : fs}
matplotlib.rc('font', **fnt)
# Important for svg output of text as 'things that can be edited in inkscape'
import pylab
pylab.rcParams['svg.fonttype'] = 'none'


## Computation of a cluster size. Matches C++ fn branch::clustersz() in branch.h
def clustersz_complex (_EphAx, ki, _EphA4, _EphA4_cisbound):
    # Take cisbound and offset it by a static amount:
    #cisbound = (_EphA4_cisbound-0.8404305)
    cisbound = (_EphA4_cisbound-1.3)
    # or (and this is harder to justify):
    #cisbound = (_EphA4_cisbound-np.min(_EphA4_cisbound))

    # If cisbound>0, set to cisbound, otherwise set to 0
    ##cisbound = np.where (cisbound>0, cisbound, 0)
    # Numerator is an exponential function of EphAx expression.
    numer = cisbound * np.exp(0.5*((_EphAx+ki)-np.min(_EphAx)))
    # Denominator is some power of EphA4. The + 2 limits how large the clustersize can become.
    denom = (10 * np.power(_EphA4, 2) + 2)
    #x = np.exp(np.linspace(1,0,21))
    #denom = (denom * x)-15
    #denom = (10 * np.power(_EphA4, 3) + 1)
    ratio = numer / denom
    cs = 3 * (0.2 + ratio)
    return (cs, numer, denom)

def clustersz_testing (_EphAx, ki, _EphA4, _EphA4_cisbound):
    cisbound = (_EphA4_cisbound-0)
    # If cisbound>0, set to cisbound, otherwise set to 0
    cisbound = np.power (cisbound, 2)
    cs = np.ones (len(_EphAx)) * cisbound
    return cs

## Computation of a cluster size, simplest scheme
def clustersz_simple (_EphAx, ki, _EphA4, _EphA4_cisbound):
    cs = np.ones (len(_EphAx)) / np.power(_EphA4, 1.2)
    return cs

## The signal function.
def signal (rcpt, cs):
    return rcpt * cs

##
## Epha3/Epha4
##
def examineRetEph():

    x = np.linspace(0,1,21)
    epha4_constant = 3.5

    # Knockin and knockdown, which should be copied from e_eph_ki-wt.json and e_eph_ki-kd.json (
    ki = 0.9
    kd = 0.7
    kiki = 3.0

    ## Binding affinity for EphAx. prop. to 1/K_D. see Monschau et al
    ## 1997 for dissociation constants K_D = 6.16e-10 for ephrinA5/EphA5
    ## and 1.44e-10 for ephrinA5/EphA3. Call it 3.8e-10 for EphAx,
    ## ignore the magnitude and let w_EphAx = 0.25.
    w_EphAx = 0.3
    ## Binding affinity parameter for EphA4. K_D = 6.22e-10 for
    ## ephrinA5, so binds slightly less well than EphA3/5. 0.611 = 3.8/6.22.
    w_EphA4 = w_EphAx * 0.611

    _exp = 0.26 * np.exp(2.3*x) + 1.05 # T exponential_expression (const T& x)

    ## Start with a certain amount of ephrinA (ephrinA5 in retina - Frisen et al)
    ephrinA = np.flip(_exp)

    ## Choose between simple knockdown, where we just take the free EphA4 and
    ## subtract a scalar, vs knocking down the original amount of EPhA4 and
    ## re-multiplying by retinal ephrinA level
    simple_knockdown = 1

    if simple_knockdown == 0:
        # Have a larger knockdown if we're not doing the 'simple' knockdown thing
        kd = 2 * kd

    ## Define an even expression of EphA4 (No nasal-temporal gradient)
    _epha4 = np.ones(len(ephrinA)) * epha4_constant
    _epha4_kd = _epha4 - kd

    ## Phosphorylised EphA4
    _p_epha4    = (w_EphA4 * _epha4    * ephrinA)
    ## Remaining epha4 could interact
    EphA4_free = _epha4 - _p_epha4

    ## Let EphA4 interact with cis ephrinA (Hornberger et al 1999)
    if simple_knockdown:
        EphA4_free_kd = EphA4_free - kd
        _p_epha4_kd = _p_epha4 + kd
    else:
        _p_epha4_kd = (w_EphA4 * _epha4_kd * ephrinA)
        # Knockdown static level of EphA4
        EphA4_free_kd = _epha4_kd - _p_epha4_kd

    ## And some expression of EphA3/x whatever
    EphAx = _exp

    # Have to create a curated legend for some plots (those with dashed coloured line)
    filled_line_wt = plt.Line2D([], [], linestyle="-", color=clr_wt)
    filled_line_ki = plt.Line2D([], [], linestyle="-", color=clr_knockin)
    filled_line_kd = plt.Line2D([], [], linestyle="-", color=clr_knockdown)
    dotted_line1 = plt.Line2D([], [], linestyle="-", color=clr_knockdown)
    dotted_line2 = plt.Line2D([], [], linestyle="--", color=clr_knockin, dashes=(5, 5))

    figsz = (15,10)
    fig, ((ax1, ax2, ax21),(ax3, ax31, ax4)) = plt.subplots(2, 3, figsize=figsz)

    # For debug denom/numer
    fig2, (b1, b2) = plt.subplots(1, 2, figsize=(8,4))

    # WT
    ax1.plot (x, EphAx, linestyle='-', color=clr_wt, label='$r_0$ (EphA wildtype cells)')
    ax1.plot (x, ephrinA, linestyle=':', color=clr_wt, label='$l_0$ (ephrinA)')
    # EphA3 knockin
    ax1.plot (x, EphAx+ki, linestyle='-', color=clr_knockin, label='$r_0 + ki$ (EphA3 ki)')

    ax1.legend()
    ax1.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax1.set_ylabel('Expression')
    ax1.set_ylim(0,5)
    ax1.text (0.55, 2.9, 'ki = {0}'.format(ki))

    ax2.plot (x, _epha4, linestyle=':', color=clr_wt, label='Total EphA4 expression')
    ax2.plot (x, _p_epha4, linestyle='--', color=clr_wt, label='$r_{A4}^{cis}$ (cis-bound, wildtype)')
    ax2.plot (x, EphA4_free, linestyle='-', color=clr_wt, label='$r_{A4}$  (un-bound, wildtype)')
    ax2.plot (x, EphA4_free_kd, linestyle='-', color=clr_knockdown, label='$r_{A4} - kd$ (un-bound, kd)')

    show_sum_epha = 0
    if show_sum_epha:
        ax21 = ax2.twinx()
        ax21.plot (x, EphA4_free + EphAx , linestyle='-.', color=clr_wt, label='')
        ax21.plot (x, EphA4_free + EphAx + ki, linestyle='-.', color=clr_knockin, label='$\Sigma$EphA ki')

        ax21.plot (x, EphA4_free + EphAx + ki - kd, linestyle='-.', color=clr_knockin)
        ax21.plot (x, EphA4_free + EphAx + ki - kd, linestyle='--', color=clr_knockdown, dashes=(5, 5))

        ax21.legend([filled_line_wt, filled_line_ki, (dotted_line1, dotted_line2)],
                        ['$\Sigma$EphA','$\Sigma$EphA ki', '$\Sigma$EphA ki/kd'])
        ax21.set_ylabel('Expression')

    ax2.legend()
    ax2.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax2.set_ylabel('Expression')
    ax2.set_ylim(0,5)
    ax2.text (0.6, 2, 'kd = {0}'.format(kd))

    ax21.plot (x, signal (EphAx,    clustersz_simple(EphAx, 0,  EphA4_free,    _p_epha4)), linestyle='-', color=clr_wt)
    ax21.plot (x, signal (EphAx+ki, clustersz_simple(EphAx, ki, EphA4_free,    _p_epha4)), linestyle='-', color=clr_knockin)
    ax21.plot (x, signal (EphAx,    clustersz_simple(EphAx, 0,  EphA4_free_kd, _p_epha4_kd)), linestyle='-', color=clr_knockdown)
    ax21.plot (x, signal (EphAx+ki, clustersz_simple(EphAx, ki, EphA4_free_kd, _p_epha4_kd)), linestyle='-', color=clr_knockdown)
    ax21.plot (x, signal (EphAx+ki, clustersz_simple(EphAx, ki, EphA4_free_kd, _p_epha4_kd)), linestyle='--', color=clr_knockin, dashes=(5, 5))
    ax21.legend([filled_line_wt, filled_line_ki, filled_line_kd, (dotted_line1, dotted_line2)], ['$r_0/r_{A4}$ (wildtype)','$(r_0 + ki)/r_{A4}$ (EphA3 ki)', '${r_0}/(r_{A4}-kd)$ (EphA4 kd)', '$(r_0+ki)/(r_{A4}-kd)$ (ki + kd)'])
    ax21.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax21.set_ylabel('Signal')
    yl = ax21.get_ylim()
    yl = (0, yl[1])
    ax21.set_ylim(yl)

    # Cluster size vs. EphAx/N/T posn
    # Or vs. EphAx:
    (cs, n, d) = clustersz_complex (EphAx, 0, 0.5, _p_epha4_kd)
    ax3.plot (EphAx, cs, label='$r_{A4} = 0.5$', color=clr_knockdown)
    (cs, n, d) = clustersz_complex (EphAx, 0, 1.5, _p_epha4)
    ax3.plot (EphAx, cs, label='$r_{A4} = 1.5$', color=clr_wt)
    ax3.legend()
    ax3.set_xlabel('EphAx')
    ax3.set_ylabel('Clustersize')
    yl = ax3.get_ylim()
    yl = (0, yl[1])
    ax3.set_ylim(yl)

    (cs, n, d) = clustersz_complex (EphAx, 0, 0.5, _p_epha4_kd)
    ax31.plot (x, cs, label='$r_{A4} = 0.5$', color=clr_knockdown)

    (cs, n, d) = clustersz_complex (EphAx, 0, 1.5, _p_epha4)
    ax31.plot (x, cs, label='$r_{A4} = 1.5$', color=clr_wt)

    (cs, n, d) = clustersz_complex (EphAx, ki, 0.5, _p_epha4_kd)
    ax31.plot (x, cs, label='ki EphA3, EphA4 = 0.5', color=clr_knockin)
    ax31.plot (x, cs, label='ki EphA3, EphA4 = 0.5', linestyle='--', color=clr_knockdown, dashes=(5, 5))

    (cs, n, d) = clustersz_complex (EphAx, ki, 1.5, _p_epha4)
    ax31.plot (x, cs, label='ki EphA3, EphA4 = 1.5', color=clr_knockin)


    ax31.legend([filled_line_kd, filled_line_wt, filled_line_ki, (dotted_line1, dotted_line2)],
                    ['$r_{A4}=0.5$ (kd)','$r_{A4}=1.5$ (wildtype)', '$r_{A4}=0.5$ (EphA3 ki/kd)', '$r_{A4}=1.5$ (EphA3 ki)'])
    ax31.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax31.set_ylabel('Clustersize')
    yl = ax31.get_ylim()
    yl = (0, yl[1])
    ax31.set_ylim(yl)

    (cs, n, d) = clustersz_complex(EphAx, 0,  EphA4_free,    _p_epha4)
    ax4.plot (x, signal (EphAx,    cs), linestyle='-', color=clr_wt)
    b1.plot (x, n, linestyle='-', color=clr_wt)
    b2.plot (x, d, linestyle='-', color=clr_wt)
    (cs, n, d) = clustersz_complex(EphAx, ki, EphA4_free,    _p_epha4)
    ax4.plot (x, signal (EphAx+ki, cs), linestyle='-', color=clr_knockin)
    b1.plot (x, n, linestyle='-', color=clr_knockin)
    b2.plot (x, d, linestyle='-', color=clr_knockin)
    (cs, n, d) = clustersz_complex(EphAx, 0,  EphA4_free_kd, _p_epha4_kd)
    ax4.plot (x, signal (EphAx,    cs), linestyle='-', color=clr_knockdown)
    b1.plot (x, n, linestyle='-', color=clr_knockdown)
    b2.plot (x, d, linestyle='-', color=clr_knockdown)
    (cs, n, d) = clustersz_complex(EphAx, ki, EphA4_free_kd, _p_epha4_kd)
    ax4.plot (x, signal (EphAx+ki, cs), linestyle='-', color=clr_knockdown)
    ax4.plot (x, signal (EphAx+ki, cs), linestyle='--', color=clr_knockin, dashes=(5, 5))
    b1.plot (x, n, linestyle='-', color=C.red)
    b2.plot (x, d, linestyle='-', color=C.red)
    ax4.legend([filled_line_wt, filled_line_ki, filled_line_kd, (dotted_line1, dotted_line2)], ['wildtype','EphA3 ki', 'EphA4 kd', 'ki + kd'])
    ax4.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax4.set_ylabel('Signal')
    yl = ax4.get_ylim()
    yl = (0, yl[1])
    ax4.set_ylim(yl)
    b1.set_ylabel('numerator')
    b2.set_ylabel('denominator')

    #plt.tight_layout()
    fn = 'examineRetEph_complex.svg'
    plt.savefig(fn)
    plt.show()

examineRetEph()
