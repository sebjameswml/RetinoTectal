#
# This script explores clustersize as a function of both EphAx, and EphA4
# WIP on thresholded r2 collapse.
#

import numpy as np
import matplotlib.pyplot as plt
from sebcolour import Colour as C
clr_knockin = C.darkorchid1
clr_knockdown = C.cyan2
clr_knockdown2 = C.midnightblue
clr_wt = C.yellowgreen
clr_blk = C.black
clr_red = C.crimson
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

## Computation of a cluster size, simplest scheme
def clustersz_simple (_EphAx, ki, _EphA4, _EphA4_cisbound):
    cs = np.ones (len(_EphAx)) / np.power(_EphA4, 1)
    return cs

## The signal function.
def signal (rcpt, cs):
    return rcpt * cs

##
## Epha3/Epha4
##
def examineRetEph():

    x = np.linspace(0,1,21)
    print ('x = {0}'.format(x))
    epha4_constant = 3.5

    # Knockin and knockdown, which should be copied from e_eph_ki-wt.json and e_eph_ki-kd.json (
    ki = .8
    kd = 2.1
    kiki = 3.2

    ## Binding affinity for EphAx. prop. to 1/K_D. see Monschau et al
    ## 1997 for dissociation constants K_D = 6.16e-10 for ephrinA5/EphA5
    ## and 1.44e-10 for ephrinA5/EphA3. Call it 3.8e-10 for EphAx,
    ## ignore the magnitude and let w_EphAx be the single parameter (value about 0.2).

    w_EphAx = 0.25 ## This must also be set in m_ee_GJ_best_1_EphA4_r2collapse.json

    ## Binding affinity parameter for EphA4. K_D = 6.22e-10 for
    ## ephrinA5, so binds slightly less well than EphA3/5.
    ## 0.611 = 3.8/6.22 = K_D(Ax)/K_D(A4)

    w_EphA4 = w_EphAx * 0.611

    ## Thresholds for EphA4 mechanisms
    h_A4 = 1.1626
    h_0 = 2

    B=0.26
    C=2.3
    D=1.05

    _exp = B * np.exp(C*x) + D # T exponential_expression (const T& x)

    ## Start with a certain amount of ephrinA (ephrinA5 in retina - Frisen et al)
    ephrinA = np.flip(_exp)

    ## Choose between simple knockdown, where we just take the free EphA4 and
    ## subtract a scalar, vs knocking down the original amount of EPhA4 and
    ## re-multiplying by retinal ephrinA level
    simple_knockdown = 0

    ## Define an even expression of EphA4 (No nasal-temporal gradient)
    _epha4 = np.ones(len(ephrinA)) * epha4_constant
    _epha4_kd = _epha4 - kd
    print ('_epha4 = {0}'.format(_epha4))

    ## Phosphorylised EphA4
    _p_epha4    = (w_EphA4 * _epha4    * ephrinA)
    print ('_p_epha4 = {0}'.format(_p_epha4))

    ## Remaining epha4 could interact
    EphA4_free = _epha4 - _p_epha4
    print ('EphA4_free = {0}'.format(EphA4_free))


    ## Let EphA4 interact with cis ephrinA (Hornberger et al 1999)
    if simple_knockdown:
        EphA4_free_kd = EphA4_free - kd
        _p_epha4_kd = _p_epha4 + kd
    else:
        # Knockdown static level of EphA4
        _p_epha4_kd = (w_EphA4 * _epha4_kd * ephrinA)
        print ('_epha4_kd = {0}'.format(_epha4_kd))
        print ('_p_epha4_kd = {0}'.format(_p_epha4_kd))
        EphA4_free_kd = _epha4_kd - _p_epha4_kd

    ## And some expression of EphA3/x whatever
    EphAx = _exp

    # Have to create a curated legend for some plots (those with dashed coloured line)
    filled_line_wt = plt.Line2D([], [], linestyle="-", color=clr_wt)
    filled_line_ki = plt.Line2D([], [], linestyle="-", color=clr_knockin)
    filled_line_kd = plt.Line2D([], [], linestyle="-", color=clr_knockdown)
    dotted_line1 = plt.Line2D([], [], linestyle="-", color=clr_knockdown)
    dotted_line2 = plt.Line2D([], [], linestyle="--", color=clr_knockin, dashes=(5, 5))

    figsz = (15,4)
    fig, (ax1, ax2, ax21) = plt.subplots(1, 3, figsize=figsz)

    # WT
    ax1.plot (x, EphAx, linestyle='-', color=clr_wt, label='$\Sigma$EphA ($r_0$)')
    ax1.plot (x, ephrinA, linestyle=':', color=clr_wt, label='ephrinA ($l_0)$ ')
    # EphA3 knockin
    ax1.plot (x, EphAx+ki, linestyle='-', color=clr_knockin, label='EphA3 knock-in ($r_0 + ki$)')
    ax1.plot ([0,1], [h_0, h_0], linestyle=':', color=clr_red, label='$\Sigma$EphA threshold ($h_0$)')

    # Could look at cis/free EphAx in future, but avoid this complexity at present
    study_EphAx_binding = 0
    if study_EphAx_binding:
        EphAx_cis = w_EphAx * EphAx * ephrinA
        ax1.plot (x, EphAx_cis, linestyle='--', color=clr_wt, label='$\Sigma$EphA cis-bound')
        ax1.plot (x, EphAx - EphAx_cis, linestyle='-.', color=clr_wt, label='$\Sigma$EphA free')
        EphAx_ki_cis = w_EphAx * (EphAx+ki) * ephrinA
        ax1.plot (x, EphAx_ki_cis, linestyle='--', color=clr_knockin, label='$\Sigma$EphA ki cis-bound')
        ax1.plot (x, (EphAx+ki) - EphAx_ki_cis, linestyle='-.', color=clr_knockin, label='$\Sigma$EphA ki free')


    ax1.legend()
    ax1.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax1.set_ylabel('Expression')
    ax1.set_ylim(0,5)
    x0,x1 = ax1.get_xlim()
    y0,y1 = ax1.get_ylim()
    ax1.set_aspect(abs(x1-x0)/abs(y1-y0))
    ax1.text (0.52, 2.67, 'ki = {0}'.format(ki))

    ax2.plot (x, _epha4, linestyle=':', color=clr_wt, label='Total EphA4, wildtype ($r_{A4}$)')
    ax2.plot (x, _epha4_kd, linestyle=':', color=clr_knockdown, label='Total EphA4, kd ($r_{A4}^{kd}$)')
    ax2.plot (x, _p_epha4, linestyle='--', color=clr_wt, label='$r_{A4}^{cis}$ (cis-bound, wt)')
    ax2.plot (x, _p_epha4_kd, linestyle='--', color=clr_knockdown, label='$r_{A4}^{cis,kd}$ (cis-bound, kd)')
    ax2.plot (x, EphA4_free, linestyle='-', color=clr_wt, label='$r_{A4}^{free}$  (un-bound, wt)')
    ax2.plot (x, EphA4_free_kd, linestyle='-', color=clr_knockdown, label='$r_{A4}^{free,kd}$ (un-bound, kd)')
    print ('EphA4_free_kd = {0}'.format(EphA4_free_kd))

    #if not simple_knockdown:
    #    ax2.plot (x, EphA4_free_kd2, linestyle='-', color=clr_knockdown2, label='$r_{A4} - kd$ (un-bound, kd/kd)')
    ax2.plot ([0,1], [h_A4, h_A4], linestyle=':', color=clr_red, label='EphA4 threshold ($h_{A4}$)')

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
    x0,x1 = ax2.get_xlim()
    y0,y1 = ax2.get_ylim()
    ax2.set_aspect(abs(x1-x0)/abs(y1-y0))
    ax2.text (0.6, 2, 'kd = {0}'.format(kd))

    ax21.plot (x, signal (EphAx,    clustersz_simple(EphAx, 0,  EphA4_free,    _p_epha4)), linestyle='-', color=clr_wt)
    ax21.plot (x, signal (EphAx+ki, clustersz_simple(EphAx, ki, EphA4_free,    _p_epha4)), linestyle='-', color=clr_knockin)
    ax21.plot (x, signal (EphAx,    clustersz_simple(EphAx, 0,  EphA4_free_kd, _p_epha4_kd)), linestyle='-', color=clr_knockdown)
    ax21.plot (x, signal (EphAx+ki, clustersz_simple(EphAx, ki, EphA4_free_kd, _p_epha4_kd)), linestyle='-', color=clr_knockdown)
    ax21.plot (x, signal (EphAx+ki, clustersz_simple(EphAx, ki, EphA4_free_kd, _p_epha4_kd)), linestyle='--', color=clr_knockin, dashes=(5, 5))

    if study_EphAx_binding:
        ax21.plot (x, signal (EphAx - EphAx_cis,    clustersz_simple(EphAx - EphAx_cis, 0,  EphA4_free,    _p_epha4)), linestyle='-.', color=clr_wt)
        ax21.plot (x, signal ((EphAx+ki) - EphAx_ki_cis,    clustersz_simple(EphAx+ki - EphAx_ki_cis, 0,  EphA4_free,    _p_epha4)), linestyle='-.', color=clr_knockin)

    ax21.legend([filled_line_wt, filled_line_ki, filled_line_kd, (dotted_line1, dotted_line2)], ['$r_0/r_{A4}^{free}$ (wildtype)','$(r_0 + ki)/r_{A4}^{free}$ (EphA3 ki)', '${r_0}/(r_{A4}^{free,kd})$ (EphA4 kd)', '$(r_0+ki)/(r_{A4}^{free,kd})$ (ki + kd)'])
    ax21.set_xlabel('N {0} retina {0} T'.format(u"\u27f6"))
    ax21.set_ylabel('Signal')
    yl = ax21.get_ylim()
    yl = (0, yl[1])
    ax21.set_ylim(yl)
    x0,x1 = ax21.get_xlim()
    y0,y1 = ax21.get_ylim()
    ax21.set_aspect(abs(x1-x0)/abs(y1-y0))

    #plt.tight_layout()
    fn = 'examineRetEph.svg'
    plt.savefig(fn)
    plt.show()

examineRetEph()
