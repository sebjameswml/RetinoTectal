/*
 * Retino-tectal system. This specialization inits all the RT projections.
 */

#include "rd_james_dncomp.h"
#include <morph/MathConst.h>
#include <cmath>
using std::floor;
using std::ceil;
using std::abs;

template <class Flt>
class RD_RetTec : public RD_James_dncomp<Flt>
{
public:
    /*!
     * Initially, I will assume a uniformly distributed retina, with no fovea. In other words, the
     * mean density of retinal neurons is independent of position on the retina.
     *
     * I don't think we need to care about the size of the retina, it has radius 1. But, we do need
     * to decide how to arrange the N origins for the RT axons. They might be arranged in some sort
     * of pie slice, with the angle of the slice being variable from 0 to 2PI. The slight might go
     * all the way from r=0 to r=1, but it might also be a small slice from 0 to ret_outer or from
     * ret_inner to ret_outer. This should allow me to carry out the various different experimental
     * manipulations carried out by Sperry and others.
     *
     * So, with these parameters, work out the area of the region which will be populated with
     * retinal neuron somas, then figure out how to arrange this->N somas in rings.
     */
    //@{
    Flt ret_inner = static_cast<Flt>(0.0);
    Flt ret_outer = static_cast<Flt>(1.0);
    Flt ret_angle = static_cast<Flt>(morph::TWO_PI_D);
    //@}

    RD_RetTec (void)
        : RD_James_dncomp<Flt>() {
    }

    virtual void allocate (void) {
        RD_James_dncomp<Flt>::allocate();
    }

    virtual void init (void) {
        RD_James_dncomp<Flt>::init();
        this->arrangeRetina();
    }

private:
    //! How many items (dots, for example) could you arrange on a circle of radius=@radius with @d
    //! between each item's centre?
    int numOnCircle (Flt radius, Flt d) const {
        if (radius == static_cast<Flt>(0.0)) {
            return 1;
        }
        Flt circum = (Flt)morph::TWO_PI_D * radius;
        return static_cast<int>(floor (circum / d));
    }

    //! How many items on a circular arc of angle @a?
    int numOnCircleArc (Flt radius, Flt d, Flt a) const {
        if (radius == static_cast<Flt>(0.0)) {
            return 1;
        }
        Flt circum = (Flt)morph::TWO_PI_D * radius;
        Flt numf = floor (circum / d);
        Flt proportion = morph::TWO_PI_D / a;
        return (int)(numf * proportion);
    }

    //! How many dots spaced by d can be placed on circlular rings with d between them?
    int numDotsOnRings (Flt minRadius, Flt maxRadius, Flt d,
                        Flt startAngle = (Flt)0.0, Flt endAngle = (Flt)morph::TWO_PI_D) const {

        int nrings = (int) floor ((maxRadius-minRadius)/d);
        if (minRadius == 0.0) {
            nrings++; // cos of centre dot.
        }

        int n_dots = 0;
        for (int r=0; r<nrings; ++r) {
            n_dots += this->numOnCircleArc (r*d, d, abs(endAngle-startAngle));
        }
        //cout << "n_dots for d=" << d << " is " << n_dots << endl;
        cout << d << "," << n_dots << endl;
        return n_dots;
    }

    int costFunction (int num, Flt d) {
        int diff = num - numDotsOnRings (this->ret_inner, this->ret_outer, d);
        int diff_abs = abs(diff);
        if (diff < 0) {
            // Then we got MORE dots than num on the rings. Bad.
        }
        return diff;
    }

    void arrangeRetina (void) {

        // First determine area to obtain an approximate distance between
        Flt d = 0.0; // d is our min separation
        for (d = 0.501; d<=0.502; d+=0.001) {
            // This works out number of dots on rings keeping d really rigid. What I want to do is
            // to allow d to vary slightly on some of the rings to allow any integer number of dots
            // to be arranged.
            this->numDotsOnRings (0.0, 1.0, d);
        }
    }

}; // RD_RetTec
