/*
 * Retino-tectal system. This specialization inits all the RT projections.
 */

#include "rd_james_dncomp.h"
#include <morph/MathConst.h>
#include <cmath>
using std::floor;
using std::ceil;
using std::abs;
using std::sin;
using std::cos;

#define scf(a) static_cast<Flt>(a)

template <class Flt>
class RD_RetTec : public RD_James_dncomp<Flt>
{
public:
    /*!
     * Initially, I will assume a uniformly distributed retina, with no fovea. In
     * other words, the mean density of retinal neurons is independent of position on
     * the retina.
     *
     * I don't think we need to care about the size of the retina, it has radius
     * 1. But, we do need to decide how to arrange the N origins for the RT
     * axons. They might be arranged in some sort of pie slice, with the angle of the
     * slice being variable from 0 to 2PI. The slight might go all the way from r=0 to
     * r=1, but it might also be a small slice from 0 to ret_outer or from ret_inner
     * to ret_outer. This should allow me to carry out the various different
     * experimental manipulations carried out by Sperry and others.
     *
     * So, with these parameters, work out the area of the region which will be
     * populated with retinal neuron somas, then figure out how to arrange this->N
     * somas in rings.
     */
    //@{
    Flt ret_inner = scf(0.3);
    Flt ret_outer = scf(1.0);
    Flt ret_startangle = scf(0.0);
    Flt ret_endangle = scf(morph::TWO_PI_D/4.0);
    //! The Cartesian coordinates of the retinal neurons. This vector is of size N.
    vector<array<Flt, 2>> ret_coords;
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
    //! How many items (dots, for example) could you arrange on a circle of
    //! radius=@radius with @d between each item's centre?
    int numOnCircle (Flt radius, Flt d) const {
        if (radius == scf(0.0)) {
            return 1;
        }
        Flt circum = (Flt)morph::TWO_PI_D * radius;
        return static_cast<int>(floor (circum / d));
    }

    //! How many items on a circular arc of angle @a?
    int numOnCircleArc (Flt radius, Flt d, Flt a) const {
        if (radius == scf(0.0)) {
            return 1;
        }
        Flt circum = scf(morph::TWO_PI_D) * radius;
        Flt numf = floor (circum / d);
        Flt proportion = a / scf(morph::TWO_PI_D);
        return (int)(numf * proportion);
    }

    //! How many dots spaced by d can be placed on circular arc rings with d between them?
    int numDotsOnRings (Flt minRadius, Flt maxRadius, Flt d,
                        Flt a = scf(morph::TWO_PI_D)) const {

        int nrings = (int) floor ((maxRadius-minRadius)/d);
        if (minRadius == 0.0) {
            nrings++; // cos of centre dot.
        }

        int n_dots = 0;
        for (int r=0; r<nrings; ++r) {
            n_dots += this->numOnCircleArc (r*d, d, a);
        }
        //cout << "n_dots for d=" << d << " is " << n_dots << endl;
        //cout << d << "," << n_dots << endl;
        return n_dots;
    }

    /*!
     * Arrange retinal neurons in some simple model. This computes coordinates for the
     * neurons from which interaction paramters can then be computed with or without
     * stochasticity.
     */
    void arrangeRetina (void) {
        this->arrangeRetinaInRings();
    }

    /*!
     * Retinal arrangement in concentric rings or arcs of concentric rings.
     */
    void arrangeRetinaInRings (void) {

        // First determine area to obtain an approximate distance between
        Flt d = scf(0); // d is our min separation
        int num = 0;
        for (d = scf(0.001); d<scf(1.0); d+=scf(0.001)) {
            // This works out number of dots on rings with a fixed d.
            num = this->numDotsOnRings (this->ret_inner, this->ret_outer, d,
                                        abs(this->ret_endangle - this->ret_startangle));
            if (num <= static_cast<int>(this->N)) {
                break;
            }
        }
        //cout << "d = " << d << " which makes "  << num << " dots" << endl;
        //cout << "Need to insert " << (this->N - num) << " extras" << endl;

        Flt r = d;
        vector<Flt> ringlens;
        vector<unsigned int> ringnums;
        unsigned int ntot = 0;
        if (this->ret_inner = scf(0.0)) {
            // First ring is radius 0, so length 0.
            ringlens.push_back (scf(0));
            // First ring num is 1, as it's a single dot
            ringnums.push_back (1);
            ntot = 1; // to count the first, centre dot
        }
        Flt tlen = 0.0;
        while (r <= 1.0) {
            //cout << "Ring/arc r=" << r << " has circumference " << (morph::TWO_PI_D * r) << endl;
            ringlens.push_back (scf(morph::TWO_PI_D) * r);
            tlen += ringlens.back();
            ringnums.push_back (this->numOnCircleArc (r, d, abs(this->ret_endangle - this->ret_startangle)));
            ntot += ringnums.back();
            r += d;
        }
        //cout << "tlen = " << tlen << ", and ntot = " << ntot <<  endl;

        Flt len_per_extra = tlen / (this->N-num);
        //cout << "Insert extra every = " << len_per_extra << endl;

        unsigned int extras = this->N - (unsigned int)num;
        typename vector<Flt>::iterator rli = ringlens.begin();
        //rli++; // Skip the zeroth
        vector<unsigned int> ringextras (ringlens.size(), 0);
        typename vector<unsigned int>::iterator rei = ringextras.begin();
        //rei++; // Skip the zeroth
        Flt l = len_per_extra;
        while (extras) {
            l -= *rli;
            if (l > 0.0) {
                // Then we don't insert on this ring and let len_per_extra remain a bit smaller
                rli++;
                rei++;
            } else {
                // Insert an extra here.
                (*rei)++;
                --extras; // record that we added it
                // Update the remaining ring len in which to distribute extras
                *rli = -l;
                // Reset l back to len_per_extra
                l = len_per_extra;
            }
        }

        // loop through ringextras and ringlens and generate the coordinates.
        this->ret_coords.resize (this->N);
        // "ring index"
        int ri = 0;
        int ci = 0;
        // Insert the first dot at the centre. ci is "coordinate index"
        //this->ret_coords[ci] = { scf(0.0), scf(0.0) };

        // The starting position for the first dot on a ring.
        Flt d_angle_start = scf(1.0);

        // For each ring:
        for (ri = 0; ri < static_cast<int>(ringlens.size()); ++ri) {

            // This ring has radius ri * d. Note re-use of Flt r.
            r = scf(ri) * d;

            // Number of dots in this ring
            unsigned int dots_in_ring = ringextras[ri]+ringnums[ri];

            cout << "Ring r=" << r << ", dots=" << dots_in_ring << endl;

            // The angle between each dot
            Flt a = this->ret_endangle - this->ret_startangle;
            Flt d_angle = scf(0.0);
            if (a == scf(morph::TWO_PI_D)) {
                d_angle = a/scf(dots_in_ring);
                // Set up the starting angle
                d_angle_start = (d_angle_start == scf(0.0)) ? (d_angle/scf(2.0)) : scf(0.0);
            } else {
                // For NOT a full pie, i.e. for a pie slice, we can allow the neurons
                // to start and end on both edges of the pie slice.
                d_angle = a/scf(dots_in_ring-1);
                // Start all rings at 0.
                d_angle_start = scf(0.0);
            }

            // For each dot in the ring:
            Flt phi = d_angle_start;
            for (unsigned int dir = 0; dir < dots_in_ring; ++dir) {
                this->ret_coords[ci] = { r*cos(phi), r*sin(phi) };
                ci++;
                phi += d_angle;
            }

        }
#define DEBUG__ 1
#ifdef DEBUG__
        cout << "Coordinates:" << endl;
        for (auto dot : this->ret_coords) {
            cout << dot[0] << "," << dot[1] << endl;
        }
#endif
    }

}; // RD_RetTec
