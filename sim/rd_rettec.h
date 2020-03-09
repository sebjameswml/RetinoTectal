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
     * slice being variable from 0 to 2PI. The slice might go all the way from r=0 to
     * r=1, but it might also be a small slice from 0 to ret_outer or from ret_inner
     * to ret_outer. This should allow me to carry out the various different
     * experimental manipulations carried out by Sperry and others.
     *
     * So, with these parameters, work out the area of the region which will be
     * populated with retinal neuron somas, then figure out how to arrange this->N
     * somas in rings.
     */
    //@{
    Flt ret_inner = scf(0.0);
    Flt ret_outer = scf(1.0);
    Flt ret_startangle = scf(0.0);
    Flt ret_endangle = scf(morph::TWO_PI_D);
    //! The Cartesian coordinates of the retinal neurons. This vector is of size N.
    vector<array<Flt, 2>> ret_coords;
    //! Store the radii which can be passed to a visualization to give a colourmap indication of the radius.
    vector<Flt> ret_coords_radii;
    //@}

    RD_RetTec (void)
        : RD_James_dncomp<Flt>() {
    }

    virtual void allocate (void) {
        RD_James_dncomp<Flt>::allocate();
    }

    virtual void init (void) {
        // Set up alpha, beta, epsilon (before RD_James_dncomp::init)
        this->setupPerNParams();
        // This populates ret_coords:
        this->arrangeRetina();
        // Now, from ret_coords, populate gammas, assuming two, orthogonal morphogen
        // fields which are noisy.
        this->setupGammas();
        // the rest of the init chain has to be called AFTER gammas are set up
        RD_James_dncomp<Flt>::init();
    }

private:

    void setupPerNParams (void) {
        // Index through thalamocortical fields, setting params:
        for (unsigned int i = 0; i < this->N; ++i) {
            this->alpha[i] = this->alpha_;
            this->beta[i] = this->beta_;
            this->epsilon[i] = this->epsilon_;
            // Sets up mask for initial branching density. This is hardcoded here.
            GaussParams<Flt> gp;
            gp.gain = 1.0;
            gp.sigma = 0.0;
            gp.x = 0.0;
            gp.y = 0.0;
            this->initmasks.push_back (gp);
        }
    }

    void setupGammas (void) {
        // Something like:
        for (unsigned int i = 0; i < this->N; ++i) {
            for (unsigned int m_idx = 0; m_idx < 2; ++m_idx) {
                cout << "Setting gamma_" << i << "," << m_idx << " to "
                     << this->G << " x " << this->ret_coords[i][m_idx] << endl;
                this->setGamma (m_idx, i, this->G * this->ret_coords[i][m_idx]);
            }
        }
        // Make a map of name to float id value
        //this->rtnames[(FLT)i/(FLT)rts.size()] = rtname.asString();
    }

    // FIXME: Send to MathAlgo
    //! How many items (dots, for example) could you arrange on a circle of
    //! radius=@radius with @d between each item's centre?
    int numOnCircle (Flt radius, Flt d) const {
        if (radius == scf(0.0)) {
            return 1;
        }
        Flt circum = (Flt)morph::TWO_PI_D * radius;
        return static_cast<int>(floor (circum / d));
    }

    // FIXME: Send to MathAlgo
    //! How many items on a circular arc of angle @a?
    int numOnCircleArc (Flt radius, Flt d, Flt a) const {
        //cout << "Called for radius == " << radius << ", d=" << d <<  endl;
        if (radius == scf(0.0)) {
            return 1;
        }
        Flt circum = scf(morph::TWO_PI_D) * radius;
        //cout << "circum = " << circum << endl;
        Flt rtn = 0;
#if 1
        // longhand, with a test for a circular arc
        if (a >= scf(morph::TWO_PI_D)) {
            rtn = floor (circum / d);
        } else {
            Flt proportion = a / scf(morph::TWO_PI_D);
            //cout << "prop = " << proportion << endl;
            Flt arclen = circum * proportion;
            //cout << "arclen = " << arclen << endl;
            rtn = floor (arclen / d);
        }
#else
        Flt proportion = a / scf(morph::TWO_PI_D);
        rtn = floor (circum * proportion / d);
#endif
        //cout << "rtn " << rtn << endl;
        return rtn;
    }

    // FIXME: Send to MathAlgo
    //! How many dots spaced by d can be placed on circular arc rings with d between them?
    int numDotsOnRings (Flt minRadius, Flt maxRadius, Flt d,
                        Flt a = scf(morph::TWO_PI_D)) const {

        // Computation of nrings differs depending on whether we have a dot and nrings, or nrings
        // from minRadius to maxRadius. Herein lies the problem!
        int n_dots = 0;
        if (minRadius == scf(0.0)) {
            int nrings = (int) floor ((maxRadius-minRadius)/d);
            if (minRadius == 0.0) {
                nrings++; // cos of centre dot.
            }

            for (int r=0; r<nrings; ++r) {
                n_dots += this->numOnCircleArc (minRadius+r*d, d, a);
            }
        } else {
            // Annulus
            int nrings = 1 + (int) floor ((maxRadius-minRadius)/d);
            cout << nrings << " rings" << endl;
            for (int r=0; r<nrings; ++r) {
                n_dots += this->numOnCircleArc (minRadius+r*d, d, a);
            }
        }
        //cout << "n_dots for d=" << d << " is " << n_dots << endl;
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
                cout << "break on num=" << num << endl;
                break;
            }
        }
        cout << "d = " << d << " which makes "  << num << " dots" << endl;
        cout << "Need to insert N - num extras or " << this->N << " - " << num << endl;
        int required_extras = ((int)this->N - num);
        cout << "Need to insert " << required_extras << " extras" << endl;
        // Right. Extras could be negative. FIXME! Or is it a BUG that it's apparently negative?
        // Don't think it's actually a bug. Just the fact that you can't keep distances equivalent
        // with one dot in the middle and a row around it with < 6 dots.

        Flt r = this->ret_inner;
        cout << "ret_inner is " << this->ret_inner << endl;
        vector<Flt> ringlens;
        vector<unsigned int> ringnums;
        unsigned int ntot = 0;
        Flt tlen = 0.0;
        // This loop puts all the lengths of the rings in ringlens adn all the numbers of dots on
        // each ring in ringnums. ntot should agree with number returned by numDotsOnRings
        while (r <= this->ret_outer) {
            Flt prop = abs(this->ret_endangle - this->ret_startangle) / scf(morph::TWO_PI_D);
            Flt alen = (prop * scf(morph::TWO_PI_D) * r);
            cout << "Ring/arc r=" << r << " has length " << alen << endl;
            ringlens.push_back (alen);
            tlen += ringlens.back();
            ringnums.push_back (this->numOnCircleArc (r, d, abs(this->ret_endangle - this->ret_startangle)));
            cout << "This ring has " << ringnums.back() << " dots on it" << endl;
            ntot += ringnums.back();
            r += d;
        }
        cout << "tlen = " << tlen << ", and ntot = " << ntot <<  endl;

        Flt len_per_extra = tlen / ((int)this->N-(int)num);
        cout << "Insert extra every = " << len_per_extra << endl;

        int extras = (int)this->N - (int)num;
        cout << "extras = " << extras << endl;
        typename vector<Flt>::iterator rli = ringlens.begin();
        cout << "Setting up ringextras with size " << ringlens.size() << endl;
        vector<unsigned int> ringextras (ringlens.size(), 0);
        typename vector<unsigned int>::iterator rei = ringextras.begin();
        Flt l = len_per_extra;
        while (extras > 0) {
            l -= *rli;
            if (l > 0.0) {
                // Then we don't insert on this ring and let len_per_extra remain a bit smaller
                rli++;
                rei++;
                // Don't drop off the last ring:
                if (rei == ringextras.end()) {
                    --rei;
                    --rli;
                }
            } else {
                // Insert an extra here.
                (*rei)++;
                --extras; // record that we added it
                cout << "added extra *rei=" << (*rei) << endl;
                // Update the remaining ring len in which to distribute extras
                *rli = -l;
                // Reset l back to len_per_extra
                l = len_per_extra;
            }
        }
        cout << "At end, l = " << l << " extras = " << extras << endl;

        // loop through ringextras and ringlens and generate the coordinates.
        this->ret_coords.resize (this->N);
        this->ret_coords_radii.resize (this->N);
        // "ring index"
        int ri = 0;
        int ci = 0;
        // The starting position for the first dot on a ring.
        Flt d_angle_start = scf(1.0);

        cout << "--------------------------" << endl;
        // For each ring:
        for (ri = 0; ri < static_cast<int>(ringlens.size()); ++ri) {

            // This ring has radius ri * d. Note re-use of Flt r.
            r = this->ret_inner + scf(ri) * d;

            // Number of dots in this ring
            unsigned int dots_in_ring = ringextras[ri]+ringnums[ri];

            cout << "Ring r=" << r << ", dots=" << dots_in_ring << " (extras was " << ringextras[ri] << ")" << endl;

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
                cout << "Setting ret_coords["<<ci<<"]\n";
                this->ret_coords[ci] = { r*cos(phi), r*sin(phi) };
                this->ret_coords_radii[ci] = r;
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
        cout << "Done." << endl;
    }

}; // RD_RetTec
