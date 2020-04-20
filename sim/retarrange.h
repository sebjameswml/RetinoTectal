/*
 * Retino-tectal common code. inits all the RT projections.
 */

#include <morph/MathAlgo.h>
using morph::MathAlgo;
#include <morph/MathConst.h>
#include <morph/Random.h>
using morph::RandNormal;
#include <cmath>
using std::floor;
using std::ceil;
using std::abs;
using std::sin;
using std::cos;

#define scf(a) static_cast<Flt>(a)

template <class Flt>
class RetArrange
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
     * populated with retinal neuron somas, then figure out how to arrange NN
     * somas in rings.
     */
    //@{
    Flt ret_inner = scf(0.0);
    Flt ret_outer = scf(1.0);
    Flt ret_startangle = scf(0.0);
    Flt ret_endangle = scf(morph::TWO_PI_D);
    //! The Cartesian coordinates of the retinal neurons. This vector is of size N.
    vector<array<Flt, 2>> ret_coords;
    //! The Cartesian coordinates of the equivalent locations of the retinal neurons
    //! on the tectum. To compute these (from ret_coords) we use ellipse_a and
    //! ellipse_b.
    vector<array<Flt, 3>> tec_coords;
    //! tec_offsets = tec_coords - reg_centroids
    vector<array<Flt, 3>> tec_offsets;
    //! Store the radii which can be passed to a visualization to give a colourmap
    //! indication of the radius.
    vector<Flt> ret_coords_radii;
    //! The ring-to-ring distance; also the approximate distance between points on
    //! each ring.
    Flt ring_d = 0;
    //@}

    //! A random distribution of noise for setting up the gamma values from retinal
    //! neuron positions.
    RandNormal<Flt>* gamma_noise = (RandNormal<Flt>*)0;
    //! The standard deviation of the gamma_noise; noise to apply to the gammas. This
    //! is the sigma of a normal distribution of mean 0.
    Flt sigma_gamma = static_cast<Flt>(0.0);

    //! A random distribution of noise to add to the determination of the gradient of
    //! rho.
    RandNormal<Flt>* rho_noise = (RandNormal<Flt>*)0;
    //! Standard deviation for rho_noise
    Flt sigma_rho = static_cast<Flt>(0.0);

    //! Number of neurons. A duplicate of the RD_Barrel class's N
    unsigned int NN;
    //! Ellipse parameters ellipse_a goes in first, ellipse_b in second.
    pair<float, float> ellipse_radii;

    //! Constructor is a no-op
    RetArrange (void) {
        cout << "RetArrange constructor" << endl;
    }

    //! Deconstructor has some memory deallocation to carry out
    ~RetArrange (void) {
        if (this->sigma_gamma != scf(0.0)) {
            delete this->gamma_noise;
        }
        if (this->sigma_rho != scf(0.0)) {
            delete this->rho_noise;
        }
    }

    /*!
     * @N_ Number of neurons, copied into this->NN.
     * @e_a ellipse radius a, copied into ellipse_radii.first
     * @e_b ellipse radius b, copied into ellipse_radii.second
     */
    virtual void initRet (unsigned int N_, float e_a, float e_b) {

        cout << "RetArrange initRet()" << endl;
        this->NN = N_;
        this->ellipse_radii.first = e_a;
        this->ellipse_radii.second = e_b;

        // This populates ret_coords:
        this->arrangeRetina();
        // Set the random number generator parameters
        if (this->sigma_gamma != scf(0.0)) {
            this->gamma_noise = new RandNormal<Flt> (scf(0.0), this->sigma_gamma);
        }
        if (this->sigma_rho != scf(0.0)) {
            this->rho_noise = new RandNormal<Flt> (scf(0.0), this->sigma_rho);
        }
    }

protected:
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
     *
     * @NN 'Number of Neurons'
     */
    void arrangeRetinaInRings (void) {

        // First determine area to obtain an approximate distance between
        Flt d = scf(0); // d is our min separation
        int num = 0;
        for (d = scf(0.001); d<scf(1.0); d+=scf(0.001)) {
            // This works out number of dots on rings with a fixed d.
            num = MathAlgo::numDotsOnRings<Flt> (this->ret_inner, this->ret_outer, d,
                                                 abs(this->ret_endangle - this->ret_startangle));
            if (num <= static_cast<int>(this->NN)) {
                cout << "break on num=" << num << endl;
                break;
            }
        }
        cout << "d = " << d << " which makes "  << num << " dots" << endl;
        cout << "Need to insert N - num extras or " << this->NN << " - " << num << endl;
        int required_extras = ((int)this->NN - num);
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
            ringnums.push_back (MathAlgo::numOnCircleArc<Flt> (r, d, abs(this->ret_endangle - this->ret_startangle)));
            cout << "This ring has " << ringnums.back() << " dots on it" << endl;
            ntot += ringnums.back();
            r += d;
        }
        cout << "tlen = " << tlen << ", and ntot = " << ntot <<  endl;

        Flt len_per_extra = tlen / ((int)this->NN - (int)num);
        cout << "Insert extra every = " << len_per_extra << endl;

        int extras = (int)this->NN - (int)num;
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
        this->ret_coords.resize (this->NN);
        this->tec_coords.resize (this->NN);
        this->tec_offsets.resize (this->NN, {0,0});
        this->ret_coords_radii.resize (this->NN);
        // "ring index"
        int ri = 0;
        int ci = 0;
        // The starting position for the first dot on a ring.
        Flt d_angle_start = scf(1.0);

        this->ring_d = d;

        cout << "--------------------------" << endl;
        vector<unsigned int> dots_in_ring(ringlens.size(), 0);
        // For each ring:
        for (ri = 0; ri < static_cast<int>(ringlens.size()); ++ri) {

            // This ring has radius ri * d. Note re-use of Flt r.
            r = this->ret_inner + scf(ri) * d;

            // Number of dots in this ring
            dots_in_ring[ri] = ringextras[ri]+ringnums[ri];

            cout << "Ring r=" << r << ", dots=" << dots_in_ring[ri] << " (extras was " << ringextras[ri] << ")" << endl;

            // The angle between each dot
            Flt a = this->ret_endangle - this->ret_startangle;
            Flt d_angle = scf(0.0);
            if (a == scf(morph::TWO_PI_D)) {
                d_angle = a/scf(dots_in_ring[ri]);
                // Set up the starting angle
                d_angle_start = (d_angle_start == scf(0.0)) ? (d_angle/scf(2.0)) : scf(0.0);
            } else {
                // For NOT a full pie, i.e. for a pie slice, we can allow the neurons
                // to start and end on both edges of the pie slice.
                d_angle = a/scf(dots_in_ring[ri]-1);
                // Start all rings at 0.
                d_angle_start = this->ret_startangle;
            }

            // For each dot in the ring:
            Flt phi = d_angle_start;
            for (unsigned int dir = 0; dir < dots_in_ring[ri]; ++dir) {
                // cout << "Setting ret_coords["<<ci<<"]\n";
                this->ret_coords[ci] = { r*cos(phi), r*sin(phi) };
                // Note that I don't set set the equivalent elliptical tectal
                // coordinate here, as it needs to be scaled by the radius of the
                // outer ring of neurons.
                this->ret_coords_radii[ci] = r;
                ci++;
                phi += d_angle;
            }
        }

        // Now I loop through the rings again, setting up the variable tec_coords from
        // the values in ret_coords. A separate loop is necessary because I scale the
        // posistions on the tectum based on the radius of the outer ring of neurons,
        // rather than the radius of the circle that I originally specified.
        //
        // cout << "Before setting coordinates, r=" << r << "1/r=" << (1.0/r) << endl;
        ci = 0;
        Flt r_mult = 1.0f/r; // A factor to frig the tectal positions a bit
        for (ri = 0; ri < static_cast<int>(ringlens.size()); ++ri) {
            for (unsigned int dir = 0; dir < dots_in_ring[ri]; ++dir) {
                this->tec_coords[ci] = { r_mult*this->ellipse_radii.first*this->ret_coords[ci][0],
                                         r_mult*this->ellipse_radii.second*this->ret_coords[ci][1], scf(0.0) };
                ci++;
            }
        }

//#define DEBUG__ 1
#ifdef DEBUG__
        cout << "Tectal coordinates:" << endl;
        for (auto dot : this->tec_coords) {
            cout << dot[0] << "," << dot[1] << endl;
        }
        cout << __FUNCTION__ << ": Done." << endl;
#endif
    }

}; // RetArrange
