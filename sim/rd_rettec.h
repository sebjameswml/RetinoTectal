/*
 * Retino-tectal system. This specialization inits all the RT projections.
 */

#include "rd_barrel.h"
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
class RD_RetTec : public RD_Barrel<Flt>
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
    //! The Cartesian coordinates of the equivalent locations of the retinal neurons
    //! on the tectum. To compute these (from ret_coords) we use ellipse_a and
    //! ellipse_b.
    vector<array<Flt, 2>> tec_coords;
    //! tec_offsets = tec_coords - reg_centroids
    vector<array<Flt, 2>> tec_offsets;
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

    RD_RetTec (void)
        : RD_Barrel<Flt>() {
    }

    ~RD_RetTec (void) {
        if (this->sigma_gamma != scf(0.0)) {
            delete this->gamma_noise;
        }
        if (this->sigma_rho != scf(0.0)) {
            delete this->rho_noise;
        }
    }

    virtual void allocate (void) {
        RD_Barrel<Flt>::allocate();
    }

    virtual void init (void) {
        // Set up alpha, beta, epsilon (before RD_Barrel::init)
        this->setupPerNParams();
        // This populates ret_coords:
        this->arrangeRetina();
        // Set the random number generator parameters
        if (this->sigma_gamma != scf(0.0)) {
            this->gamma_noise = new RandNormal<Flt> (scf(0.0), this->sigma_gamma);
        }
        if (this->sigma_rho != scf(0.0)) {
            this->rho_noise = new RandNormal<Flt> (scf(0.0), this->sigma_rho);
        }
        // Now, from ret_coords, populate gammas, assuming two, orthogonal morphogen
        // fields which are noisy.
        this->setupGammas();
        // the rest of the init chain has to be called AFTER gammas are set up
        RD_Barrel<Flt>::init();
    }

protected:

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
        for (unsigned int i = 0; i < this->N; ++i) {
            for (unsigned int m_idx = 0; m_idx < 2; ++m_idx) {
                // With noise on gammas?
                if (this->sigma_gamma == scf(0.0)) {
                    this->setGamma (m_idx, i, this->G * (this->ret_coords[i][m_idx]));
                } else {
                    this->setGamma (m_idx, i, this->G * (this->ret_coords[i][m_idx] + this->gamma_noise->get()));
                }
            }
        }
        // Make a map of name to float id value
        //this->rtnames[(FLT)i/(FLT)rts.size()] = rtname.asString();
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
            num = MathAlgo<Flt>::numDotsOnRings (this->ret_inner, this->ret_outer, d,
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
            ringnums.push_back (MathAlgo<Flt>::numOnCircleArc (r, d, abs(this->ret_endangle - this->ret_startangle)));
            cout << "This ring has " << ringnums.back() << " dots on it" << endl;
            ntot += ringnums.back();
            r += d;
        }
        cout << "tlen = " << tlen << ", and ntot = " << ntot <<  endl;

        Flt len_per_extra = tlen / ((int)this->N - (int)num);
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
        this->tec_coords.resize (this->N);
        this->tec_offsets.resize (this->N, {0,0});
        this->ret_coords_radii.resize (this->N);
        // "ring index"
        int ri = 0;
        int ci = 0;
        // The starting position for the first dot on a ring.
        Flt d_angle_start = scf(1.0);

        this->ring_d = d;

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
                d_angle_start = this->ret_startangle;
            }

            // For each dot in the ring:
            Flt phi = d_angle_start;
            for (unsigned int dir = 0; dir < dots_in_ring; ++dir) {
                //cout << "Setting ret_coords["<<ci<<"]\n";
                this->ret_coords[ci] = { r*cos(phi), r*sin(phi) };
                // Also set the equivalent elliptical tectal coordinate
                this->tec_coords[ci] = { r*this->ellipse_a*cos(phi), r*this->ellipse_b*sin(phi) };
                this->ret_coords_radii[ci] = r;
                ci++;
                phi += d_angle;
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

public:
    //! Override step() to recompute g on each sim step
    virtual void step (void) {
        this->stepCount++;
        // 0. Rebuild g, with a new set of random samples, if necessary
        if (this->sigma_rho != scf(0.0)) {
            this->build_g();
            this->compute_divg_over3d();
        }
        // 1. Compute Karb2004 Eq 3. (coupling between connections made by each TC type)
        this->compute_n();
        // 2. Call Runge Kutta numerical integration code
        this->integrate_a();
        this->summation_a();
        this->integrate_c();
        this->spatialAnalysisComputed = false;
    }

protected:
    /*!
     * Build g from the gradient of rho and the gammas, with an option for a noisy calculation.
     */
    void build_g (void) {

        // First zero g.
        this->zero_vector_vector_array_vector (this->g, this->N, this->M);

        if (this->sigma_rho == scf(0.0)) {
            // No noise
            for (unsigned int i=0; i<this->N; ++i) {
                for (auto h : this->hg->hexen) {
                    for (unsigned int m = 0; m<this->M; ++m) {
                        this->g[m][i][0][h.vi] += (this->gamma[m][i]
                                                   * this->grad_rho[m][0][h.vi]
                                                   * this->bSig[h.vi]);
                        this->g[m][i][1][h.vi] += (this->gamma[m][i]
                                                   * this->grad_rho[m][1][h.vi]
                                                   * this->bSig[h.vi]);
                    }
                }
            }
        } else {
            // Adds noise
            for (unsigned int i=0; i<this->N; ++i) {
                for (auto h : this->hg->hexen) {
                    for (unsigned int m = 0; m<this->M; ++m) {
                        this->g[m][i][0][h.vi] += (this->gamma[m][i]
                                                   * (this->rho_noise->get()
                                                      + this->grad_rho[m][0][h.vi])
                                                   * this->bSig[h.vi]);
                        this->g[m][i][1][h.vi] += (this->gamma[m][i]
                                                   * (this->rho_noise->get()
                                                      + this->grad_rho[m][1][h.vi])
                                                   * this->bSig[h.vi]);
                    }
                }
            }
        }
    }

    /*!
     * Computes the "flux of axonal branches" term, J_i(x) (Eq 4)
     *
     * Inputs: this->g, fa (which is this->a[i] or a q in the RK algorithm), this->D,
     * @a i, the RT type.  Helper functions: spacegrad2D().  Output: this->divJ
     *
     * This is an overload of the instance in rd_barrel.h. It adds noise to the
     * direction of the contribution provided by the static vector field g (which is
     * due to the molecular signalling gradients on the tectum).
     */
    virtual void compute_divJ (vector<Flt>& fa, unsigned int i) {

        // Compute gradient of a_i(x), for use computing the third term, below.
        this->spacegrad2D (fa, this->grad_a[i]);

        // Three terms to compute; see Eq. 17 in methods_notes.pdf
#pragma omp parallel for //schedule(static) // This was about 10% faster than schedule(dynamic,50).
        for (unsigned int hi=0; hi<this->nhex; ++hi) {

            // 1. The D Del^2 a_i term. Eq. 18.
            // Compute the sum around the neighbours
            Flt thesum = -6 * fa[hi];

            thesum += fa[(HAS_NE(hi)  ? NE(hi)  : hi)];
            thesum += fa[(HAS_NNE(hi) ? NNE(hi) : hi)];
            thesum += fa[(HAS_NNW(hi) ? NNW(hi) : hi)];
            thesum += fa[(HAS_NW(hi)  ? NW(hi)  : hi)];
            thesum += fa[(HAS_NSW(hi) ? NSW(hi) : hi)];
            thesum += fa[(HAS_NSE(hi) ? NSE(hi) : hi)];

            // Multiply sum by 2D/3d^2 to give term1
            Flt term1 = this->twoDover3dd * thesum;

            // 2. The (a div(g)) term.
            Flt term2 = 0.0;

            // 3. Third term is this->g . grad a_i. Should not contribute to J, as
            // g(x) decays towards boundary.
            Flt term3 = 0.0;

            // In previous code, divg_over3d was computed once only, at the start of
            // the sim. Here, we need to compute it on every step.
            for (unsigned int m =0 ; m < this->M; ++m) {
                if (this->stepCount >= this->guidance_time_onset[m]) {
                    // g contributes to term2
                    term2 += fa[hi] * this->divg_over3d[m][i][hi];
                    // and to term3
                    term3 += this->g[m][i][0][hi] * this->grad_a[i][0][hi] + (this->g[m][i][1][hi] * this->grad_a[i][1][hi]);
                }
            }

            this->divJ[i][hi] = term1 - term2 - term3;
        }
    }

public:
    /*!
     * Overrides RD_Barrel spatial analysis, adding the tec_coords-reg_centroids
     * computation
     */
    virtual void spatialAnalysis (void) {

        // Don't recompute unnecessarily
        if (this->spatialAnalysisComputed == true) {
            DBG ("analysis already computed, no need to recompute.");
            return;
        }

        // Clear out previous results from an earlier timestep
        this->regions.clear();
        this->regions = morph::ShapeAnalysis<Flt>::dirichlet_regions (this->hg, this->c);
        this->reg_centroids = morph::ShapeAnalysis<Flt>::region_centroids (this->hg, this->regions);

        // Can now do diffs between reg_centroids and tec_coords, placing the
        // resulting vectors in tec_offsets. To get a scalar value for the pattern's
        // match to the origin pattern, can simply sum the squares of the tec_offset
        // vector lengths.
        for (unsigned int i = 0; i < this->N; ++i) {
            array<Flt, 2> tc = this->tec_coords[i];
            pair<Flt, Flt> rc = this->reg_centroids[(Flt)i/(Flt)this->N];
            //cout << "Comparing Tectal coordinate (" << tc[0] << "," << tc[1]
            //     << ") with centroid coord (" << rc.first << "," << rc.second << endl;
            array<Flt, 2> vec = tc;
            vec[0] -= rc.first;
            vec[1] -= rc.second;
            // tec_offsets contains vectors pointing FROM reg_centroids TO tc_coords
            this->tec_offsets[i] = vec;
        }

        this->spatialAnalysisComputed = true;
    }

    void saveSpatial (void) {
        stringstream fname;
        fname << this->logpath << "/spatial_";
        fname.width(5);
        fname.fill('0');
        fname << this->stepCount << ".h5";
        HdfData data(fname.str());

        data.add_contained_vals ("/tec_coords", this->tec_coords);
        // Need to implement map<Flt, pair<Flt,Flt>> for this one:
        //data.add_contained_vals ("/reg_centroids", this->reg_centroids);
        data.add_contained_vals ("/tec_offsets", this->tec_offsets);

        // "identified regions" here (as same size as n, c etc)
        data.add_contained_vals ("/dr", this->regions);
    }

}; // RD_RetTec
