/*
 * A population model system which derives from RD_Base and implements a retino-tectal
 * system which lacks competition, but which instead incorporates a stopping
 * mechanism, wherein the diffusion and interaction with morphogens decays to zero in
 * the region where the morphogen gradient reaches a correct level.
 */
#include <morph/RD_Base.h>

#include <morph/ShapeAnalysis.h>

/*!
 * Enumerates the way that the guidance molecules are set up
 */
enum class FieldShape {
    Gauss1D,
    Gauss2D,
    Exponential1D,
    Sigmoid1D,
    Linear1D,
    CircLinear2D
};

/*!
 * A small collection of parameters to define width and location of a symmetric
 * (i.e. circular) 2D Gaussian.
 */
template <class Flt>
struct GaussParams {
    Flt gain;
    Flt sigma;
    Flt x;
    Flt y;
};

/*!
 * Reaction diffusion system. Based on Karbowski 2004, but with a removal of the Fgf8,
 * Pax6, Emx2 system, and instead an option to define several guidance molecules and
 * thalamocortical types (i.e. configurable N and M).
 *
 * This class also has a mechanism for providing normalization of the a variable.
 *
 * This is a template class using 'Flt' for the float type, this can either be single
 * precision (float) or double precision (double).
 */
template <class Flt>
class RD_RetNoComp : public morph::RD_Base<Flt>
{
public:

    /*!
     * how many retino-tectal axons are there?
     */
    alignas(Flt) unsigned int N = 5;

    /*!
     * M is the number of guidance molecules to use.
     */
    alignas(Flt) unsigned int M = 3;

    /*!
     * These are the c_i(x,t) variables from the Karb2004 paper. x is a vector in
     * two-space.
     */
    alignas(alignof(vector<vector<Flt> >))
    vector<vector<Flt> > c;

    /*!
     * To record dci/dt, as this is used in the computation of a.
     */
    alignas(alignof(vector<vector<Flt> >))
    vector<vector<Flt> > dc;

    /*!
     * These are the a_i(x,t) variables from the Karb2004 paper. x is a vector in
     * two-space. The first vector is over the different TC axon types, enumerated by
     * i, the second vector are the a_i values, indexed by the vi in the Hexes in
     * HexGrid.
     */
    alignas(alignof(vector<vector<Flt> >))
    vector<vector<Flt> > a;

    /*!
     * f is the function, unique for each of N axons (or axon groups) which is 1
     * everywhere except for the axon's destination, where it tends to 0. It
     * represents the axon's ability to sense the absolute morphogen levels of all
     * morphogens and to stop branching when morphogen levels match a set of
     * 'stopping' values for that particular axon.
     *
     * f is a scalar function of gamma_i,0, gamma_i,1 and rho(x). See notes or paper
     * for its form, which is a sigmoid-squashed 2D Gaussian. It is the same size as
     * the variable c; that is it is a vector of size N, each element of which is a
     * vector of nhex Flts.
     */
    alignas(alignof(vector<vector<Flt> >))
    vector<vector<Flt> > f;

    /*!
     * The `sharpness' of f.
     */
    alignas(Flt) Flt s = 15.0;

protected:
    /*!
     * The width of the hill defined by the function f
     */
    alignas(Flt) Flt w = 0.2; // determine by subtracting adjacent gammas?
    alignas(Flt) Flt two_w_sq = 4.0 * 0.2 * 0.2; // determine by subtracting adjacent gammas?

public:
    /*!
     * For each TC axon type, this holds the two components of the gradient field of
     * the scalar value a(x,t) (where this x is a vector in two-space)
     */
    alignas(alignof(vector<array<vector<Flt>, 2> >))
    vector<array<vector<Flt>, 2> > grad_a;

    /*!
     * Contains the chemo-attractant modifiers which are applied to a_i(x,t) in Eq 4.
     */
    alignas(alignof(vector<vector<array<vector<Flt>, 2> > >))
    vector<vector<array<vector<Flt>, 2> > > g;

    /*!
     * To hold div(g) / 3d, a static scalar field. There are M vectors of N of these
     * vectors of Flts
     */
    alignas(alignof(vector<vector<vector<Flt> > >))
    vector<vector<vector<Flt> > > divg_over3d;

    /*!
     * n(x,t) variable from the Karb2004 paper.
     */
    alignas(alignof(vector<Flt>))
    vector<Flt> n;

    /*!
     * J_i(x,t) variables - the "flux current of axonal branches of type i". This is a
     * vector field.
     */
    alignas(alignof(vector<array<vector<Flt>, 2> >))
    vector<array<vector<Flt>, 2> > J;

    /*!
     * Holds the divergence of the J_i(x)s
     */
    alignas(alignof(vector<vector<Flt> >))
    vector<vector<Flt> > divJ;

    /*!
     * The power to which a_i(x,t) is raised in Eqs 1 and 2 in the paper.
     */
    alignas(Flt) Flt k = 3.0;

protected:
    /*!
     * The diffusion parameter.
     */
    alignas(Flt) Flt D = 0.1;

    alignas(Flt) Flt twoDover3dd = this->D+this->D / 3*this->d*this->d;

public:

    /*!
     * alpha_i parameters
     */
    alignas(alignof(vector<Flt>))
    vector<Flt> alpha;
    //! To set a single alpha for all N
    alignas(Flt) Flt alpha_;

    /*!
     * beta_i parameters
     */
    alignas(alignof(vector<Flt>))
    vector<Flt> beta;
    //! To set a single beta for all N
    alignas(Flt) Flt beta_;

    /*!
     * Parameters for initial 2D Gaussian masks over the initial branching levels.
     */
    alignas(alignof(vector<GaussParams<Flt> >))
    vector<GaussParams<Flt> > initmasks;

    /*!
     * The string identifiers for each RT axon. This is of size N. Can be populated
     * from the config file. This allows me to look up the name, as given in the
     * config file, from the floating point index, obtained from (integer index / N)
     */
    alignas(alignof(map<Flt, string>)) map<Flt, string> rtnames;

protected: // We have a setter for gamma.
    /*!
     * gamma_A/B/C_i (etc) parameters from Eq 4. There are M vectors of Flts in here.
     */
    //@{
    alignas(alignof(vector<vector<Flt> >))
    vector<vector<Flt> > gamma;

    /*!
     * Used for group-based competition. One of the sets of gammas is used as a group
     * identifier.
     */
    alignas(alignof(vector<Flt>)) vector<Flt> group;
    alignas(alignof(set<Flt>)) set<Flt> groupset;
    //@}

public:

    /*!
     * Gain for setting gammas, if required.
     */
    Flt G = static_cast<Flt>(0.0);

    /*!
     * A vector of parameters for the direction of the guidance molecules. This is an
     * angle in Radians.
     */
    alignas(alignof(vector<Flt>))
    vector<Flt> guidance_phi;

    /*!
     * Guidance molecule parameters for the width of the function
     */
    alignas(alignof(vector<Flt>))
    vector<Flt> guidance_width;

    /*!
     * Width in orthogonal direction, for 2D fields.
     */
    alignas(alignof(vector<Flt>))
    vector<Flt> guidance_width_ortho;

    /*!
     * Guidance molecule parameters for the offset of the function
     */
    alignas(alignof(vector<Flt>))
    vector<Flt> guidance_offset;

    /*!
     * Guidance molecule parameters to be the gains of the functions
     */
    alignas(alignof(vector<Flt>))
    vector<Flt> guidance_gain;

    /*!
     * Guidance molecule parameters to be the time (i.e. step) at which each guidance
     * gradient is switched on.
     */
    alignas(alignof(vector<unsigned int>))
    vector<unsigned int> guidance_time_onset;

    /*!
     * Rho variables in Eq 4 - the concentrations of axon guidance molecules A, B, C,
     * etc. In Karbowski 2004, these are time independent and we will treat them as
     * such, populating them at initialisation.
     *
     * There are M vector<Flts> in rho.
     */
    //@{
    alignas(alignof(vector<vector<Flt> >))
    vector<vector<Flt> > rho;
    //@}

    /*!
     * Into grad_rho put the two components of the gradient of rho computed across the
     * HexGrid surface.
     *
     * There are M gradient fields stored in this variable.
     */
    //@{
    alignas(alignof(vector<array<vector<Flt>, 2> >))
    vector<array<vector<Flt>, 2> > grad_rho;
    //@}

    /*!
     * Memory to hold an intermediate result
     */
    alignas(alignof(vector<vector<Flt> >))
    vector<vector<Flt> > betaterm;

    /*!
     * Holds an intermediate value for the computation of Eqs 1 and 2.
     */
    alignas(alignof(vector<vector<Flt> >))
    vector<vector<Flt> > alpha_c;

    /*!
     * The contour threshold. For contour plotting [see plot_contour()], the field is
     * normalised, then the contour is plotted where the field crosses this threshold.
     */
    alignas(Flt) Flt contour_threshold = 0.5;

    alignas(Flt) Flt aNoiseGain = 0.1;
    alignas(Flt) Flt aInitialOffset = 0.8;

    /*!
     * An N element vector holding the sum of a_i for each TC type.
     */
    alignas(vector<Flt>) vector<Flt> sum_a;

    /*!
     * An N element vector holding the initial sum of a_i for each TC type.
     */
    alignas(vector<Flt>) vector<Flt> sum_a_init;

    /*!
     * Data containers for summed n, c and a.
     */
    //@{
    alignas(vector<Flt>) vector<Flt> v_nsum;
    alignas(vector<Flt>) vector<Flt> v_csum;
    alignas(vector<Flt>) vector<Flt> v_asum;
    //@}

    /*!
     * The fall-off multiplier used towards the boundary of the field
     */
    alignas(vector<Flt>) vector<Flt> bSig;

    //! A per-hex variable. Holds the ID (as a Flt) of the index into c which has the
    //! maximal value of c in each hex.
    alignas(alignof(vector<Flt>)) vector<Flt> regions;

    //! The centroids of the regions. key is the "ID" of the region - a Flt between 0
    //! and 1, with values separated by 1/N.
    map<Flt, pair<Flt, Flt> > reg_centroids;
    //! The area of each region, by Flt ID (area in number of hexes).
    map<Flt, int> region_areas;
    //! Set true when the spatial analysis has been computed
    bool spatialAnalysisComputed = false;
    /*!
     * ALIGNAS REGION ENDS.
     *
     * Below here, there's no need to worry about alignas keywords.
     */

    /*!
     * Sets the function of the guidance molecule method
     */
    vector<FieldShape> rhoMethod;

    /*!
     * Simple constructor; no arguments. Just calls RD_Base constructor
     */
    RD_RetNoComp (void)
        : morph::RD_Base<Flt>() {
    }

    /*!
     * Initialise this vector of vectors with noise. This is a model-specific
     * function.
     *
     * I apply a sigmoid to the boundary hexes, so that the noise drops away towards
     * the edge of the domain.
     */
    virtual void noiseify_vector_vector (vector<vector<Flt> >& vv, vector<GaussParams<Flt> >& gp) {
        for (unsigned int i = 0; i<this->N; ++i) {
            for (auto h : this->hg->hexen) {
                // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
                // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
                // normal value. Close to boundary, noise is less.
                vv[i][h.vi] = morph::Tools::randF<Flt>() * this->aNoiseGain + this->aInitialOffset;
                if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                    Flt bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary-this->boundaryFalloffDist)) );
                    vv[i][h.vi] = vv[i][h.vi] * bSig * gp[i].gain; // New: apply gain here (and not
                                                                   // in the Gaussian mask).
                }
            }
        }
    }

    /*!
     * Apply a mask to the noise in a vector of vectors. This masks with a 2D Gaussian
     * for each a (there are N TC type, so for each i in N, apply a different Gaussian
     * mask, probably with the same width, but different centre).
     *
     * This allows me to initialise the system in a more biologically realistic manner.
     */
    void mask_a (vector<vector<Flt> >& vv, vector<GaussParams<Flt> >& gp) {

        // Once-only parts of the calculation of the Gaussian.
        const Flt root_2_pi = 2.506628275;

        for (unsigned int i = 0; i<this->N && i < gp.size(); ++i) {

            if (!(gp[i].sigma > 0.0)) {
                continue;
            }

            Flt one_over_sigma_root_2_pi = 1 / gp[i].sigma * root_2_pi;
            Flt two_sigma_sq = 2 * gp[i].sigma * gp[i].sigma;

            for (auto h : this->hg->hexen) {

                Flt rx = gp[i].x - h.x;
                Flt ry = gp[i].y - h.y;
                Flt r = sqrt (rx*rx + ry*ry);
                // Note that the gain of the gauss (gp[i].gain) has already been
                // applied in noiseify_vector_vector()
                Flt gauss = (one_over_sigma_root_2_pi
                             * exp ( static_cast<Flt>(-(r*r))
                                     / two_sigma_sq ));
                vv[i][h.vi] *= gauss;
            }
        }
    }

    /*!
     * Perform memory allocations, vector resizes and so on.
     */
    virtual void allocate (void) {

        morph::RD_Base<Flt>::allocate();

        // Resize and zero-initialise the various containers
        this->resize_vector_vector (this->c, this->N);
        this->resize_vector_vector (this->dc, this->N);
        this->resize_vector_vector (this->a, this->N);
        this->resize_vector_vector (this->f, this->N);
        this->resize_vector_vector (this->betaterm, this->N);
        this->resize_vector_vector (this->alpha_c, this->N);
        this->resize_vector_vector (this->divJ, this->N);
        this->resize_vector_vector_vector (this->divg_over3d, this->N, this->M);

        this->resize_vector_variable (this->n);
        this->resize_vector_vector (this->rho, this->M);

        this->resize_vector_param (this->alpha, this->N);
        this->resize_vector_param (this->beta, this->N);
        this->resize_vector_param (this->group, this->N);
        this->resize_vector_vector_param (this->gamma, this->N, this->M);

        this->resize_vector_array_vector (this->grad_rho, this->M);

        // Resize grad_a and other vector-array-vectors
        this->resize_vector_array_vector (this->grad_a, this->N);
        this->resize_vector_vector_array_vector (this->g, this->N, this->M);
        this->resize_vector_array_vector (this->J, this->N);

        this->resize_vector_param (this->sum_a, this->N);
        this->resize_vector_param (this->sum_a_init, this->N);

        this->resize_vector_variable (this->bSig);

        // rhomethod is a vector of size M
        this->rhoMethod.resize (this->M);
        for (unsigned int j=0; j<this->M; ++j) {
            // Set up with Sigmoid1D as default
            this->rhoMethod[j] = FieldShape::Sigmoid1D;
        }

        // Initialise alpha, beta
        for (unsigned int i=0; i<this->N; ++i) {
            this->alpha[i] = 3;
            this->beta[i] = 3;
        }
    }

    /*!
     * Initialise variables and parameters. Carry out one-time computations required
     * of the model. This should be able to re-initialise a finished simulation as
     * well as initialise the first time.
     */
    virtual void init (void) {

        this->stepCount = 0;
        this->spatialAnalysisComputed = false;

        // Zero c and n and other temporary variables
        this->zero_vector_vector (this->c, this->N);
        //this->zero_vector_vector (this->a); // gets noisified below
        this->zero_vector_vector (this->f, this->N);
        this->zero_vector_vector (this->betaterm, this->N);
        this->zero_vector_vector (this->alpha_c, this->N);
        this->zero_vector_vector (this->divJ, this->N);
        this->zero_vector_vector_vector (this->divg_over3d, this->N, this->M);

        this->zero_vector_variable (this->n);
        this->zero_vector_vector (this->rho, this->M);

        this->zero_vector_array_vector (this->grad_rho, this->M);

        // Resize grad_a and other vector-array-vectors
        this->zero_vector_array_vector (this->grad_a, this->N);
        this->zero_vector_vector_array_vector (this->g, this->N, this->M);
        this->zero_vector_array_vector (this->J, this->N);

        this->zero_vector_variable (this->bSig);

        // Initialise a with noise
        cout << "initmasks.size(): " << this->initmasks.size() << endl;
        this->noiseify_vector_vector (this->a, this->initmasks);

        // Mask the noise off (set sigmas to 0 to ignore the masking)
        this->mask_a (this->a, this->initmasks);

        // If client code didn't initialise the guidance molecules, then do so
        if (this->guidance_phi.empty()) {
            for (unsigned int m=0; m<this->M; ++m) {
                this->guidance_phi.push_back(0.0);
            }
        }
        if (this->guidance_width.empty()) {
            for (unsigned int m=0; m<this->M; ++m) {
                this->guidance_width.push_back(1.0);
            }
        }
        if (this->guidance_width_ortho.empty()) {
            for (unsigned int m=0; m<this->M; ++m) {
                this->guidance_width_ortho.push_back(1.0);
            }
        }
        if (this->guidance_offset.empty()) {
            for (unsigned int m=0; m<this->M; ++m) {
                this->guidance_offset.push_back(0.0);
            }
        }
        if (this->guidance_gain.empty()) {
            for (unsigned int m=0; m<this->M; ++m) {
                this->guidance_gain.push_back(1.0);
            }
        }
        if (this->guidance_time_onset.empty()) {
            for (unsigned int m=0; m<this->M; ++m) {
                this->guidance_time_onset.push_back(0);
            }
        }

        for (unsigned int m=0; m<this->M; ++m) {
            if (this->rhoMethod[m] == FieldShape::Gauss1D) {
                // Construct Gaussian-waves rather than doing the full-Karbowski shebang.
                this->gaussian1D_guidance (m);

            } else if (this->rhoMethod[m] == FieldShape::Gauss2D) {
                // Construct 2 dimensional gradients
                this->gaussian2D_guidance (m);

            } else if (this->rhoMethod[m] == FieldShape::Exponential1D) {
                // Construct an 'exponential wave'
                this->exponential_guidance (m);

            } else if (this->rhoMethod[m] == FieldShape::Sigmoid1D) {
                this->sigmoid_guidance (m);

            } else if (this->rhoMethod[m] == FieldShape::Linear1D) {
                this->linear_guidance (m);

            } else if (this->rhoMethod[m] == FieldShape::CircLinear2D) {
                this->circlinear_guidance (m);
            }
        }

        // Compute gradients of guidance molecule concentrations once only
        for (unsigned int m = 0; m<this->M; ++m) {
            this->spacegrad2D (this->rho[m], this->grad_rho[m]);
        }

        // Having computed gradients, build this->g; has to be done once only (though
        // in derived classes, g may incorporate noise and thus need to be built on
        // each timestep). Note that a sigmoid is applied so that g(x) drops to zero
        // around the boundary of the domain.
        this->build_bSig();
        this->build_g();
        this->compute_divg_over3d();

        // Also compute f once only (given no noise)
        this->compute_f();

        // Now compute sum of a and record this as sum_a_init.
        for (unsigned int i = 0; i < this->N; ++i) {
            this->sum_a_computation(i);
        }
        this->sum_a_init = this->sum_a;
    }

protected:
    /*!
     * Given a RT id string @idstr, look it up in rtnames and find the Flt ID that it
     * corresponds to. Client code should have set up rtnames.
     */
    Flt rt_name_to_id (const string& idstr) {
        Flt theid = -1.0;
        typename std::map<Flt, string>::iterator rtn = this->rtnames.begin();
        while (rtn != this->rtnames.end()) {
            DBG2 ("Compare " << rtn->second << " and " << idstr << "...");
            if (rtn->second == idstr) {
                DBG ("ID string " << idstr << " matches; set theid to " << rtn->first);
                theid = rtn->first;
                break;
            }
            ++rtn;
        }
        return theid;
    }

    /*!
     * Require private setter for d. Slightly different from the base class version.
     */
    //@{
    void set_d (Flt d_) {
        morph::RD_Base<Flt>::set_d (d_);
        this->updateTwoDover3dd();
    }
    //@}

public:
    /*!
     * Public accessors for D, as it requires another attribute to be updated at the
     * same time.
     */
    //@{
    void set_D (Flt D_) {
        this->D = D_;
        this->updateTwoDover3dd();
    }
    Flt get_D (void) {
        return this->D;
    }
    //@}

protected:
    /*!
     * Compute 2D/3d^2 (and 1/3d^2 too)
     */
    void updateTwoDover3dd (void) {
        this->twoDover3dd = (this->D+this->D) / (3*this->d*this->d);
    }

public:
    /*!
     * Parameter setter methods
     */
    //@{
    /*!
     * setGamma for the guidance molecule index m_idx and the RT index n_idx to
     * @value. If group_m==m_idx, then set this->group[n_idx]=@value
     */
    int setGamma (unsigned int m_idx, unsigned int n_idx, Flt value, unsigned int group_m = 0) {
        if (gamma.size() > m_idx) {
            if (gamma[m_idx].size() > n_idx) {
                // Ok, we can set the value
                //cout << "Setting gamma[m="<<m_idx<<"][n="<<n_idx<<"] to " << value << endl;
                //cout << m_idx << "," << n_idx << "," << value << endl;
                this->gamma[m_idx][n_idx] = value;
                if (group_m == m_idx) {
                    this->group[n_idx] = value;
                    this->groupset.insert (value);
                }
            } else {
                cerr << "WARNING: DID NOT SET GAMMA (too few RT axon types for n_idx=" << n_idx << ")" << endl;
                return 1;
            }
        } else {
            cerr << "WARNING: DID NOT SET GAMMA (too few guidance molecules for m_idx=" << m_idx << ")" << endl;
            return 2;
        }
        return 0;
    }

    /*!
     * Set the w parameter
     */
    void set_w (Flt _w) {
        this->w = _w;
        // 2w * 2w:
        this->two_w_sq = 4.0 * _w * _w;
    }
    Flt get_w (void) const {
        return this->w;
    }
    //@}

    /*!
     * HDF5 file saving/loading methods.
     */
    //@{

    /*!
     * Save the c, a and n variables.
     */
    virtual void save (void) {
        stringstream fname;
        fname << this->logpath << "/c_";
        fname.width(5);
        fname.fill('0');
        fname << this->stepCount << ".h5";
        HdfData data(fname.str());
        for (unsigned int i = 0; i<this->N; ++i) {
            stringstream path;
            // The c variables
            path << "/c" << i;
            data.add_contained_vals (path.str().c_str(), this->c[i]);
            // The a variable
            path.str("");
            path.clear();
            path << "/a" << i;
            data.add_contained_vals (path.str().c_str(), this->a[i]);
            // divJ
            path.str("");
            path.clear();
            path << "/j" << i;
            data.add_contained_vals (path.str().c_str(), this->divJ[i]);
        }
        data.add_contained_vals ("/n", this->n);
    }

    void saveHG (void) {
        stringstream hgname;
        hgname << this->logpath << "/hexgrid.h5";
        this->hg->save(hgname.str().c_str());
    }

#if 1
    // Save out spatial analysis results
    void saveRegions (void) {
        stringstream fname;
        fname << this->logpath << "/regions_";
        fname.width(5);
        fname.fill('0');
        fname << this->stepCount << ".h5";
        HdfData data(fname.str());
        // Save the region centroids
        vector<Flt> keys;
        vector<Flt> x_;
        vector<Flt> y_;
        // Hopefully, this ensures that we always save N centroids, even if some
        // default to 0,0.
        for (unsigned int i = 0; i<this->N; ++i) {
            Flt k = (Flt)i/this->N;
            keys.push_back (k);
            x_.push_back (this->reg_centroids[k].first);
            y_.push_back (this->reg_centroids[k].second);

        }
        data.add_contained_vals ("/reg_centroids_id", keys);
        data.add_contained_vals ("/reg_centroids_x", x_);
        data.add_contained_vals ("/reg_centroids_y", y_);
        data.add_val ("/N", this->N);
    }
#endif

    /*!
     * Save asum, nsum and csum. Call once at end of simulation.
     */
    void savesums (void) {
        stringstream fname;
        fname << this->logpath << "/sums.h5";
        HdfData data(fname.str());
        data.add_contained_vals ("/csum", this->v_csum);
        data.add_contained_vals ("/asum", this->v_asum);
        data.add_contained_vals ("/nsum", this->v_nsum);
    }

    /*!
     * Save the guidance molecules to a file (guidance.h5)
     *
     * Also save the experimental ID map in this file, as this is something that needs
     * saving once only.
     */
    void saveGuidance (void) {
        stringstream fname;
        fname << this->logpath << "/guidance.h5";
        HdfData data(fname.str());
        for (unsigned int m = 0; m<this->M; ++m) {
            stringstream path;
            path << "/rh" << m;
            string pth(path.str());
            data.add_contained_vals (pth.c_str(), this->rho[m]);
            pth[1] = 'g'; pth[2] = 'x';
            data.add_contained_vals (pth.c_str(), this->grad_rho[m][0]);
            pth[2] = 'y';
            data.add_contained_vals (pth.c_str(), this->grad_rho[m][1]);
            for (unsigned int i = 0; i<this->N; ++i) {
                stringstream path;
                path << "/divg_" << m << "_" << i;
                string pth(path.str());
                data.add_contained_vals (pth.c_str(), this->divg_over3d[m][i]);
            }
        }
    }

    /*!
     * Computation methods
     */
    //@{

    /*!
     * Compute the values of c, the connection density
     */
    virtual void integrate_c (void) {
        // 3. Do integration of c
        for (unsigned int i=0; i<this->N; ++i) {

#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; h++) {
                // Note: betaterm used in compute_dci_dt()
                this->betaterm[i][h] = this->beta[i] * this->n[h] * static_cast<Flt>(pow (this->a[i][h], this->k));
            }

            // Runge-Kutta integration for C (or ci)
            vector<Flt> qq(this->nhex,0.);
            vector<Flt> k1 = this->compute_dci_dt (this->c[i], i);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; h++) {
                qq[h] = this->c[i][h] + k1[h] * this->halfdt;
            }

            vector<Flt> k2 = this->compute_dci_dt (qq, i);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; h++) {
                qq[h] = this->c[i][h] + k2[h] * this->halfdt;
            }

            vector<Flt> k3 = this->compute_dci_dt (qq, i);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; h++) {
                qq[h] = this->c[i][h] + k3[h] * this->dt;
            }

            vector<Flt> k4 = this->compute_dci_dt (qq, i);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; h++) {
                this->dc[i][h] = (k1[h]+2. * (k2[h] + k3[h]) + k4[h]) * this->sixthdt;
                Flt c_cand = this->c[i][h] + this->dc[i][h];
                // Avoid over-saturating c_i and make sure dc is similarly modified.
                this->dc[i][h] = (c_cand > 1.0) ? (1.0 - this->c[i][h]) : this->dc[i][h];
                this->c[i][h] = (c_cand > 1.0) ? 1.0 : c_cand;
            }
        }
    }

    /*!
     * a Computation methods
     */
    //@{

    /*!
     * A possibly normalization-function specific task to carry out
     * once after the sum of a has been computed.
     */
    virtual void sum_a_computation (const unsigned int _i) {
        // Compute the sum of a[i] across the sheet.
        this->sum_a[_i] = 0.0;
        Flt sum_tmp = 0.0;
#pragma omp parallel for reduction(+:sum_tmp)
        for (unsigned int h=0; h<this->nhex; ++h) {
            sum_tmp += this->a[_i][h];
        }
        this->sum_a[_i] = sum_tmp;
    }

    /*!
     * The normalization/transfer function.
     */
    virtual inline Flt transfer_a (const Flt& _a, const unsigned int _i) {
#ifdef NORMALIZE_TO_ONE
        // Divisive normalization to one
        Flt a_rtn = this->nhex * _a / this->sum_a[_i];
#else
        // Divisive norm with initial sum multiplier
        Flt a_rtn = this->sum_a_init[_i] * _a / this->sum_a[_i];
#endif
        // Prevent a from becoming negative, necessary only when competition is implemented:
        return (a_rtn < 0.0) ? 0.0 : a_rtn;
    }

    /*!
     * Compute the values of a, the branching density
     */
    virtual void integrate_a (void) {

        // 2. Do integration of a (RK in the 1D model). Involves computing axon branching flux.

        // Pre-compute:
        // 1) The intermediate val alpha_c.
        for (unsigned int i=0; i<this->N; ++i) {
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                this->alpha_c[i][h] = this->alpha[i] * this->c[i][h];
            }
        }

        // Runge-Kutta:
        // No OMP here - there are only N(<10) loops, which isn't
        // enough to load the threads up.
        for (unsigned int i=0; i<this->N; ++i) {

            // Runge-Kutta integration for A
            vector<Flt> qq(this->nhex, 0.0);
            this->compute_divJ (this->a[i], i); // populates divJ[i]

            vector<Flt> k1(this->nhex, 0.0);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                k1[h] = this->divJ[i][h] - this->dc[i][h];
                qq[h] = this->a[i][h] + k1[h] * this->halfdt;
            }

            vector<Flt> k2(this->nhex, 0.0);
            this->compute_divJ (qq, i);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                k2[h] = this->divJ[i][h] - this->dc[i][h];
                qq[h] = this->a[i][h] + k2[h] * this->halfdt;
            }

            vector<Flt> k3(this->nhex, 0.0);
            this->compute_divJ (qq, i);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                k3[h] = this->divJ[i][h] - this->dc[i][h];
                qq[h] = this->a[i][h] + k3[h] * this->dt;
            }

            vector<Flt> k4(this->nhex, 0.0);
            this->compute_divJ (qq, i);

#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                k4[h] = this->divJ[i][h] - this->dc[i][h];
                this->a[i][h] += (k1[h] + 2.0 * (k2[h] + k3[h]) + k4[h]) * this->sixthdt;
            }
        }
    }

    /*!
     * Compute n
     */
    virtual void compute_n (void) {

        Flt nsum = 0.0;
        Flt csum = 0.0;
#pragma omp parallel for reduction(+:nsum,csum)
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            this->n[hi] = 0;
            // First, use n[hi] so sum c over all i:
            for (unsigned int i=0; i<this->N; ++i) {
                this->n[hi] += this->c[i][hi];
            }
            // Prevent sum of c being too large:
            this->n[hi] = (this->n[hi] > 1.0) ? 1.0 : this->n[hi];
            csum += this->c[0][hi];
            // Now compute n for real:
            this->n[hi] = 1. - this->n[hi];
            nsum += this->n[hi];
        }

#ifdef DEBUG__
        if (this->stepCount % 100 == 0) {
            DBG ("System computed " << this->stepCount << " times so far...");
            DBG ("sum of n+c is " << nsum+csum);
        }
#endif
    }

    //! Sum up the integration and pass through the transfer function (i.e. the normalization)
    virtual void summation_a (void) {
        for (unsigned int i=0; i<this->N; ++i) {
            // Do any necessary computation which involves summing a here
            this->sum_a_computation (i);
            // Now apply the transfer function
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                this->a[i][h] = this->transfer_a (this->a[i][h], i);
            }
        }
    }

    //! One step of the simulation
    virtual void step (void) {
        this->stepCount++;
        // 1. Compute Karb2004 Eq 3. (coupling between connections made by each TC type)
        this->compute_n();
        // 2. Call Runge Kutta numerical integration code
        this->integrate_a();
        this->summation_a();
        this->integrate_c();
        this->spatialAnalysisComputed = false;
    }

    /*!
     * Examine the value in each Hex of the hexgrid of the scalar field f. If
     * abs(f[h]) exceeds the size of dangerThresh, then output debugging information.
     */
    void debug_values (vector<Flt>& f, Flt dangerThresh) {
        for (auto h : this->hg->hexen) {
            if (abs(f[h.vi]) > dangerThresh) {
                DBG ("Blow-up threshold exceeded at Hex.vi=" << h.vi << " ("<< h.ri <<","<< h.gi <<")" <<  ": " << f[h.vi]);
                unsigned int wait = 0;
                while (wait++ < 120) {
                    usleep (1000000);
                }
            }
        }
    }

    /*!
     * Does: f = (alpha * f) + betaterm. c.f. Karb2004, Eq 1. f is c[i] or q from the
     * RK algorithm.
     */
    vector<Flt> compute_dci_dt (vector<Flt>& f, unsigned int i) {
        vector<Flt> dci_dt (this->nhex, 0.0);
#pragma omp parallel for
        for (unsigned int h=0; h<this->nhex; h++) {
            dci_dt[h] = (this->betaterm[i][h] - this->alpha[i] * f[h]);
        }
        return dci_dt;
    }

    // Compute bSig vector, which need be carried out once only.
    void build_bSig (void) {
        for (auto h : this->hg->hexen) {
            // Sigmoid/logistic fn params: 100 sharpness, 0.02 dist offset from boundary
            this->bSig[h.vi] = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary-this->boundaryFalloffDist)) );
        }
    }

    /*!
     * Compute the function f(gamma, rho) for each axon i. See \ref{eq:gf} in the
     * notes/paper.
     */
    virtual void compute_f (void) {

        // f is computed from \vec{\gamma} and \vec{\rho} (with params w and s).
        // Components of \vec{\gamma} are this->gamma[m][i] with m=0,1.
        // Components of \vec{\rho} are (\rho_0 [+ \xi_\rho]) and (\rho_1 [+ \xi_\rho])

        // (vecrho - vecgamma)^2 = (vecrho - vecgamma) . (vecrho - vecgamma) = |vecrho-vecgamma|^2
        // for each hex:
        Flt reduce_factor = 0.4;
        for (unsigned int i=0; i<this->N; ++i) {
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                // NB: If noise is to be included, replace this->rho with (this->rho + noise)
                Flt vr_minus_vg_x = this->rho[0][h] - (reduce_factor*this->gamma[0][i]);
                Flt vr_minus_vg_y = this->rho[1][h] - (reduce_factor*this->gamma[1][i]);
                Flt vr_len_sq = vr_minus_vg_x*vr_minus_vg_x + vr_minus_vg_y*vr_minus_vg_y; // This is |vecrho-vecgamma|^2
                Flt gausshump = exp (-vr_len_sq/this->two_w_sq);
                this->f[i][h] = static_cast<Flt>(2.0) - (static_cast<Flt>(2.0)/(static_cast<Flt>(1.0) + exp (-s*gausshump)));
            }
        }
    }

    /*!
     * Build g from the gradient of rho and the gammas.
     */
    virtual void build_g (void) {

        // First zero g out
        this->zero_vector_vector_array_vector (this->g, this->N, this->M);

        for (unsigned int i=0; i<this->N; ++i) {
            for (auto h : this->hg->hexen) {
                for (unsigned int m = 0; m<this->M; ++m) {
                    this->g[m][i][0][h.vi] += (this->gamma[m][i] * this->grad_rho[m][0][h.vi]) * this->bSig[h.vi];
                    this->g[m][i][1][h.vi] += (this->gamma[m][i] * this->grad_rho[m][1][h.vi]) * this->bSig[h.vi];
                }
            }
        }
    }

    /*!
     * Compute the divergence of g and divide by 3d. Used in computation of term2 in
     * compute_divJ().
     *
     * This computation is based on Gauss's theorem.
     */
    void compute_divg_over3d (void) {

        // Change to have one for each m in M? They should then sum, right?

        for (unsigned int i = 0; i < this->N; ++i) {

#pragma omp parallel for schedule(static)
            for (unsigned int hi=0; hi<this->nhex; ++hi) {

                vector<Flt> divg(this->M, 0.0);
                // Sum up over each gradient.
                for (unsigned int m = 0; m<this->M; ++m) {
                    // First sum
                    if (HAS_NE(hi)) {
                        divg[m] += /*cos (0)*/ (this->g[m][i][0][NE(hi)] + this->g[m][i][0][hi]);
                    } else {
                        // Boundary condition _should_ be satisfied by sigmoidal roll-off of g
                        // towards the boundary, so add only g[i][0][hi]
                        divg[m] += /*cos (0)*/ (this->g[m][i][0][hi]);
                    }
                    if (HAS_NNE(hi)) {
                        divg[m] += /*cos (60)*/ 0.5 * (this->g[m][i][0][NNE(hi)] + this->g[m][i][0][hi])
                            +  (/*sin (60)*/ this->R3_OVER_2 * (this->g[m][i][1][NNE(hi)] + this->g[m][i][1][hi]));
                    } else {
                        //divg += /*cos (60)*/ (0.5 * (this->g[i][0][hi]))
                        //    +  (/*sin (60)*/ this->R3_OVER_2 * (this->g[i][1][hi]));
                    }
                    if (HAS_NNW(hi)) {
                        divg[m] += -(/*cos (120)*/ 0.5 * (this->g[m][i][0][NNW(hi)] + this->g[m][i][0][hi]))
                            +    (/*sin (120)*/ this->R3_OVER_2 * (this->g[m][i][1][NNW(hi)] + this->g[m][i][1][hi]));
                    } else {
                        //divg += -(/*cos (120)*/ 0.5 * (this->g[i][0][hi]))
                        //    +    (/*sin (120)*/ this->R3_OVER_2 * (this->g[i][1][hi]));
                    }
                    if (HAS_NW(hi)) {
                        divg[m] -= /*cos (180)*/ (this->g[m][i][0][NW(hi)] + this->g[m][i][0][hi]);
                    } else {
                        divg[m] -= /*cos (180)*/ (this->g[m][i][0][hi]);
                    }
                    if (HAS_NSW(hi)) {
                        divg[m] -= (/*cos (240)*/ 0.5 * (this->g[m][i][0][NSW(hi)] + this->g[m][i][0][hi])
                                 + ( /*sin (240)*/ this->R3_OVER_2 * (this->g[m][i][1][NSW(hi)] + this->g[m][i][1][hi])));
                    } else {
                        divg[m] -= (/*cos (240)*/ 0.5 * (this->g[m][i][0][hi])
                                 + (/*sin (240)*/ this->R3_OVER_2 * (this->g[m][i][1][hi])));
                    }
                    if (HAS_NSE(hi)) {
                        divg[m] += /*cos (300)*/ 0.5 * (this->g[m][i][0][NSE(hi)] + this->g[m][i][0][hi])
                            - ( /*sin (300)*/ this->R3_OVER_2 * (this->g[m][i][1][NSE(hi)] + this->g[m][i][1][hi]));
                    } else {
                        divg[m] += /*cos (300)*/ 0.5 * (this->g[m][i][0][hi])
                            - ( /*sin (300)*/ this->R3_OVER_2 * (this->g[m][i][1][hi]));
                    }

                    this->divg_over3d[m][i][hi] = divg[m] * this->oneover3d;
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
     * Stable with dt = 0.0001;
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

            for (unsigned int m=0 ; m<this->M; ++m) {
                if (this->stepCount >= this->guidance_time_onset[m]) {
                    // g contributes to term2
                    term2 += fa[hi] * this->divg_over3d[m][i][hi];
                    // and to term3
                    term3 += this->g[m][i][0][hi] * this->grad_a[i][0][hi] + (this->g[m][i][1][hi] * this->grad_a[i][1][hi]);
                }
            }

            this->divJ[i][hi] = (term1 - term2 - term3) * this->f[i][hi];
        }
    }

    /*!
     * Generate Gaussian profiles for the chemo-attractants.
     *
     * Instead of using the Karbowski equations, just make some gaussian 'waves'
     *
     * @m The molecule id
     */
    void gaussian1D_guidance (unsigned int m) {
        for (auto h : this->hg->hexen) {
            Flt cosphi = (Flt) cos (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt sinphi = (Flt) sin (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt x_ = (h.x * cosphi) + (h.y * sinphi);
            this->rho[m][h.vi] = guidance_gain[m] * exp(-((x_-guidance_offset[m])*(x_-guidance_offset[m])) / guidance_width[m]);
        }
    }

    /*!
     * Circular symmetric 2D Gaussian
     *
     * @m The molecule id
     */
    void gaussian2D_guidance (unsigned int m) {

        // Centre of the Gaussian is offset from 0 by guidance_offset, then rotated by
        // guidance_phi
        Flt x_ = (Flt)this->guidance_offset[m];
        Flt y_ = (Flt)0.0;

        // Rotate the initial location of the 2D Gaussian
        Flt cosphi = (Flt) cos (this->TWOPI_OVER_360 * this->guidance_phi[m]);
        Flt sinphi = (Flt) sin (this->TWOPI_OVER_360 * this->guidance_phi[m]);
        Flt x_gCentre = (x_ * cosphi) + (y_ * sinphi);
        Flt y_gCentre = - (x_ * sinphi) + (y_ * cosphi);

        for (auto h : this->hg->hexen) {

            Flt rx = x_gCentre - h.x;
            Flt ry = y_gCentre - h.y;
            Flt r = sqrt (rx*rx + ry*ry);
            this->rho[m][h.vi] = guidance_gain[m] * exp (static_cast<Flt>( -(r*r) / (2.0 * guidance_width[m])) );
        }
    }

    /*!
     * An exponential wave
     *
     * @m The molecule id
     */
    void exponential_guidance (unsigned int m) {
        for (auto h : this->hg->hexen) {
            Flt cosphi = (Flt) cos (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt sinphi = (Flt) sin (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt x_ = (h.x * cosphi) + (h.y * sinphi);
            this->rho[m][h.vi] = exp (this->guidance_gain[m] * (x_-guidance_offset[m]));
        }
    }

    /*!
     * @m The molecule id
     */
    void sigmoid_guidance (unsigned int m) {
        for (auto h : this->hg->hexen) {
            Flt cosphi = (Flt) cos (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt sinphi = (Flt) sin (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            DBG2 ("phi= " << this->guidance_phi[m] << ". cosphi: " << cosphi << " sinphi: " << sinphi);
            Flt x_ = (h.x * cosphi) + (h.y * sinphi);
            DBG2 ("x_[" << h.vi << "] = " << x_);
            this->rho[m][h.vi] = guidance_gain[m] / (1.0 + exp(-(x_-guidance_offset[m])/this->guidance_width[m]));
        }
    }

    /*!
     * @m The molecule id
     */
    void linear_guidance (unsigned int m) {
        for (auto h : this->hg->hexen) {
            Flt cosphi = (Flt) cos (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt sinphi = (Flt) sin (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt x_ = (h.x * cosphi) + (h.y * sinphi);
            this->rho[m][h.vi] = (x_-guidance_offset[m]) * this->guidance_gain[m];
        }
    }

    /*!
     * @m The molecule id
     */
    void circlinear_guidance (unsigned int m) {
        for (auto h : this->hg->hexen) {
            // Initial position is guidance_offset * cosphi/sinphi
            Flt cosphi = (Flt) cos (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt sinphi = (Flt) sin (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt x_centre = guidance_offset[m] * cosphi;
            Flt y_centre = guidance_offset[m] * sinphi;

            Flt x_ = (h.x - x_centre);
            Flt y_ = (h.y - y_centre);
            Flt r_ = sqrt(x_*x_ + y_*y_);
            this->rho[m][h.vi] = (this->guidance_gain[m] - r_) * this->guidance_gain[m];
        }
    }

#if 1
    /*!
     * Carry out any sensible spatial analysis required
     */
    virtual void spatialAnalysis (void) {

        // Don't recompute unnecessarily
        if (this->spatialAnalysisComputed == true) {
            DBG ("analysis already computed, no need to recompute.");
            return;
        }

        // Clear out previous results from an earlier timestep
        this->regions.clear();
        //this->vertices.clear(); // Not interested in a Dirichlet analysis for this
        // work, but it could be done Find regions. Based on an 'ID field'. Note that
        // although this is called dirichlet_regions, there's nothing specifically
        // Dirichlet-analysis about the function. It just finds the
        // i-that-give-the-max-of arg number 2 (this->c).
        this->regions = morph::ShapeAnalysis<Flt>::dirichlet_regions (this->hg, this->c);
        // Compute centroids of regions; used to determine aligned-ness of the barrels
        this->reg_centroids = morph::ShapeAnalysis<Flt>::region_centroids (this->hg, this->regions);

        this->spatialAnalysisComputed = true;
    }
#endif
}; // RD_RetNoComp
