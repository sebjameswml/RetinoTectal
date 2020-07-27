/*
 * A reaction diffusion system which derives from RD_Base and implements the *comp2*
 * BarrelEmerge model. This merges together code from RD_James and RD_James_comp2 (see
 * https://github.com/ABRG-Models/BarrelEmerge) into this single file.
 *
 * This is used as an axon guidance model which makes use of axon-axon competition for
 * branching.
 */
#include <morph/RD_Base.h>
#include <morph/ShapeAnalysis.h>
#include <stdexcept>

/*!
 * Enumerates the way that the guidance molecules are set up
 */
enum class FieldShape
{
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
struct GaussParams
{
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
class RD_Barrel : public morph::RD_Base<Flt>
{
public:

    //! N is the number of thalamo-cortical axon types.
    alignas(Flt) unsigned int N = 5;

    //! M is the number of guidance molecules to use.
    alignas(Flt) unsigned int M = 3;

    //! These are the c_i(x,t) variables from the Karb2004 paper.
    alignas(alignof(std::vector<std::vector<Flt> >))
    std::vector<std::vector<Flt> > c;

    //! To record dci/dt, as this is used in the computation of a.
    alignas(alignof(std::vector<std::vector<Flt> >))
    std::vector<std::vector<Flt> > dc;

    /*!
     * These are the a_i(x,t) variables from the Karb2004 paper. x is a vector in
     * two-space. The first vector is over the different TC axon types, enumerated by
     * i, the second vector are the a_i values, indexed by the vi in the Hexes in
     * HexGrid.
     */
    alignas(alignof(std::vector<std::vector<Flt> >))
    std::vector<std::vector<Flt> > a;

    /*!
     * For each TC axon type, this holds the two components of the gradient field of
     * the scalar value a(x,t) (where this x is a vector in two-space)
     */
    alignas(alignof(std::vector<std::array<std::vector<Flt>, 2> >))
    std::vector<std::array<std::vector<Flt>, 2> > grad_a;

    //! Contains the chemo-attractant modifiers which are applied to a_i(x,t) in Eq 4.
    alignas(alignof(std::vector<std::vector<std::array<std::vector<Flt>, 2> > >))
    std::vector<std::vector<std::array<std::vector<Flt>, 2> > > g;

    //! \hat{a}_i. the sum of the branching densities of all axon types except i.
    alignas(alignof(std::vector<Flt>)) std::vector<Flt> ahat;

    //! Gradient of ahat
    alignas(alignof(std::array<std::vector<Flt>, 2>))
    std::array<std::vector<Flt>, 2> grad_ahat;

    //! divergence of \hat{a}_i(x,t).
    alignas(alignof(std::vector<Flt>)) std::vector<Flt> div_ahat;

    /*!
     * To hold div(g) / 3d, a static scalar field. There are M vectors of N of these
     * vectors of Flts
     */
    alignas(alignof(std::vector<std::vector<std::vector<Flt> > >))
    std::vector<std::vector<std::vector<Flt> > > divg_over3d;

    //! n(x,t) variable from the Karb2004 paper.
    alignas(alignof(std::vector<Flt>))
    std::vector<Flt> n;

    /*!
     * J_i(x,t) variables - the "flux current of axonal branches of type i". This is a
     * vector field.
     */
    alignas(alignof(std::vector<std::array<std::vector<Flt>, 2> >))
    std::vector<std::array<std::vector<Flt>, 2> > J;

    //! Holds the divergence of the J_i(x)s
    alignas(alignof(std::vector<std::vector<Flt> >))
    std::vector<std::vector<Flt> > divJ;

    //! The power to which a_i(x,t) is raised in Eqs 1 and 2 in the paper.
    alignas(Flt) Flt k = 3.0;

protected:
    //! The diffusion parameter.
    alignas(Flt) Flt D = 0.1;

    alignas(Flt) Flt twoDover3dd = this->D+this->D / 3*this->d*this->d;

public:

    //! alpha_i parameters
    alignas(alignof(std::vector<Flt>))
    std::vector<Flt> alpha;
    //! To set a single alpha for all N
    alignas(Flt) Flt alpha_;

    //! beta_i parameters
    alignas(alignof(std::vector<Flt>))
    std::vector<Flt> beta;
    //! To set a single beta for all N
    alignas(Flt) Flt beta_;

    //! Parameters for initial 2D Gaussian masks over the initial branching levels.
    alignas(alignof(std::vector<GaussParams<Flt> >))
    std::vector<GaussParams<Flt> > initmasks;

    /*!
     * The string identifiers for each RT axon. This is of size N. Can be populated
     * from the config file. This allows me to look up the name, as given in the
     * config file, from the floating point index, obtained from (integer index / N)
     */
    alignas(alignof(std::map<Flt, std::string>)) std::map<Flt, std::string> rtnames;

protected: // We have a setter for gamma.
    //! gamma_A/B/C_i (etc) parameters from Eq 4. There are M vectors of Flts in here.
    alignas(alignof(std::vector<std::vector<Flt> >))
    std::vector<std::vector<Flt> > gamma;
    /*!
     * Used for group-based competition. One of the sets of gammas is used as a group
     * identifier.
     */
    alignas(alignof(std::vector<Flt>)) std::vector<Flt> group;
    alignas(alignof(std::set<Flt>)) std::set<Flt> groupset;

public:

    //! Gain for setting gammas, if required.
    Flt G = static_cast<Flt>(0.0);

    /*!
     * A vector of parameters for the direction of the guidance molecules. This is an
     * angle in Radians.
     */
    alignas(alignof(std::vector<Flt>))
    std::vector<Flt> guidance_phi;

    //! Guidance molecule parameters for the width of the function
    alignas(alignof(std::vector<Flt>))
    std::vector<Flt> guidance_width;

    //! Width in orthogonal direction, for 2D fields.
    alignas(alignof(std::vector<Flt>))
    std::vector<Flt> guidance_width_ortho;

    //! Guidance molecule parameters for the offset of the function
    alignas(alignof(std::vector<Flt>))
    std::vector<Flt> guidance_offset;

    //! Guidance molecule parameters to be the gains of the functions
    alignas(alignof(std::vector<Flt>))
    std::vector<Flt> guidance_gain;

    /*!
     * Guidance molecule parameters to be the time (i.e. step) at which each guidance
     * gradient is switched on.
     */
    alignas(alignof(std::vector<unsigned int>))
    std::vector<unsigned int> guidance_time_onset;

    /*!
     * Rho variables in Eq 4 - the concentrations of axon guidance molecules A, B, C,
     * etc. In Karbowski 2004, these are time independent and we will treat them as
     * such, populating them at initialisation.
     *
     * There are M vector<Flts> in rho.
     */
    //@{
    alignas(alignof(std::vector<std::vector<Flt> >))
    std::vector<std::vector<Flt> > rho;
    //@}

    /*!
     * Into grad_rho put the two components of the gradient of rho computed across the
     * HexGrid surface.
     *
     * There are M gradient fields stored in this variable.
     */
    //@{
    alignas(alignof(std::vector<std::array<std::vector<Flt>, 2> >))
    std::vector<std::array<std::vector<Flt>, 2> > grad_rho;
    //@}

    //! Memory to hold an intermediate result
    alignas(alignof(std::vector<std::vector<Flt> >))
    std::vector<std::vector<Flt> > betaterm;

    //! Holds an intermediate value for the computation of Eqs 1 and 2.
    alignas(alignof(std::vector<std::vector<Flt> >))
    std::vector<std::vector<Flt> > alpha_c;

    /*!
     * The contour threshold. For contour plotting [see plot_contour()], the field is
     * normalised, then the contour is plotted where the field crosses this threshold.
     */
    alignas(Flt) Flt contour_threshold = 0.5;

    alignas(Flt) Flt aNoiseGain = 0.1;
    alignas(Flt) Flt aInitialOffset = 0.8;

    /*!
     * Noise for Guidance molecules. Note this is a common parameter, and not a
     * per-guidance molecule parameter.
     */
    alignas(Flt) Flt mNoiseGain = Flt{0};
    alignas(Flt) Flt mNoiseSigma = 0.09; // hex to hex d is usually 0.03


    //! An N element vector holding the sum of a_i for each TC type.
    alignas(std::vector<Flt>) std::vector<Flt> sum_a;

    //! An N element vector holding the initial sum of a_i for each TC type.
    alignas(std::vector<Flt>) std::vector<Flt> sum_a_init;

    //! Data containers for summed n, c and a.
    alignas(std::vector<Flt>) std::vector<Flt> v_nsum;
    alignas(std::vector<Flt>) std::vector<Flt> v_csum;
    alignas(std::vector<Flt>) std::vector<Flt> v_asum;

    //! The fall-off multiplier used towards the boundary of the field
    alignas(std::vector<Flt>) std::vector<Flt> bSig;

    //! epsilon: axon competition parameter
    alignas(Flt) Flt epsilon = 0.2;
    alignas(Flt) Flt epsilonOverNm1 = 0.0;

    //! A per-hex variable. Holds the ID (as a Flt) of the index into c which has the
    //! maximal value of c in each hex.
    alignas(alignof(std::vector<Flt>)) std::vector<Flt> regions;

    //! The centroids of the regions. key is the "ID" of the region - a Flt between 0
    //! and 1, with values separated by 1/N.
    std::map<Flt, std::pair<Flt, Flt> > reg_centroids;
    //! The area of each region, by Flt ID (area in number of hexes).
    std::map<Flt, int> region_areas;
    //! Set true when the spatial analysis has been computed
    bool spatialAnalysisComputed = false;
    /*!
     * ALIGNAS REGION ENDS.
     *
     * Below here, there's no need to worry about alignas keywords.
     */

    //! Sets the function of the guidance molecule method
    std::vector<FieldShape> rhoMethod;

    //! A basic constructor
    RD_Barrel() : morph::RD_Base<Flt>() {}

    /*!
     * Initialise this vector of vectors with noise. This is a model-specific
     * function.
     *
     * I apply a sigmoid to the boundary hexes, so that the noise drops away towards
     * the edge of the domain.
     */
    virtual void noiseify_vector_vector (std::vector<std::vector<Flt> >& vv, std::vector<GaussParams<Flt> >& gp)
    {
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
     * Similar to the above, but just adds noise to v (with a gain only) to \a vv. Has no boundary sigmoid.
     */
    virtual void addnoise_vector (std::vector<Flt>& v)
    {
        std::cout << "Add noise to vector?...";
        if (this->mNoiseGain == Flt{0}) {
            // No noise
            std::cout << "NO.\n";
            return;
        }
        std::cout << "Yes.\n";

        // First, fill a duplicate vector with noise
        morph::RandUniform<Flt> rng(-this->mNoiseGain/Flt{2}, this->mNoiseGain/Flt{2});
        std::vector<Flt> noise (v.size(), Flt{0});
        for (unsigned int h = 0; h<v.size(); ++h) {
            noise[h] = rng.get();
        }

        // Set up the Gaussian convolution kernel on a circular HexGrid.
        morph::HexGrid kernel(this->hextohex_d, Flt{20}*this->mNoiseSigma, 0, morph::HexDomainShape::Boundary);
        kernel.setCircularBoundary (Flt{6}*this->mNoiseSigma);
        std::vector<Flt> kerneldata (kernel.num(), 0.0f);
        // Once-only parts of the calculation of the Gaussian.
        Flt one_over_sigma_root_2_pi = 1 / this->mNoiseSigma * 2.506628275;
        Flt two_sigma_sq = 2.0f * this->mNoiseSigma * this->mNoiseSigma;
        Flt gsum = 0;
        for (auto& k : kernel.hexen) {
            Flt gauss = (one_over_sigma_root_2_pi * std::exp ( -(k.r*k.r) / two_sigma_sq ));
            kerneldata[k.vi] = gauss;
            gsum += gauss;
        }
        // Renormalise
        for (size_t k = 0; k < kernel.num(); ++k) { kerneldata[k] /= gsum; }

        // A vector for the result
        std::vector<Flt> convolved (v.size(), 0.0f);

        // Call the convolution method from HexGrid:
        this->hg->convolve (kernel, kerneldata, noise, convolved);

        // Now add the noise to the vector:
        for (size_t h = 0; h < v.size(); ++h) {
            v[h] += convolved[h];
        }
    }

    /*!
     * Apply a mask to the noise in a vector of vectors. This masks with a 2D Gaussian
     * for each a (there are N TC type, so for each i in N, apply a different Gaussian
     * mask, probably with the same width, but different centre).
     *
     * This allows me to initialise the system in a more biologically realistic manner.
     */
    void mask_a (std::vector<std::vector<Flt> >& vv, std::vector<GaussParams<Flt> >& gp)
    {
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

    //! Perform memory allocations, vector resizes and so on.
    virtual void allocate()
    {
        morph::RD_Base<Flt>::allocate();

        // Resize and zero-initialise the various containers
        this->resize_vector_vector (this->c, this->N);
        this->resize_vector_vector (this->dc, this->N);
        this->resize_vector_vector (this->a, this->N);
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

        this->resize_vector_variable (this->ahat);
        this->resize_gradient_field (this->grad_ahat);
        this->resize_vector_variable (this->div_ahat);

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
    virtual void init()
    {
        this->stepCount = 0;
        this->spatialAnalysisComputed = false;

        // Zero c and n and other temporary variables
        this->zero_vector_vector (this->c, this->N);
        //this->zero_vector_vector (this->a); // gets noisified below
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

        this->zero_vector_variable (this->ahat);
        this->zero_gradient_field (this->grad_ahat);
        this->zero_vector_variable (this->div_ahat);

        // Initialise a with noise
        std::cout << "initmasks.size(): " << this->initmasks.size() << std::endl;
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

        // Having computed gradients, build this->g; has to be done once only. Note
        // that a sigmoid is applied so that g(x) drops to zero around the boundary of
        // the domain.
        this->build_bSig();
        this->build_g();
        this->compute_divg_over3d();

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
    Flt rt_name_to_id (const std::string& idstr)
    {
        Flt theid = -1.0;
        typename std::map<Flt, std::string>::iterator rtn = this->rtnames.begin();
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

    //! Require private setter for d. Slightly different from the base class version.
    void set_d (Flt d_) {
        morph::RD_Base<Flt>::set_d (d_);
        this->updateTwoDover3dd();
    }

public:
    /*!
     * Public accessors for D, as it requires another attribute to be updated at the
     * same time.
     */
    void set_D (Flt D_)
    {
        this->D = D_;
        this->updateTwoDover3dd();
    }
    Flt get_D()
    {
        return this->D;
    }

protected:
    //! Compute 2D/3d^2 (and 1/3d^2 too)
    void updateTwoDover3dd() {
        this->twoDover3dd = (this->D+this->D) / (3*this->d*this->d);
    }

public:
    /*!
     * setGamma for the guidance molecule index m_idx and the RT index n_idx to
     * @value. If group_m==m_idx, then set this->group[n_idx]=@value
     */
    int setGamma (unsigned int m_idx, unsigned int n_idx, Flt value, unsigned int group_m = 0)
    {
        if (gamma.size() > m_idx) {
            if (gamma[m_idx].size() > n_idx) {
                // Ok, we can set the value
                //std::cout << "Setting gamma[m="<<m_idx<<"][n="<<n_idx<<"] to " << value << std::endl;
                //std::cout << m_idx << "," << n_idx << "," << value << std::endl;
                this->gamma[m_idx][n_idx] = value;
                if (group_m == m_idx) {
                    this->group[n_idx] = value;
                    this->groupset.insert (value);
                }
            } else {
                std::cerr << "WARNING: DID NOT SET GAMMA (too few RT axon types for n_idx=" << n_idx << ")" << std::endl;
                return 1;
            }
        } else {
            std::cerr << "WARNING: DID NOT SET GAMMA (too few guidance molecules for m_idx=" << m_idx << ")" << std::endl;
            return 2;
        }
        return 0;
    }

    //! Save the c, a and n variables.
    virtual void save()
    {
        std::stringstream fname;
        fname << this->logpath << "/c_";
        fname.width(5);
        fname.fill('0');
        fname << this->stepCount << ".h5";
        morph::HdfData data(fname.str());
        for (unsigned int i = 0; i<this->N; ++i) {
            std::stringstream path;
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

    //! Save the hexgrid information into a separate file
    void saveHG()
    {
        std::stringstream hgname;
        hgname << this->logpath << "/hexgrid.h5";
        this->hg->save(hgname.str().c_str());
    }

#if 1
    // Save out spatial analysis results
    void saveRegions()
    {
        std::stringstream fname;
        fname << this->logpath << "/regions_";
        fname.width(5);
        fname.fill('0');
        fname << this->stepCount << ".h5";
        morph::HdfData data(fname.str());
        // Save the region centroids
        std::vector<Flt> keys;
        std::vector<Flt> x_;
        std::vector<Flt> y_;
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

    //! Save asum, nsum and csum. Call once at end of simulation.
    void savesums()
    {
        std::stringstream fname;
        fname << this->logpath << "/sums.h5";
        morph::HdfData data(fname.str());
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
    void saveGuidance()
    {
        std::stringstream fname;
        fname << this->logpath << "/guidance.h5";
        morph::HdfData data(fname.str());
        for (unsigned int m = 0; m<this->M; ++m) {
            std::stringstream path;
            path << "/rh" << m;
            std::string pth(path.str());
            data.add_contained_vals (pth.c_str(), this->rho[m]);
            pth[1] = 'g'; pth[2] = 'x';
            data.add_contained_vals (pth.c_str(), this->grad_rho[m][0]);
            pth[2] = 'y';
            data.add_contained_vals (pth.c_str(), this->grad_rho[m][1]);
            for (unsigned int i = 0; i<this->N; ++i) {
                std::stringstream path;
                path << "/divg_" << m << "_" << i;
                std::string pth(path.str());
                data.add_contained_vals (pth.c_str(), this->divg_over3d[m][i]);
            }
            std::stringstream gpath;
            gpath << "/gamma" << m;
            data.add_contained_vals (gpath.str().c_str(), this->gamma[m]);
        }
    }

    //! Compute the values of c, the connection density
    virtual void integrate_c()
    {
        // 3. Do integration of c
        for (unsigned int i=0; i<this->N; ++i) {

#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; h++) {
                // Note: betaterm used in compute_dci_dt()
                this->betaterm[i][h] = this->beta[i] * this->n[h] * static_cast<Flt>(pow (this->a[i][h], this->k));
            }

            // Runge-Kutta integration for C (or ci)
            std::vector<Flt> qq(this->nhex,0.);
            std::vector<Flt> k1 = this->compute_dci_dt (this->c[i], i);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; h++) {
                qq[h] = this->c[i][h] + k1[h] * this->halfdt;
            }

            std::vector<Flt> k2 = this->compute_dci_dt (qq, i);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; h++) {
                qq[h] = this->c[i][h] + k2[h] * this->halfdt;
            }

            std::vector<Flt> k3 = this->compute_dci_dt (qq, i);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; h++) {
                qq[h] = this->c[i][h] + k3[h] * this->dt;
            }

            std::vector<Flt> k4 = this->compute_dci_dt (qq, i);
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
     * A possibly normalization-function specific task to carry out once after the sum
     * of a has been computed.
     */
    virtual void sum_a_computation (const unsigned int _i)
    {
        // Compute the sum of a[i] across the sheet.
        this->sum_a[_i] = 0.0;
        Flt sum_tmp = 0.0;
#pragma omp parallel for reduction(+:sum_tmp)
        for (unsigned int h=0; h<this->nhex; ++h) {
            sum_tmp += this->a[_i][h];
        }
        this->sum_a[_i] = sum_tmp;
    }

    //! The normalization/transfer function.
    virtual inline Flt transfer_a (const Flt& _a, const unsigned int _i)
    {
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

    //! Compute the values of a, the branching density
    virtual void integrate_a()
    {
        // 2. Do integration of a (RK in the 1D model). Involves computing axon
        // branching flux.

        // Pre-compute:
        // 1) The intermediate val alpha_c.
        for (unsigned int i=0; i<this->N; ++i) {
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                this->alpha_c[i][h] = this->alpha[i] * this->c[i][h];
            }
        }

        // Runge-Kutta: No OMP here - there are only N(<10) loops, which isn't enough
        // to load the threads up.
        for (unsigned int i=0; i<this->N; ++i) {

            // --- START specific to comp2 method ---
            // Compute "the sum of all a_j for which j!=i"
            this->zero_vector_variable (this->ahat);
            for (unsigned int j=0; j<this->N; ++j) {
                if (j==i) { continue; }
#pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; ++h) {
                    this->ahat[h] += this->a[j][h];
                }
            }
            // 1.1 Compute divergence and gradient of ahat
            this->compute_divahat();
            this->spacegrad2D (this->ahat, this->grad_ahat);
            // --- END specific to comp2 method ---

            // Runge-Kutta integration for A
            std::vector<Flt> qq(this->nhex, 0.0);
            this->compute_divJ (this->a[i], i); // populates divJ[i]

            std::vector<Flt> k1(this->nhex, 0.0);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                k1[h] = this->divJ[i][h] - this->dc[i][h];
                qq[h] = this->a[i][h] + k1[h] * this->halfdt;
            }

            std::vector<Flt> k2(this->nhex, 0.0);
            this->compute_divJ (qq, i);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                k2[h] = this->divJ[i][h] - this->dc[i][h];
                qq[h] = this->a[i][h] + k2[h] * this->halfdt;
            }

            std::vector<Flt> k3(this->nhex, 0.0);
            this->compute_divJ (qq, i);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                k3[h] = this->divJ[i][h] - this->dc[i][h];
                qq[h] = this->a[i][h] + k3[h] * this->dt;
            }

            std::vector<Flt> k4(this->nhex, 0.0);
            this->compute_divJ (qq, i);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                k4[h] = this->divJ[i][h] - this->dc[i][h];
                this->a[i][h] += (k1[h] + 2.0 * (k2[h] + k3[h]) + k4[h]) * this->sixthdt;
            }
        }
    }

    virtual void compute_n()
    {
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
    virtual void summation_a()
    {
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
    virtual void step()
    {
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
    void debug_values (std::vector<Flt>& f, Flt dangerThresh)
    {
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
    std::vector<Flt> compute_dci_dt (std::vector<Flt>& f, unsigned int i)
    {
        std::vector<Flt> dci_dt (this->nhex, 0.0);
#pragma omp parallel for
        for (unsigned int h=0; h<this->nhex; h++) {
            dci_dt[h] = (this->betaterm[i][h] - this->alpha[i] * f[h]);
        }
        return dci_dt;
    }

    //! Compute bSig vector, which need be carried out once only.
    void build_bSig()
    {
        for (auto h : this->hg->hexen) {
            // Sigmoid/logistic fn params: 100 sharpness, 0.02 dist offset from boundary
            this->bSig[h.vi] = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary-this->boundaryFalloffDist)) );
        }
    }

    //! Build g from the gradient of rho and the gammas.
    virtual void build_g()
    {
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
    void compute_divg_over3d()
    {
        // Change to have one for each m in M? They should then sum, right?

        for (unsigned int i = 0; i < this->N; ++i) {

#pragma omp parallel for schedule(static)
            for (unsigned int hi=0; hi<this->nhex; ++hi) {

                std::vector<Flt> divg(this->M, 0.0);
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
     * Inputs: this->g, fa (which is this->a[i] or a q in the RK algorithm), this->D, \a a
     * i, the TC type.  Helper functions: spacegrad2D().  Output: this->divJ
     *
     * In the competition term, it's possible to set \bar{a} equal to either sigmoid
     * transfer function of {a} with SIGMOID_ROLLOFF_FOR_A (which maxes out at 2.0) or a
     * linear transfer function of {a} with a maximum of 2.0 with LINEAR_MAX. Do so when
     * compiling with, e.g. -DLINEAR_MAX. Initially, I thought a transfer function was
     * necessary, but it is not (though use of a transfer function does extend the range
     * of parameters for which this model is stable).
     */
    virtual void compute_divJ (std::vector<Flt>& fa, unsigned int i)
    {
        // Compute gradient of a_i(x), for use computing the third term, below.
        this->spacegrad2D (fa, this->grad_a[i]);

        if (this->N > 0) {
            this->epsilonOverNm1 = this->epsilon/(this->N-1);
        } else {
            this->epsilonOverNm1 = 0.0;
        }

        // _Five_ terms to compute; see Eq. 17 in methods_notes.pdf. Copy comp3.
        volatile bool breakflag = false;
#pragma omp parallel for shared(breakflag) //schedule(static) // This was about 10% faster than schedule(dynamic,50).
        for (unsigned int hi=0; hi<this->nhex; ++hi) {

            // In OpenMP, iterations must be finished.
            if (breakflag) { continue; }

            // 1. The D Del^2 a_i term. Eq. 18.
            // 1a. Or D Del^2 Sum(a_i) (new)
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
            if (isnan(term1)) {
                std::cerr << "term1 isnan" << std::endl;
                std::cerr << "thesum is " << thesum << " fa[hi=" << hi << "] = " << fa[hi] << std::endl;
                breakflag = true;
            }

            Flt term1_1 = Flt{0};

            // Term 1.1 is F/N-1 a div(ahat)
            term1_1 = this->epsilonOverNm1 * fa[hi] * this->div_ahat[hi];

            if (isnan(term1_1)) {
                std::cerr << "term1_1 isnan" << std::endl;
                std::cerr << "fa[hi="<<hi<<"] = " << fa[hi] << ", this->div_ahat[hi] = " << this->div_ahat[hi] << std::endl;
                breakflag = true;
            }

            Flt term1_2 = Flt{0};
            // Term 1.2 is eps/N-1 grad(ahat) . grad(a)
            term1_2 = this->epsilonOverNm1 * (this->grad_ahat[0][hi] * this->grad_a[i][0][hi]
                                              + this->grad_ahat[1][hi] * this->grad_a[i][1][hi]);

            if (isnan(term1_2)) {
                std::cerr << "term1_2 isnan at hi=" << hi << std::endl;
                if (isnan(this->grad_ahat[0][hi])) {
                    std::cerr << "grad_ahat[0][hi] isnan\n";
                    breakflag = true;
                }
                if (isnan(this->grad_ahat[1][hi])) {
                    std::cerr << "grad_ahat[1][hi] isnan\n";
                    breakflag = true;
                }
            }

            // 2. The (a div(g)) term.
            Flt term2 = 0.0;

            // 3. Third term is this->g . grad a_i. Should not contribute to J, as
            // g(x) decays towards boundary.
            Flt term3 = 0.0;

            for (unsigned int m =0 ; m < this->M; ++m) {
                if (this->stepCount >= this->guidance_time_onset[m]) {
                    // g contributes to term2
                    term2 += fa[hi] * this->divg_over3d[m][i][hi];
                    // and to term3
                    term3 += this->g[m][i][0][hi] * this->grad_a[i][0][hi] + (this->g[m][i][1][hi] * this->grad_a[i][1][hi]);
                }
            }

            // - term1_1/2 or + term1_1/2? It's + in the supp.tex
            this->divJ[i][hi] = term1 + term1_1 + term1_2 - term2 - term3;
        }
        if (breakflag == true) {
            throw std::runtime_error ("compute_divJ: There was a nan");
        }
    }

    //! Compute divergence of \hat{a}_i
    void compute_divahat()
    {
        volatile bool breakflag = false;
#pragma omp parallel for shared(breakflag)
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            if (breakflag) { continue; }
            Flt thesum = -6 * this->ahat[hi];
            thesum += this->ahat[(HAS_NE(hi)  ? NE(hi)  : hi)];
            thesum += this->ahat[(HAS_NNE(hi) ? NNE(hi) : hi)];
            thesum += this->ahat[(HAS_NNW(hi) ? NNW(hi) : hi)];
            thesum += this->ahat[(HAS_NW(hi)  ? NW(hi)  : hi)];
            thesum += this->ahat[(HAS_NSW(hi) ? NSW(hi) : hi)];
            thesum += this->ahat[(HAS_NSE(hi) ? NSE(hi) : hi)];
            this->div_ahat[hi] = this->twoover3dd * thesum;
            if (isnan(this->div_ahat[hi])) {
                std::cerr << "div ahat isnan" << std::endl;
                breakflag = true;
            }
        }
        if (breakflag == true) {
            throw std::runtime_error ("div ahat isnan");
        }
    }

    /*!
     * Generate Gaussian profiles for the chemo-attractants.
     *
     * Instead of using the Karbowski equations, just make some gaussian 'waves'
     *
     * \a m The molecule id
     */
    void gaussian1D_guidance (unsigned int m)
    {
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
     * \a m The molecule id
     */
    void gaussian2D_guidance (unsigned int m)
    {

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
     * \a m The molecule id
     */
    void exponential_guidance (unsigned int m)
    {
        for (auto h : this->hg->hexen) {
            Flt cosphi = (Flt) cos (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt sinphi = (Flt) sin (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt x_ = (h.x * cosphi) + (h.y * sinphi);
            this->rho[m][h.vi] = exp (this->guidance_gain[m] * (x_-guidance_offset[m]));
        }
    }

    //! \a m The molecule id
    void sigmoid_guidance (unsigned int m)
    {
        for (auto h : this->hg->hexen) {
            Flt cosphi = (Flt) cos (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt sinphi = (Flt) sin (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            DBG2 ("phi= " << this->guidance_phi[m] << ". cosphi: " << cosphi << " sinphi: " << sinphi);
            Flt x_ = (h.x * cosphi) + (h.y * sinphi);
            DBG2 ("x_[" << h.vi << "] = " << x_);
            this->rho[m][h.vi] = guidance_gain[m] / (1.0 + exp(-(x_-guidance_offset[m])/this->guidance_width[m]));
        }
    }

    //! \a m The molecule id
    void linear_guidance (unsigned int m)
    {
        for (auto h : this->hg->hexen) {
            Flt cosphi = (Flt) cos (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt sinphi = (Flt) sin (this->TWOPI_OVER_360 * this->guidance_phi[m]);
            Flt x_ = (h.x * cosphi) + (h.y * sinphi);
            this->rho[m][h.vi] = (x_-guidance_offset[m]) * this->guidance_gain[m];
        }
    }

    //! \a m The molecule id
    void circlinear_guidance (unsigned int m)
    {
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
    //! Carry out any sensible spatial analysis required
    virtual void spatialAnalysis()
    {
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
}; // RD_Barrel
