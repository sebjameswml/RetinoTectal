/*
 * Retino-tectal non-competition system. This specialization inits all the RT
 * projections using RetArrange. RD_RetNoComp is the base class, which implements the
 * no-competition base model.
 */

#include "rd_retnocomp.h"
#include "retarrange.h"
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
#include <iostream>
using std::cout;
using std::endl;

#define scf(a) static_cast<Flt>(a)

template <class Flt>
class RD_RetTec_NoComp : public RD_RetNoComp<Flt>, public RetArrange<Flt>
{
public:
    RD_RetTec_NoComp (void) {}
    ~RD_RetTec_NoComp (void) {}

    virtual void allocate (void) {
        RD_RetNoComp<Flt>::allocate();
    }

    virtual void init (void) {
        cout << "RD_RetTec_NoComp init()" << endl;
        // Set up alpha, beta, epsilon (before RD_RetNoComp::init)
        this->setupPerNParams();
        // Because we derive from RetArrange, we have tec_coords and ret_coords as members.
        RetArrange<Flt>::initRet (this->N, this->ellipse_a, this->ellipse_b);
        // Now, from ret_coords, populate gammas, assuming two, orthogonal morphogen
        // fields which are noisy.
        this->setupGammas();
        // the rest of the init chain has to be called AFTER gammas are set up
        RD_RetNoComp<Flt>::init();
    }

protected:
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
    }

    void setupPerNParams (void) {
        // Index through thalamocortical fields, setting params:
        for (unsigned int i = 0; i < this->N; ++i) {
            this->alpha[i] = this->alpha_;
            this->beta[i] = this->beta_;
            // Note: No epsilon setup here
            // Sets up mask for initial branching density. This is hardcoded here.
            GaussParams<Flt> gp;
            gp.gain = 1.0;
            gp.sigma = 0.0;
            gp.x = 0.0;
            gp.y = 0.0;
            this->initmasks.push_back (gp);
        }
    }

public:
    //! Override step() to recompute g on each sim step (for noise)
    virtual void step (void) {
        this->stepCount++;
        // 0. Rebuild g, with a new set of random samples, if necessary
        if (this->sigma_rho != scf(0.0)) {
            this->build_g();
            this->compute_divg_over3d();
        }
        this->compute_n();
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

public:
    //! spatial analysis, adding the tec_coords-reg_centroids computation
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

}; // RD_RetTecMI
