/*
 * Retinotectal model resembling one presented by Hugh Simpson and Geoffrey
 * Goodhill in "A simple model can unify a broad range of phenomena in retinotectal map
 * development", Biol Cybern (2011) 104:9-29
 *
 * I'm bringing the idea of variable interaction with signalling gradients and
 * competition to try to get rid of the non-biological part of Simpson & Goodhill's
 * work.
 */

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <set>
#include <chrono>

#include <morph/Config.h>
#include <morph/Random.h>
#include <morph/HdfData.h>

#ifdef VISUALISE
# include <morph/Visual.h>
# include <morph/GraphVisual.h>
#endif

#include "branch.h"
#include "branch_geb.h"
#include "net.h"
#include "tissue.h"

#ifdef VISUALISE
# include "branchvisual.h"
# include "netvisual.h"
# include "tissuevisual.h"
#endif

template<typename T, size_t N, typename B=branch<T, N>>
struct Agent1
{
    Agent1 (morph::Config* cfg, morph::Config* mcfg)
    {
        this->conf = cfg;
        this->mconf = mcfg;
        this->init();
    }
    ~Agent1()
    {
        delete this->ret;
        delete this->tectum;
    }

    static constexpr unsigned int showevery = 100;
    static constexpr unsigned int visevery = 5;

    //! Just show the tissue. Don't use at same time as run()
    void showtissue()
    {
#ifdef VISUALISE
        this->tvisinit();
        this->tvv->render();
        this->tvv->keepOpen();
#else
        throw std::runtime_error ("showtissue() only works if VISUALISE was defined during compilation.");
#endif
    }

    //! Run this model!
    void run()
    {
#ifdef VISUALISE
        this->visinit();
#endif
        std::chrono::steady_clock::time_point laststep = std::chrono::steady_clock::now();

        typename std::vector<B>::iterator pending_br_it = this->pending_branches.begin();
        typename morph::vVector<size_t>::iterator pb_sz_it = this->pb_sizes.begin();
        // How often to introduce groups of axons? 0 means 'all at once'
        unsigned int intro_every = this->mconf->getUInt ("intro_every", 0);
        if (intro_every == 0) {
            this->branches.resize (this->pending_branches.size());
            this->branches.swap (this->pending_branches);
        }

        // Are we running the random 'model'?
        if (this->mconf->getString ("model", "axgrad") == "random") {
            std::cout << "branches.size(): " << this->branches.size() << std::endl;
            this->steprandom();
            std::cout << "RMS error of axon centroids: " << this->ax_centroids.rms() << std::endl;
#ifdef VISUALISE
            this->vis(1);
            this->v->keepOpen();
#endif
            return;
        }

        for (unsigned int i = 0; i < this->conf->getUInt ("steps", 1000); ++i) {

            if (intro_every > 0 && pb_sz_it != this->pb_sizes.end() && i%intro_every == 0) {
                // Introduce some of pending_branches into branches
                typename std::vector<B>::iterator brit = this->branches.end();
                std::cout << "Adding " << (*pb_sz_it) << " to this->branches...\n";
                this->branches.resize (this->branches.size() + *pb_sz_it);
                brit = this->branches.end();
                brit -= *pb_sz_it;
                std::copy (pending_br_it, pending_br_it+*pb_sz_it, brit);
                pending_br_it += *pb_sz_it;
                pb_sz_it++;
            }

            this->step();
#ifdef VISUALISE
            if (i%visevery == 0) { this->vis(i); } // Visualize every 10 doubles time required for program
#endif
            if (i%showevery == 0) {
                std::chrono::steady_clock::duration since = std::chrono::steady_clock::now() - laststep;
                std::cout << "step " << i << ". Per step: "
                          << std::chrono::duration_cast<std::chrono::milliseconds>(since).count()/showevery << " ms\n";
                std::cout << "RMS error of axon centroids: " << this->ax_centroids.rms() << std::endl;
                laststep = std::chrono::steady_clock::now();
            }
        }
        std::cout << "Done simulating\n";
#ifdef VISUALISE
        this->vis(this->conf->getUInt ("steps", 1000));
        this->v->keepOpen();
#endif
    }

    //! Save any relevant results of the simulation to an HdfData object.
    void save (const std::string& outfile)
    {
        morph::HdfData d(outfile, morph::FileAccess::TruncateWrite);
        d.add_val ("/sos", this->ax_centroids.sos());
        d.add_val ("/rms", this->ax_centroids.rms());
        d.add_val ("/steps", this->conf->getUInt ("steps", 1000));
    }

#ifdef VISUALISE
    //! Update the visualization
    unsigned int framenum = 0;
    void vis (unsigned int stepnum)
    {
        // crash with e_single occurs at either of these glfw calls
        if (this->goslow == true) {
            glfwWaitEventsTimeout (0.1); // to add artificial slowing
        } else {
            glfwPollEvents();
        }
        this->bv->reinit();
        this->av->reinit();
        this->cv->reinit();
        this->tcv->reinit();
        this->gv->append ((float)stepnum, this->ax_centroids.sos(), 0);
        this->v->render();
        //this->tvv->render();
        if (this->conf->getBool ("movie", false)) {
            std::stringstream frame;
            frame << "log/agent/";
            frame.width(4);
            frame.fill('0');
            frame << framenum++;
            frame << ".png";
            this->v->saveImage(frame.str());
        }
    }
#endif

    //! One step of the simulation in which branches postitions are randomly set
    void steprandom()
    {
        morph::RandUniform<T, std::mt19937> rng(T{0}, T{1.0});
#if 1
        // Randomise axon growth cones uniformly.
        for (auto& b : this->branches) {
            b.next = {rng.get(), rng.get()};
            b.current = b.next;
        }
        // Update centroids, which should end up as a Poisson distribution
        for (unsigned int i = 0; i < this->branches.size()/this->bpa; ++i) { this->ax_centroids.p[i] = {T{0}, T{0}, T{0}}; }
        for (auto& b : this->branches) {
            this->ax_centroids.p[b.aid][0] += b.current[0] / static_cast<T>(this->bpa);
            this->ax_centroids.p[b.aid][1] += b.current[1] / static_cast<T>(this->bpa);
        }
#else
        // Only update centroids of ax_centroids
        for (unsigned int i = 0; i < this->ax_centroids.p.size(); ++i) {
            this->ax_centroids.p[i] = {rng.get(), rng.get(), rng.get()};
        }
#endif
   }

    //! Perform one step of the simulation
    void step()
    {
        // Compute the next position for each branch:
#ifdef __OSX__
        // Mac compiler didn't like omp parallel for in front of a for(auto...
#pragma omp parallel for
        for (unsigned int i = 0; i < this->branches.size(); ++i) {
            this->branches[i].compute_next (this->branches, this->tectum, this->m);
        }
#else
#pragma omp parallel for
        for (auto& b : this->branches) { b.compute_next (this->branches, this->tectum, this->m); }
#endif
        // Update centroids
        for (unsigned int i = 0; i < this->branches.size()/this->bpa; ++i) { this->ax_centroids.p[i] = {T{0}, T{0}, T{0}}; }
        for (auto& b : this->branches) {
            this->ax_centroids.p[b.aid][0] += b.next[0] / static_cast<T>(this->bpa);
            this->ax_centroids.p[b.aid][1] += b.next[1] / static_cast<T>(this->bpa);
        }
        // Record centroid history for the axons to be seen
        for (size_t ax = 0; ax < this->ax_centroids.p.size(); ++ax) {
            if (this->seeaxons.count(ax)) { this->ax_history[ax].push_back (ax_centroids.p[ax]); }
        }

        // Once 'next' has been updated for each branch, run through & copy it to 'current'
        for (auto& b : this->branches) { b.current = b.next; }
    }

#ifdef VISUALISE
    //! Create a tissue visual, to reduce boilerplate code in init()
    tissuevisual<float, N>* createTissueVisual (morph::Vector<T,3>& offset, guidingtissue<T, N>* gtissue,
                                                const std::string& tag,
                                                expression_view exview, size_t pair_to_view, bool alt_cmap=false)
    {
        tissuevisual<float, N>* tv = new tissuevisual<float, N>(this->tvv->shaderprog, this->tvv->tshaderprog, gtissue, offset);
        tv->view = exview;
        tv->pair_to_view = pair_to_view;
        tv->cm.setType (morph::ColourMapType::Duochrome);
        if (alt_cmap == true) {
            tv->cm.setHueCM();
        } else {
            tv->cm.setHueRG();
        }
        std::stringstream ss;
        if (exview == expression_view::receptor_exp) {
            ss << tag << " receptor expression " << (pair_to_view*2) << "/" << (1+pair_to_view*2);
        } else if (exview == expression_view::receptor_grad_x) {
            ss << tag << " receptor gradient " << (pair_to_view*2) << "/" << (1+pair_to_view*2) << " x";
        } else if (exview == expression_view::receptor_grad_y) {
            ss << tag << " receptor gradient " << (pair_to_view*2) << "/" << (1+pair_to_view*2) << " y";
        } else if (exview == expression_view::ligand_exp) {
            ss << tag << " ligand expression " << (pair_to_view*2) << "/" << (1+pair_to_view*2);
        } else if (exview == expression_view::ligand_grad_x) {
            ss << tag << " ligand gradient " << (pair_to_view*2) << "/" << (1+pair_to_view*2) << " x";
        } else if (exview == expression_view::ligand_grad_y) {
            ss << tag << " ligand gradient " << (pair_to_view*2) << "/" << (1+pair_to_view*2) << " y";
        }
        tv->addLabel (ss.str(), {0.0f, 1.1f, 0.0f});
        tv->finalize();

        return tv;
    }
#endif

    static constexpr size_t singleaxon_idx = 210;
    // Set true to use orthographic projection
    static constexpr bool use_ortho = false;

    //! Simulation init
    void init()
    {
        /*
         * Create tectum and retina tissue objects
         */

        this->rgcside = this->mconf->getUInt ("rgcside", this->rgcside);
        T gr_denom = rgcside-1;
        T gr = T{1}/gr_denom; // gr is grid element length

        this->tectum = new guidingtissue<T, N>(this->rgcside, this->rgcside, {gr, gr}, {0.0f, 0.0f},
                                               (expression_form)this->mconf->getUInt ("tectum_form", 3),
                                               (expression_form)this->mconf->getUInt ("tectum_form", 3));

        this->ret = new guidingtissue<T, N>(this->rgcside, this->rgcside, {gr, gr}, {0.0f, 0.0f},
                                            (expression_form)this->mconf->getUInt ("retina_form", 2),
                                            (expression_form)this->mconf->getUInt ("retina_form", 2));

        /*
         * Apply any manipulations to retina or tectum
         */

        // Use a variable to prevent multiple manipulations from being applied. May
        // need to tweak this to allow selected combinations of manipulations.
        bool manipulated = false;

        // Graft swap manipulation
        Json::Value gs_coords = this->conf->getValue ("graftswap_coords");
        //std::cout << gs_coords << std::endl;
        Json::Value l1 = gs_coords.get ("locn1", "[0,0]");
        Json::Value l2 = gs_coords.get ("locn2", "[0,0]");
        Json::Value ps = gs_coords.get ("patchsize", "[0,0]");
        morph::Vector<size_t, 2> l1v;
        morph::Vector<size_t, 2> l2v;
        morph::Vector<size_t, 2> psv;
        // now, those are only used for rotation and graftswap, so only error on that combination
        if ((this->conf->getBool ("tectal_graftswap", false)
             || this->conf->getBool ("retinal_graftswap", false)
             || (this->conf->getUInt ("retinal_rotations", 0) > 0)
             || (this->conf->getUInt ("tectal_rotations", 0) > 0) )
            && (l1.size() < 2 || l2.size() < 2 || ps.size() < 2)) {
            throw std::runtime_error ("Bad values for locn1, locn2 or patchsize.");
        }

        if (this->conf->getBool ("tectal_graftswap", false)
            || this->conf->getBool ("retinal_graftswap", false)
            || (this->conf->getUInt ("retinal_rotations", 0) > 0)
            || (this->conf->getUInt ("tectal_rotations", 0) > 0) ) {
            l1v = { l1[0].asUInt(), l1[1].asUInt() };
            l2v = { l2[0].asUInt(), l2[1].asUInt() };
            psv = { ps[0].asUInt(), ps[1].asUInt() };
            std::cout << "Location 1: " << l1v << " Location 2: " << l2v << " patch size: " << psv << std::endl;
        }

        if (this->conf->getBool ("tectal_graftswap", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->tectum->graftswap (l1v, psv, l2v);
            manipulated = true;
        }
        if (this->conf->getBool ("retinal_graftswap", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->ret->graftswap (l1v, psv, l2v);
            manipulated = true;
        }

        // Graft and rotate manipulation
        if (l1.size() > 1 && ps.size() > 1) {
            l1v = { l1[0].asUInt(), l1[1].asUInt() };
            psv = { ps[0].asUInt(), ps[1].asUInt() };
            if (this->conf->getUInt ("retinal_rotations", 0) > 0) {
                if (psv[0] != psv[1]) {
                    throw std::runtime_error ("patch size has to be square for rotations");
                }
                this->ret->graftrotate (l1v, psv[0], (size_t)this->conf->getUInt ("retinal_rotations", 0));
                manipulated = true;
            }
            if (this->conf->getUInt ("tectal_rotations", 0) > 0) {
                if (psv[0] != psv[1]) {
                    throw std::runtime_error ("patch size has to be square for rotations");
                }
                this->tectum->graftrotate (l1v, psv[0], (size_t)this->conf->getUInt ("tectal_rotations", 0));
                manipulated = true;
            }
        }

        // Compound retina manipulation
        if (this->conf->getBool ("compound_retina", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->ret->compound_tissue(tissue_region::left_half);
            manipulated = true;
        }

        // Various tissue ablation manipulations are possible, we only ablate ret left or tec top.
        if (this->conf->getBool ("ablate_ret_left", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->ret->ablate_left_half();
            manipulated = true;
        }
        if (this->conf->getBool ("ablate_tec_top", false)) {
            // Note: we ARE allowed to have ablate_tec_top along with ablate_ret_left (the 'mismatch' manipulation)
            if (manipulated && !this->conf->getBool ("ablate_ret_left", false)) {
                throw std::runtime_error ("Code is only tested for one manipulation at a time!");
            }
            this->tectum->ablate_top_half();
            manipulated = true;
        }

        // Knockin will increase expression for half of all cells. There's code in tissue.h to make this randomised or regular
        T affected = T{0.5};
        T amount = T{1};
        if (this->conf->getBool ("reber", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            // Knockin/knockdown receptor 0:
            this->ret->receptor_knockin (0, affected, amount);
            this->ret->receptor_knockdown (0, amount);
            manipulated = true;
        }

        std::cout << "Retina has " << this->ret->num() << " cells\n";
        std::cout << "Tectum has " << this->tectum->num() << " cells\n";
        std::cout << "Retina is " << this->ret->w << " wide and " << this->ret->h << " high\n";

        /*
         * Set up the axon branches (Agent1::pending_branches, which will be copied into
         * Agent1::branches) and the ax_centroids net object used for visualisation/analysis.
         */

        this->bpa = this->mconf->getUInt ("bpa", 8);

        size_t num_branches = this->ret->num() * this->bpa;

        // Can override num branches for single/sparse experiments
        if (conf->getBool ("singleaxon", false)) { num_branches = this->bpa; }

        this->pending_branches.resize(num_branches);

        if (conf->getBool ("singleaxon", false)) {
            this->ax_centroids.init (1, 1);
        } else {
            // Although the axons arrange themselves on the tectum, use retina w/h to
            // initialise as it is the NUMBER of axons on the retina which counts.
            this->ax_centroids.init (this->ret->w, this->ret->h);
            // To help visualization
            this->ax_centroids.domain_w = this->tectum->w;
            this->ax_centroids.domain_h = this->tectum->h;
            this->ax_centroids.dx = this->tectum->dx;
        }

        // Axon initial positions x and y can be uniformly randomly selected...
        morph::RandUniform<T, std::mt19937> rng_x(T{0}, T{1.0});
        morph::RandUniform<T, std::mt19937> rng_y(T{-0.2}, T{0});
        // ...or set from the ideal position plus a random perturbation
        morph::RandNormal<T, std::mt19937> rng_p0(T{0}, T{0.1});
        // A normally distributed perturbation is added for each branch. SD=0.1.
        morph::RandNormal<T, std::mt19937> rng_p(T{0}, T{0.1});
        // Generate random number sequences all at once
        size_t axc_sz = this->ax_centroids.p.size();
        std::vector<T> rn_x = rng_x.get (axc_sz); // ax_centroids size?
        std::vector<T> rn_y = rng_y.get (axc_sz);
        std::vector<T> rn_p = rng_p.get (axc_sz * 2 * this->bpa);
        std::vector<T> rn_p0 = rng_p0.get (axc_sz * 2 * this->bpa);
        bool totally_random = this->mconf->getBool ("totally_random_init", true);
        std::string branch_model = this->mconf->getString ("branch_model", "james_agent");

        // A loop to set up each branch object in pending_branches.
        for (unsigned int i = 0; i < this->pending_branches.size(); ++i) {
            // Set the branch's termination zone
            unsigned int ri = i/bpa; // retina index
            this->pending_branches[i].aid = (int)ri; // axon index
            if (conf->getBool ("singleaxon", false)) {
                this->pending_branches[i].rcpt = this->ret->rcpt[singleaxon_idx]; // FIXME: Use seeaxons
                this->pending_branches[i].lgnd = this->ret->lgnd[singleaxon_idx];
                this->pending_branches[i].target = this->ret->posn[singleaxon_idx];
            } else {
                this->pending_branches[i].rcpt = this->ret->rcpt[ri];
                this->pending_branches[i].lgnd = this->ret->lgnd[ri];
                this->pending_branches[i].target = this->ret->posn[ri];
            }
            // Call the first interaction parameter 'EphA'
            rcpt_max =  this->pending_branches[i].rcpt[0] > rcpt_max ? pending_branches[i].rcpt[0] : rcpt_max;
            rcpt_min =  this->pending_branches[i].rcpt[0] < rcpt_min ? pending_branches[i].rcpt[0] : rcpt_min;

            // Set as in the S&G paper - starting at bottom in region x=(0,tectum->w), y=(-0.2,0)
            morph::Vector<T, 3> initpos;
            if (totally_random == true) {
                if (branch_model == "gebhardt") {
                    // In Gebhardt model, arrange randomly along x axis
                    initpos = { rn_x[ri] + rn_p[2*i], 0, 0 };
                    initpos[0] = initpos[0] > 1 ? 1 : initpos[0];
                    initpos[0] = initpos[0] < 0 ? 0 : initpos[0];
                } else {
                    initpos = { rn_x[ri] + rn_p[2*i], rn_y[ri] + rn_p[2*i+1], 0 };
                }
            } else {
                morph::Vector<T, 2> init_offset = { T{0}, T{-0.5} };
                morph::Vector<T, 2> init_mult = { T{1}, T{0.2} };
                initpos.set_from ((this->pending_branches[i].target*init_mult) + init_offset);
                initpos[0] += rn_p0[2*i];
                initpos[1] += rn_p0[2*i+1];
            }

            this->ax_centroids.p[ri] += initpos / static_cast<T>(bpa);

            // The target for axon centroids is defined by their origin location on the
            // retina. However, their target on the tectum is the retinal position
            // *transformed*. That means mirroring the retinal origins about the line
            // y=x. ALSO, if experimental manipulations have been made, then the target
            // positions will need to be modified accordingly.
            morph::Vector<T,2> tpos = this->ret->posn[ri];
            tpos.rotate(); // Achieves swapping of x and y coordinates
            this->ax_centroids.targ[ri].set_from (tpos);
            this->ax_centroids.mirrored = true;

            this->pending_branches[i].current = initpos.less_one_dim();
            this->pending_branches[i].id = i;
        }

        /*
         * Now arrange the order of the pending branches, so that the first ones in the
         * list are those closest to the centre of the retina. That means extracting the
         * branches by retina index or, more easily, by .target
         */
        // If retina has been ablated, then have to change xbreaks and y breaks
        if (this->conf->getBool ("ablate_ret_left", false) || this->conf->getBool ("ablate_tec_top", false)) {
            // Don't change pending branches at all
        } else {
            std::vector<B> pending_branches_reordered;
            // FIXME. Should this be 8 or bpa?
            morph::vVector<size_t> x_breaks = {this->ret->w/8, 2*(this->ret->w/8), 3*(this->ret->w/8), this->ret->w/2};
            size_t xstart = this->ret->w/2;

            for (auto xx : x_breaks) { std::cout << "x_break: " << xx << std::endl; }
            morph::vVector<size_t> y_breaks = {this->ret->h/8, 2*(this->ret->h/8), 3*(this->ret->h/8), this->ret->h/2};
            for (size_t i = 0; i < 4; ++i) {
                // copy elements from pending_branches in a square from
                // (w/2-x_breaks[i],h/2-y_breaks[i]) to (w/2+x_breaks[i],h/2+y_breaks[i])
                // into pending_branches_reordered.
                for (size_t yy = this->ret->h/2-y_breaks[i]; yy < this->ret->h/2+y_breaks[i]; ++yy) {
                    for (size_t xx = xstart-x_breaks[i]; xx < xstart+x_breaks[i]; ++xx) {
                        // Try to find branches with target xx,yy
                        morph::Vector<T,2> coord = this->ret->coord (xx, yy);
                        // Now go through pending_branches finding bpa branches to add to pending_branches_reordered
                        typename std::vector<B>::iterator it = this->pending_branches.begin();
                        while (it != this->pending_branches.end()) {
                            if ((it->target - coord).length() < 0.00001) {
                                pending_branches_reordered.push_back (*it);
                                it = this->pending_branches.erase (it);
                            } else {
                                ++it;
                            }
                        }
                    }
                }
                if (this->pb_sizes.empty()) {
                    this->pb_sizes.push_back (pending_branches_reordered.size());
                } else {
                    this->pb_sizes.push_back (pending_branches_reordered.size() - this->pb_sizes.sum());
                }
            }

            this->pending_branches.resize (pending_branches_reordered.size());
            this->pending_branches.swap (pending_branches_reordered);
        }

        /*
         * Setting upthe 'experiment suggests' information
         *
         * This requires some manual setup in the code below - the organisation of the
         * tissue after a manipulation is based on an interpretation of the various
         * papers in the literature.
         *
         * After initialising the ax_centroid targets from the retinal locations, we
         * have to modify ax_centroid's targ attribute, which gives the 'experiment
         * suggests' arrangement, if any of the experimental manipulations have been
         * applied.
         */

        // Have we had a graft swap?
        if (this->conf->getBool ("tectal_graftswap", false)) {
            std::cout << "tectal_graftswap..." << std::endl;
            this->ax_centroids.targ_graftswap (l1v, psv, l2v);
        }

        if (this->conf->getBool ("retinal_graftswap", false)) {
            std::cout << "WARNING: graft swap applied to retina, but axon centroids net was not updated with a prediction\n";
        }

        // Has tectum been rotated? If so, then modify the projected target - a
        // patch at l1v of size psv is rotated.
        unsigned int rots = 0;
        if ((rots = this->conf->getUInt ("tectal_rotations", 0)) > 0) {
            this->ax_centroids.targ_graftrotate (l1v, psv[0], rots);
        }

        if ((rots = this->conf->getUInt ("retinal_rotations", 0)) > 0) {
            std::cout << "WARNING: graft rotation applied to retina, but axon centroids net was not updated with a prediction\n";
        }

        if (this->conf->getBool ("compound_retina", false)) {
            // copy tissue::compound_tissue
            this->ax_centroids.targ_compound_tissue();
        }

        if (this->conf->getBool ("reber", false)) {
            this->ax_centroids.targ_reber();
        }

        if (conf->getBool ("singleaxon", false)) {
            // Now put the selected axon's expected final location in.
            if (!ax_centroids.targ.empty()) {
                ax_centroids.targ[0] = this->ret->posn[singleaxon_idx].plus_one_dim();
            }
        }

        if (this->conf->getBool ("ablate_ret_left", false)) {
            // expected layout has half as many locations, but they're stretched out into the full area?
            this->ax_centroids.targ_expand_topdown();
        }

        if (this->conf->getBool ("ablate_tec_top", false)) {
            this->ax_centroids.targ_squish_topdown();
        }

        /*
         * Now initialise the Visualisation code
         */

        // The min/max of rcpt[0] is used below to set a morph::Scale in branchvisual
        std::cout << "Receptor expression range: " << rcpt_min << " to " << rcpt_max << std::endl;

        // Which axons to see?
        Json::Value seelist = this->conf->getArray ("seeaxons");
        if (seelist.size() > 0) {
            seeaxons.clear();
            for (unsigned int i = 0; i < seelist.size(); ++i) {
                seeaxons.insert (seelist[i].asUInt());
            }
        }

        // Parameters settable from json
        this->m[0] = this->mconf->getDouble ("m1", 0.001);
        this->m[1] = this->mconf->getDouble ("m2", 0.2);
        this->m[2] = this->mconf->getDouble ("m3", 0.0);
        this->m[3] = this->mconf->getDouble ("mborder", 0.5);

        // Finally, set any addition parameters that will be needed with calling Agent1::run
        this->goslow = this->conf->getBool ("goslow", false);
    }

    static constexpr float widthtoheight = 0.5625f;
    // Initialise tissue visualisation
    void tvisinit()
    {
#ifdef VISUALISE
        const unsigned int ww = this->conf->getUInt ("win_width", 1200);
        this->tvv = new morph::Visual (ww, ww, "Retinal and Tectal expression");
        if constexpr (use_ortho) {
            this->tvv->ptype = morph::perspective_type::orthographic;
            this->tvv->ortho_bl = {-2,-2*widthtoheight};
            this->tvv->ortho_tr = {2,2*widthtoheight};
        }

        this->tvv->setCurrent();

        morph::Vector<float> offset = { -1.5f, -0.5f, 0.0f };
        morph::Vector<float> offset2 = offset;
        size_t show_pair = 0;

        // Retina
        offset2[1] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_exp, show_pair));
        offset2[0] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::ligand_exp, show_pair));
#ifdef SHOW_RET_GRADS
        offset2[0] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_grad_x, show_pair));
        offset2[1] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_grad_y, show_pair));
        offset2[1] -= 1.3f;
#endif
        // Tectum
        offset2[0] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_exp, show_pair));
        offset2[0] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_x, show_pair));
        offset2[1] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_y, show_pair));
        offset2[1] -= 1.3f;

        show_pair = 1;
        offset2[0] = offset[0];
        offset2[1] = offset[1] + 2.6f;
        // Retina
        offset2[1] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_exp, show_pair, true));
        offset2[0] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::ligand_exp, show_pair, true));
#ifdef SHOW_RET_GRADS
        offset2[0] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_grad_x, show_pair));
        offset2[1] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_grad_y, show_pair));
        offset2[1] -= 1.3f;
#endif
        // Tectum
        offset2[0] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_exp, show_pair, true));
        offset2[0] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_x, show_pair));
        offset2[1] += 1.3f;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_y, show_pair));
        offset2[1] -= 1.3f;
#endif
    }

    void visinit()
    {
#ifdef VISUALISE
        // morph::Visual init
        const unsigned int ww = this->conf->getUInt ("win_width", 1200);
        unsigned int wh = static_cast<unsigned int>(widthtoheight * (float)ww);
        std::cout << "New morph::Visual with width/height: " << ww << "/" << wh << std::endl;
        this->v = new morph::Visual (ww, ww, "Agent based retinotectal model");
        if constexpr (use_ortho) {
            this->v->ptype = morph::perspective_type::orthographic;
            this->v->ortho_bl = {-2,-2*widthtoheight};
            this->v->ortho_tr = {2,2*widthtoheight};
        }
        if (this->conf->getBool ("lighting", false)) { this->v->lightingEffects(); }

        // Offset for visuals
        morph::Vector<float> offset = { -1.5f, -0.5f, 0.0f };

        // Adding to the Main Visual second
        this->v->setCurrent();

        // Visualise the branches with a custom VisualModel
        this->bv = new BranchVisual<T, N, B> (v->shaderprog, v->tshaderprog, offset, &this->branches, &this->ax_history);
        this->bv->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->bv->target_scale.compute_autoscale (0, 1);
        this->bv->finalize();
        this->bv->addLabel ("Growth cones", {0.0f, 1.1f, 0.0f});
        v->addVisualModel (this->bv);

        // This one gives an 'axon view'
        offset[0] += 1.3f;
        this->av = new BranchVisual<T, N, B> (v->shaderprog, v->tshaderprog, offset, &this->branches, &this->ax_history);
        this->av->axonview = true;
        for (auto sa : this->seeaxons) { this->av->seeaxons.insert(sa); }
        this->av->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->av->target_scale.compute_autoscale (0, 1);
        this->av->finalize();
        this->av->addLabel ("Selected axons", {0.0f, 1.1f, 0.0f});
        v->addVisualModel (this->av);

        // Centroids of branches viewed with a NetVisual
        offset[0] += 1.3f;
        this->cv = new NetVisual<T> (v->shaderprog, v->tshaderprog, offset, &this->ax_centroids);
        this->cv->maxlen = this->conf->getDouble ("maxnetline", 1.0);
        this->cv->viewmode = netvisual_viewmode::actual;
        this->cv->finalize();
        this->cv->addLabel ("axon centroids", {0.0f, 1.1f, 0.0f});
        v->addVisualModel (this->cv);

        // Another NetVisual view showing the target locations
        offset[1] -= 1.3f;
        this->tcv = new NetVisual<T> (v->shaderprog, v->tshaderprog, offset, &this->ax_centroids);
        this->tcv->viewmode = netvisual_viewmode::targetplus;
        this->tcv->finalize();
        this->tcv->addLabel ("experiment suggests", {0.0f, 1.1f, 0.0f});
        v->addVisualModel (this->tcv);
        offset[1] += 1.3f;

        // A graph of the SOS diffs between axon position centroids and target positions from retina
        offset[0] += 1.5f;
        this->gv = new morph::GraphVisual<T> (v->shaderprog, v->tshaderprog, offset);
        this->gv->twodimensional = false;
        this->gv->setlimits (0, this->conf->getFloat ("steps", 1000),
                             0, this->conf->getFloat("graph_ymax", 200.0f));
        this->gv->policy = morph::stylepolicy::lines;
        this->gv->ylabel = "SOS";
        this->gv->xlabel = "Sim time";
        this->gv->prepdata ("SOS");
        this->gv->finalize();
        v->addVisualModel (this->gv);
#endif
    }

    // The axons to see - these will have their path information stored
    std::set<size_t> seeaxons = {21, 38, 189, 378, 361};
    // Branches per axon
    unsigned int bpa = 8;
    // Number of RGCs on a side
    unsigned int rgcside = 20;
    // If true, then slow things down a bit in the visualization
    bool goslow = false;
    // How many steps to store history. Note, we might choose not to show all of these
    // in a visualisation?
    static constexpr size_t history = 2000;
    // Access to a parameter configuration object for experiment settings
    morph::Config* conf;
    // A model parameter config object - experiment with splitting the config up
    morph::Config* mconf;
    // rgcside^2 RGCs, each with bpa axon branches growing.
    guidingtissue<T, N>* ret;
    // Same sized tectum tissue
    guidingtissue<T, N>* tectum;
    // Parameters vecto (See Table 2 in the paper)
    morph::Vector<T, 4> m = { T{0.02}, T{0.2}, T{0.15}, T{0.1} };
    // The centre coordinate
    morph::Vector<T,2> centre = { T{0.5}, T{0.5} }; // FIXME bit of a hack, this.
    // (rgcside^2 * bpa) branches, as per the paper
    std::vector<B> branches;
    // Branches are initialised in pending_branches, and introduced into branches in groups
    std::vector<B> pending_branches;
    // If pending_branches contains 'groups' of axons to introduce, then the sizes of each group are given in this container
    morph::vVector<size_t> pb_sizes;
    // Centroid of the branches for each axon
    rgcnet<T> ax_centroids;
    // Path history is a map indexed by axon id. 3D as it's used for vis.
    std::map<size_t, morph::vVector<morph::Vector<T, 3>>> ax_history;
    // Receptor max and min - used across init() and visinit()
    T rcpt_max = -1e9;
        T rcpt_min = 1e9;

#ifdef VISUALISE
    // A visual environment for the sim
    morph::Visual* v;
    // A visual environment specifically for the tissue visualisation
    morph::Visual* tvv;
    // Specialised visualization of agents as spheres with a little extra colour patch on top
    BranchVisual<T, N, B>* bv;
    // Another visualization to show axon paths with only a few axons
    BranchVisual<T, N, B>* av;
    // Centroid visual
    NetVisual<T>* cv;
    // Centroid visual for targets
    NetVisual<T>* tcv;
    // A graph for the SOS metric
    morph::GraphVisual<T>* gv;
#endif
};
