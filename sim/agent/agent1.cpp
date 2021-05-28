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
#include <morph/Visual.h>
#include <morph/GraphVisual.h>

#include "branchvisual.h"
#include "branch.h"
#include "netvisual.h"
#include "net.h"
#include "tissue.h"
#include "tissuevisual.h"

template<typename T, size_t N>
struct Agent1
{
    Agent1 (morph::Config* cfg)
    {
        this->conf = cfg;
        this->init();
    }
    ~Agent1() { delete this->ret; }

    static constexpr unsigned int showevery = 100;
    static constexpr unsigned int visevery = 5;

    //! Run this model!
    void run()
    {
        std::chrono::steady_clock::time_point laststep = std::chrono::steady_clock::now();

        typename std::vector<branch<T, N>>::iterator pending_br_it = this->pending_branches.begin();
        typename morph::vVector<size_t>::iterator pb_sz_it = this->pb_sizes.begin();
        // How often to introduce groups of axons? 0 means 'all at once'
        unsigned int intro_every = this->conf->getUInt ("intro_every", 0);
        if (intro_every == 0) {
            this->branches.resize (this->pending_branches.size());
            this->branches.swap (this->pending_branches);
        }

        for (unsigned int i = 0; i < this->conf->getUInt ("steps", 1000); ++i) {

            if (intro_every > 0 && pb_sz_it != this->pb_sizes.end() && i%intro_every == 0) {
                // Introduce some of pending_branches into branches
                typename std::vector<branch<T, N>>::iterator brit = this->branches.end();
                std::cout << "Adding " << (*pb_sz_it) << " to this->branches...\n";
                this->branches.resize (this->branches.size() + *pb_sz_it);
                brit = this->branches.end();
                brit -= *pb_sz_it;
                std::copy (pending_br_it, pending_br_it+*pb_sz_it, brit);
                pending_br_it += *pb_sz_it;
                pb_sz_it++;
            }

            this->step();
            if (i%visevery == 0) { this->vis(i); } // Visualize every 10 doubles time required for program
            if (i%showevery == 0) {
                std::chrono::steady_clock::duration since = std::chrono::steady_clock::now() - laststep;
                std::cout << "step " << i << ". Per step: "
                          << std::chrono::duration_cast<std::chrono::milliseconds>(since).count()/showevery << " ms\n";
                std::cout << "SOS of axon centroids: " << this->ax_centroids.sos() << std::endl;
                laststep = std::chrono::steady_clock::now();
            }
        }
        std::cout << "Done simulating\n";
        this->vis(this->conf->getUInt ("steps", 1000));
        this->v->keepOpen();
    }

    //! Update the visualization
    void vis (unsigned int stepnum)
    {
        if (this->goslow == true) {
            glfwWaitEventsTimeout (0.1); // to add artificial slowing
        } else {
            glfwPollEvents();
        }
        this->bv->reinit();
        this->av->reinit();
        this->cv->reinit();
        this->gv->append ((float)stepnum, this->ax_centroids.sos(), 0);
        this->v->render();
        if (this->conf->getBool ("movie", false)) {
            std::stringstream frame;
            frame << "log/agent/";
            frame.width(4);
            frame.fill('0');
            frame << stepnum;
            frame << ".png";
            this->v->saveImage(frame.str());
        }
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
        for (unsigned int i = 0; i < this->ret->num(); ++i) { this->ax_centroids.p[i] = {T{0}, T{0}, T{0}}; }
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

    //! Create a tissue visual, to reduce boilerplate code in init()
    tissuevisual<float, N>* createTissueVisual (morph::Vector<T,3>& offset, guidingtissue<T, N>* gtissue,
                                                const std::string& tag,
                                                expression_view exview, size_t pair_to_view)
    {
        tissuevisual<float, N>* tv = new tissuevisual<float, N>(v->shaderprog, v->tshaderprog, gtissue, offset);
        tv->view = exview;
        tv->pair_to_view = pair_to_view;
        tv->cm.setType (morph::ColourMapType::Duochrome);
        tv->cm.setHueRG(); // BG is nice
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

    //! Simulation init
    void init()
    {
        /*
         * Create tectum and retina tissue objects
         */

        this->rgcside = this->conf->getUInt ("rgcside", this->rgcside);
        T gr_denom = rgcside-1;
        T gr = T{1}/gr_denom; // gr is grid element length

        this->tectum = new guidingtissue<T, N>(this->rgcside, this->rgcside, {gr, gr}, {0.0f, 0.0f},
                                               (expression_form)this->conf->getUInt ("tectum_form", 0),
                                               (expression_form)this->conf->getUInt ("tectum_form", 0));

        this->ret = new guidingtissue<T, N>(this->rgcside, this->rgcside, {gr, gr}, {0.0f, 0.0f},
                                            (expression_form)this->conf->getUInt ("retina_form", 0),
                                            (expression_form)this->conf->getUInt ("retina_form", 0));

        /*
         * Apply any manipulations to retina or tectum
         */

        // Use a variable to prevent multiple manipulations from being applied. May
        // need to tweak this to allow selected combinations of manipulations.
        bool manipulated = false;

        // Graft swap manipulation
        Json::Value gs_coords = this->conf->getValue ("graftswap_coords");
        Json::Value l1 = gs_coords.get ("locn1", "[0,0]");
        Json::Value l2 = gs_coords.get ("locn2", "[0,0]");
        Json::Value ps = gs_coords.get ("patchsize", "[0,0]");
        if (l1.size() < 2 || l2.size() < 2 || ps.size() < 2) {
            throw std::runtime_error ("Bad values for locn1, locn2 or patchsize.");
        }

        morph::Vector<size_t, 2> l1v = { l1[0].asUInt(), l1[1].asUInt() };
        morph::Vector<size_t, 2> l2v = { l2[0].asUInt(), l2[1].asUInt() };
        morph::Vector<size_t, 2> psv = { ps[0].asUInt(), ps[1].asUInt() };
        std::cout << "Location 1: " << l1v << " Location 2: " << l2v << " patch size: " << psv << std::endl;

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

        // Various tissue ablation manipulations
        if (this->conf->getBool ("ablate_ret_right", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->ret->ablate_right_half();
            manipulated = true;
        }
        if (this->conf->getBool ("ablate_ret_left", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->ret->ablate_left_half();
            manipulated = true;
        }
        if (this->conf->getBool ("ablate_ret_top", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->ret->ablate_top_half();
            manipulated = true;
        }
        if (this->conf->getBool ("ablate_ret_bot", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->ret->ablate_bottom_half();
            manipulated = true;
        }

        if (this->conf->getBool ("ablate_tec_right", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->tectum->ablate_right_half();
            manipulated = true;
        }
        if (this->conf->getBool ("ablate_tec_left", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->tectum->ablate_left_half();
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
        if (this->conf->getBool ("ablate_tec_bot", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->tectum->ablate_bottom_half();
            manipulated = true;
        }

        // retinal/tectal receptor/ligand knockout manipulations
        int ko = this->conf->getInt ("knockout_ret_rcpt", -1);
        if (ko > -1) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->ret->receptor_knockout ((size_t)ko);
            manipulated = true;
        }
        ko = this->conf->getInt ("knockout_tec_rcpt", -1);
        if (ko > -1) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->tectum->receptor_knockout ((size_t)ko);
            manipulated = true;
        }
        ko = this->conf->getInt ("knockout_ret_lgnd", -1);
        if (ko > -1) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->ret->ligand_knockout ((size_t)ko);
            manipulated = true;
        }
        ko = this->conf->getInt ("knockout_tec_lgnd", -1);
        if (ko > -1) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            this->tectum->ligand_knockout ((size_t)ko);
            manipulated = true;
        }

        std::cout << "Retina has " << this->ret->num() << " cells\n";
        std::cout << "Tectum has " << this->tectum->num() << " cells\n";
        std::cout << "Retina is " << this->ret->w << " wide and " << this->ret->h << " high\n";

        /*
         * Set up the axon branches (Agent1::pending_branches, which will be copied into
         * Agent1::branches) and the ax_centroids net object used for visualisation/analysis.
         */

        this->bpa = this->conf->getUInt ("bpa", 8);
        this->pending_branches.resize(this->ret->num() * bpa);

        this->ax_centroids.init (this->ret->w, this->ret->h);
        // Axon initial positions x and y can be uniformly randomly selected...
        morph::RandUniform<T, std::mt19937> rng_x(T{0}, T{1.0});
        morph::RandUniform<T, std::mt19937> rng_y(T{-0.2}, T{0});
        // ...or set from the ideal position plus a random perturbation
        morph::RandNormal<T, std::mt19937> rng_p0(T{0}, T{0.1});
        // A normally distributed perturbation is added for each branch. SD=0.1.
        morph::RandNormal<T, std::mt19937> rng_p(T{0}, T{0.1});
        // Generate random number sequences all at once
        std::vector<T> rn_x = rng_x.get (this->ret->num());
        std::vector<T> rn_y = rng_y.get (this->ret->num());
        std::vector<T> rn_p = rng_p.get (this->ret->num() * 2 * bpa);
        std::vector<T> rn_p0 = rng_p0.get (this->ret->num() * 2 * bpa);
        T rcpt_max = -1e9;
        T rcpt_min = 1e9;
        bool totally_random = this->conf->getBool ("totally_random_init", false);

        // A loop to set up each branch object in pending_branches.
        for (unsigned int i = 0; i < this->pending_branches.size(); ++i) {
            // Set the branch's termination zone
            unsigned int ri = i/bpa; // retina index
            this->pending_branches[i].aid = (int)ri; // axon index
            this->pending_branches[i].rcpt = this->ret->rcpt[ri];
            this->pending_branches[i].target = this->ret->posn[ri];
            // Call the first interaction parameter 'EphA'
            rcpt_max =  this->pending_branches[i].rcpt[0] > rcpt_max ? pending_branches[i].rcpt[0] : rcpt_max;
            rcpt_min =  this->pending_branches[i].rcpt[0] < rcpt_min ? pending_branches[i].rcpt[0] : rcpt_min;

            // Set as in the S&G paper - starting at bottom in region x=(0,tectum->w), y=(-0.2,0)
            morph::Vector<T, 3> initpos;
            if (totally_random == true) {
                initpos = { rn_x[ri] + rn_p[2*i], rn_y[ri] + rn_p[2*i+1], 0 };
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

        std::vector<branch<T, N>> pending_branches_reordered;
        std::vector<size_t> x_breaks = {this->ret->w/8, 2*(this->ret->w/8), 3*(this->ret->w/8), this->ret->w/2};
        for (auto xx : x_breaks) { std::cout << "x_break: " << xx << std::endl; }
        std::vector<size_t> y_breaks = {this->ret->h/8, 2*(this->ret->h/8), 3*(this->ret->h/8), this->ret->h/2};
        for (size_t i = 0; i < 4; ++i) {
            // copy elements from pending_branches in a square from
            // (w/2-x_breaks[i],h/2-y_breaks[i]) to (w/2+x_breaks[i],h/2+y_breaks[i])
            // into pending_branches_reordered.
            for (size_t yy = this->ret->h/2-y_breaks[i]; yy < this->ret->h/2+y_breaks[i]; ++yy) {
                for (size_t xx = this->ret->w/2-x_breaks[i]; xx < this->ret->w/2+x_breaks[i]; ++xx) {
                    // Try to find branches with target xx,yy
                    morph::Vector<T,2> coord = this->ret->coord (xx, yy);
                    // Now go through pending_branches finding bpa branches to add to pending_branches_reordered
                    typename std::vector<branch<T, N>>::iterator it = this->pending_branches.begin();
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

        /*
         * After initialising the ax_centroid targets from the retinal locations, we
         * have to modify ax_centroid's targ attribute if any of the experimental
         * manipulations have been applied.
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
            std::cout << "Implement me!\n";
        }

        if (this->conf->getBool ("ablate_ret_left", false)) {
            // expected layout has half as many locations, but they're stretched out into the full area?
            std::cout << "Implement me!\n";
        }

        if (this->conf->getBool ("ablate_tec_top", false)) {
            std::cout << "Implement me!\n";
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
        this->m[0] = this->conf->getDouble ("m1", 0.02);
        this->m[1] = this->conf->getDouble ("m2", 0.2);
        this->m[2] = this->conf->getDouble ("m3", 0.15);
        this->m[3] = this->conf->getDouble ("mborder", 0.1);

        // morph::Visual init
        const unsigned int ww = this->conf->getUInt ("win_width", 1200);
        unsigned int wh = static_cast<unsigned int>(0.5625f * (float)ww);
        std::cout << "New morph::Visual with width/height: " << ww << "/" << wh << std::endl;
        this->v = new morph::Visual (ww, wh, "Seb's agent based retinotectal model");
        this->v->backgroundWhite();
        if (this->conf->getBool ("lighting", false)) { this->v->lightingEffects(); }

        // Offset for visuals
        morph::Vector<float> offset = { -1.5f, -0.5f, 0.0f };

        // Show a vis of the retina, to compare positions/colours
        morph::Vector<float> offset2 = offset;
        size_t show_pair = 0;

        // Retina
        offset2[1] += 1.3f;
        v->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_exp, show_pair));
        offset2[0] += 1.3f;
        v->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_grad_x, show_pair));
        offset2[1] += 1.3f;
        v->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_grad_y, show_pair));
        offset2[1] -= 1.3f;
        // Tectum
        offset2[0] += 1.3f;
        v->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_exp, show_pair));
        offset2[0] += 1.3f;
        v->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_x, show_pair));
        offset2[1] += 1.3f;
        v->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_y, show_pair));
        offset2[1] -= 1.3f;

        show_pair = 1;
        offset2[0] = offset[0];
        offset2[1] = offset[1] + 2.6f;
        // Retina
        offset2[1] += 1.3f;
        v->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_exp, show_pair));
        offset2[0] += 1.3f;
        v->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_grad_x, show_pair));
        offset2[1] += 1.3f;
        v->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_grad_y, show_pair));
        offset2[1] -= 1.3f;
        // Tectum
        offset2[0] += 1.3f;
        v->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_exp, show_pair));
        offset2[0] += 1.3f;
        v->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_x, show_pair));
        offset2[1] += 1.3f;
        v->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_y, show_pair));
        offset2[1] -= 1.3f;

        // Visualise the branches with a custom VisualModel
        this->bv = new BranchVisual<T, N> (v->shaderprog, v->tshaderprog, offset, &this->branches, &this->ax_history);
        this->bv->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->bv->target_scale.compute_autoscale (0, 1);
        this->bv->finalize();
        this->bv->addLabel ("Growth cones", {0.0f, 1.1f, 0.0f});
        v->addVisualModel (this->bv);

        // This one gives an 'axon view'
        offset[0] += 1.3f;
        this->av = new BranchVisual<T, N> (v->shaderprog, v->tshaderprog, offset, &this->branches, &this->ax_history);
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
        this->cv->viewmode = netvisual_viewmode::actual;
        this->cv->finalize();
        this->cv->addLabel ("axon centroids", {0.0f, 1.1f, 0.0f});
        v->addVisualModel (this->cv);

        // Another NetVisual view showing the target locations
        offset[1] -= 1.3f;
        this->tcv = new NetVisual<T> (v->shaderprog, v->tshaderprog, offset, &this->ax_centroids);
        this->tcv->viewmode = netvisual_viewmode::target;
        this->tcv->finalize();
        this->tcv->addLabel ("experiment predicts", {0.0f, 1.1f, 0.0f});
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

        // Finally, set any addition parameters that will be needed with calling Agent1::run
        this->goslow = this->conf->getBool ("goslow", false);
    }

    // The axons to see - these will have their path information stored
    std::set<size_t> seeaxons = {21, 38, 189, 378, 361};
    // Branches per axon
    unsigned int bpa = 8;
    // Number of RGCs on a side
    unsigned int rgcside = 21;
    // If true, then slow things down a bit in the visualization
    bool goslow = false;
    // How many steps to store history. Note, we might choose not to show all of these
    // in a visualisation?
    static constexpr size_t history = 2000;
    // Access to a parameter configuration object
    morph::Config* conf;
    // rgcside^2 RGCs, each with bpa axon branches growing.
    guidingtissue<T, N>* ret;
    // Same sized tectum tissue
    guidingtissue<T, N>* tectum;
    // Parameters vecto (See Table 2 in the paper)
    morph::Vector<T, 4> m = { T{0.02}, T{0.2}, T{0.15}, T{0.1} };
    // The centre coordinate
    morph::Vector<T,2> centre = { T{0.5}, T{0.5} }; // FIXME bit of a hack, this.
    // (rgcside^2 * bpa) branches, as per the paper
    std::vector<branch<T, N>> branches;
    // Branches are initialised in pending_branches, and introduced into branches in groups
    std::vector<branch<T, N>> pending_branches;
    // If pending_branches contains 'groups' of axons to introduce, then the sizes of each group are given in this container
    morph::vVector<size_t> pb_sizes;
    // Centroid of the branches for each axon
    rgcnet<T> ax_centroids;
    // Path history is a map indexed by axon id. 3D as it's used for vis.
    std::map<size_t, morph::vVector<morph::Vector<T, 3>>> ax_history;
    // A visual environment
    morph::Visual* v;
    // Specialised visualization of agents as spheres with a little extra colour patch on top
    BranchVisual<T, N>* bv;
    // Another visualization to show axon paths with only a few axons
    BranchVisual<T, N>* av;
    // Centroid visual
    NetVisual<T>* cv;
    // Centroid visual for targets
    NetVisual<T>* tcv;
    // A graph for the SOS metric
    morph::GraphVisual<T>* gv;
};

int main (int argc, char **argv)
{
    // Set up config object
    std::string paramsfile;
    if (argc >= 2) {
        paramsfile = std::string(argv[1]);
    } else {
        // Create an empty/default json file
        paramsfile = "./a1.json";
        morph::Tools::copyStringToFile ("{}\n", paramsfile);
    }

    morph::Config conf(paramsfile);
    if (!conf.ready) {
        std::cerr << "Failed to read config " << paramsfile << ". Exiting.\n";
        return 1;
    }

    size_t num_guiders = conf.getInt("num_guiders", 4);
    if (num_guiders == 4) {
        Agent1<float, 4> model (&conf);
        model.run();
    } else if (num_guiders == 2) {
        Agent1<float, 2> model (&conf);
        model.run();
    }

    return 0;
}
