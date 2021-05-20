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

        for (unsigned int i = 0; i < this->conf->getUInt ("steps", 1000); ++i) {
            this->step();
            if (i%visevery == 0) { this->vis(i); } // Visualize every 10 doubles time required for program
            if (i%showevery == 0) {
                std::chrono::steady_clock::duration since = std::chrono::steady_clock::now() - laststep;
                std::cout << "step " << i << ". Per step: "
                          << std::chrono::duration_cast<std::chrono::milliseconds>(since).count()/showevery << " ms\n";
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

    void init()
    {
        // Simulation init
        this->rgcside = this->conf->getUInt ("rgcside", this->rgcside);
        this->bpa = this->conf->getUInt ("bpa", 8);
        this->goslow = this->conf->getBool ("goslow", false);
        // gr is grid element length
        T gr_denom = rgcside-1;
        T gr = T{1}/gr_denom;
        std::cout << "Grid element length " << gr << std::endl;

        this->tectum = new guidingtissue<T, N>(this->rgcside, this->rgcside, {gr, gr}, {0.0f, 0.0f},
                                               (expression_form)this->conf->getUInt ("tectum_rcpt_form", 0),
                                               (expression_form)this->conf->getUInt ("tectum_lgnd_form", 0));

        this->ret = new guidingtissue<T, N>(this->rgcside, this->rgcside, {gr, gr}, {0.0f, 0.0f},
                                            (expression_form)this->conf->getUInt ("retina_rcpt_form", 0),
                                            (expression_form)this->conf->getUInt ("retina_lgnd_form", 0));

        // Extract the graftswap coordinates from the config
        Json::Value gs_coords = this->conf->getValue ("graftswap_coords");
        Json::Value l1 = gs_coords.get ("locn1", "[0,0]");
        Json::Value l2 = gs_coords.get ("locn2", "[0,0]");
        Json::Value ps = gs_coords.get ("patchsize", "[0,0]");
        if (l1.size() > 1 && l2.size() > 1 && ps.size() > 1) {

            morph::Vector<size_t, 2> l1v = { l1[0].asUInt(), l1[1].asUInt() };
            morph::Vector<size_t, 2> l2v = { l2[0].asUInt(), l2[1].asUInt() };
            morph::Vector<size_t, 2> psv = { ps[0].asUInt(), ps[1].asUInt() };

            std::cout << "Location 1: " << l1v << " Location 2: " << l2v << " patch size: " << psv << std::endl;

            if (this->conf->getBool ("tectal_graftswap", false)) {
                this->tectum->graftswap (l1v, psv, l2v);
            }

            if (this->conf->getBool ("retinal_graftswap", false)) {
                this->ret->graftswap (l1v, psv, l2v);
            }
        }

        std::cout << "Retina has " << this->ret->num() << " cells\n";
        this->branches.resize(this->ret->num() * bpa);

        std::cout << "Retina is " << this->ret->w << " wide and " << this->ret->h << " high\n";
        this->ax_centroids.init (this->ret->w, this->ret->h);
        // Axon initial positions x and y are uniformly randomly selected
        morph::RandUniform<T, std::mt19937> rng_x(T{0}, T{1.0});
        morph::RandUniform<T, std::mt19937> rng_y(T{-0.2}, T{0});
        // A normally distributed perturbation is added for each branch. SD=0.1.
        morph::RandNormal<T, std::mt19937> rng_p(T{0}, T{0.1});
        // Generate random number sequences all at once
        std::vector<T> rn_x = rng_x.get (this->ret->num());
        std::vector<T> rn_y = rng_y.get (this->ret->num());
        std::vector<T> rn_p = rng_p.get (this->ret->num() * 2 * bpa);
        T EphA_max = -1e9;
        T EphA_min = 1e9;
        for (unsigned int i = 0; i < this->branches.size(); ++i) {
            // Set the branch's termination zone
            unsigned int ri = i/bpa; // retina index
            this->branches[i].aid = (int)ri; // axon index
            this->branches[i].rcpt = this->ret->rcpt[ri];
            this->branches[i].target = this->ret->posn[ri];
            // Call the first interaction parameter 'EphA'
            EphA_max =  this->branches[i].rcpt[0] > EphA_max ? branches[i].rcpt[0] : EphA_max;
            EphA_min =  this->branches[i].rcpt[0] < EphA_min ? branches[i].rcpt[0] : EphA_min;
            // Set as in the S&G paper - starting at bottom in region x=(0,1), y=(-0.2,0)
            morph::Vector<T, 3> initpos = { rn_x[ri] + rn_p[2*i], rn_y[ri] + rn_p[2*i+1], 0 };
            morph::Vector<T, 2> initpos2 = { rn_x[ri] + rn_p[2*i], rn_y[ri] + rn_p[2*i+1] };
            this->ax_centroids.p[ri] += initpos / static_cast<T>(bpa);
            this->branches[i].current = initpos2;
            this->branches[i].id = i;
        }
        // The min/max of EphA (rcpt[0]) is used below to set a morph::Scale in branchvisual
        std::cout << "EphA range: " << EphA_min << " to " << EphA_max << std::endl;

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

        // Visualization init
        const unsigned int ww = this->conf->getUInt ("win_width", 1200);
        unsigned int wh = static_cast<unsigned int>(0.5625f * (float)ww);
        std::cout << "New morph::Visual with width/height: " << ww << "/" << wh << std::endl;
        this->v = new morph::Visual (ww, wh, "Seb's agent based retinotectal model");
        this->v->backgroundWhite();
        this->v->lightingEffects();

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
        this->bv->rcpt_scale.compute_autoscale (EphA_min, EphA_max);
        this->bv->target_scale.compute_autoscale (0, 1);
        this->bv->finalize();
        this->bv->addLabel ("Growth cones", {0.0f, 1.1f, 0.0f});
        v->addVisualModel (this->bv);

        // This one gives an 'axon view'
        offset[0] += 1.3f;
        this->av = new BranchVisual<T, N> (v->shaderprog, v->tshaderprog, offset, &this->branches, &this->ax_history);
        this->av->axonview = true;
        this->av->bpa = this->bpa;
        for (auto sa : this->seeaxons) { this->av->seeaxons.insert(sa); }
        this->av->rcpt_scale.compute_autoscale (EphA_min, EphA_max);
        this->av->target_scale.compute_autoscale (0, 1);
        this->av->finalize();
        this->av->addLabel ("Selected axons", {0.0f, 1.1f, 0.0f});
        v->addVisualModel (this->av);

        // Centroids of branches viewed with a NetVisual
        offset[0] += 1.3f;
        this->cv = new NetVisual<T> (v->shaderprog, v->tshaderprog, offset, &this->ax_centroids);
        this->cv->finalize();
        this->cv->addLabel ("axon centroids", {0.0f, 1.1f, 0.0f});
        v->addVisualModel (this->cv);
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
    morph::Vector<T,2> centre = { T{0.5}, T{0.5} }; // FIXME get from CartGrid
    // (rgcside^2 * bpa) branches, as per the paper
    std::vector<branch<T, N>> branches;
    // Centroid of the branches for each axon
    net<T> ax_centroids;
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
