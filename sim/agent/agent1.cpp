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
#include <vector>
#include <array>

#include <morph/CartGrid.h>
#include <morph/Config.h>
#include <morph/Random.h>
#include <morph/Visual.h>
#include <morph/CartGridVisual.h>

#include "branchvisual.h"
#include "branch.h"
#include "netvisual.h"
#include "net.h"
#include "tissue.h"
#include "retvisual.h"

template<typename T>
struct Agent1
{
    Agent1 (morph::Config* cfg)
    {
        this->conf = cfg;
        this->init();
    }
    ~Agent1() { delete this->ret; }

    void run()
    {
        for (unsigned int i = 0; i < this->conf->getUInt ("steps", 1000); ++i) {
            this->step();
            this->vis(i);
            if (i%100 == 0) { std::cout << "step " << i << "\n"; }
        }
        std::cout << "Done simulating\n";
        this->v->keepOpen();
    }

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

    void step()
    {
        // Compute the next position for each branch:
#ifdef __OSX__
        // Mac compiler didn't like omp parallel for in front of a for(auto...
#pragma omp parallel for
        for (unsigned int i = 0; i < this->branches.size(); ++i) {
            this->branches[i].compute_next (this->branches, this->m);
        }
#else
#pragma omp parallel for
        for (auto& b : this->branches) { b.compute_next (this->branches, this->m); }
#endif
        // Update centroids
        for (unsigned int i = 0; i < this->ret->num(); ++i) { this->ax_centroids.p[i] = {T{0}, T{0}, T{0}}; }
        for (auto& b : this->branches) {
            this->ax_centroids.p[b.aid][0] += b.next[0] / static_cast<T>(this->bpa);
            this->ax_centroids.p[b.aid][1] += b.next[1] / static_cast<T>(this->bpa);
        }
        // Once 'next' has been updated, add next to path:
        for (auto& b : this->branches) {
            b.path.push_back (b.next);
            if (b.path.size() > this->history) { b.pathfront++; }
        }
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

        this->ret = new retina<T>(this->rgcside, this->rgcside, {gr, gr}, {0.0f, 0.0f});

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
            this->branches[i].tz = {this->ret->posn[ri][0], this->ret->posn[ri][1]};
            // Set its ephrin interaction parameters (though these may be related to the tz)
            this->branches[i].EphA = T{1.05} + (T{0.26} * std::exp (T{2.3} * this->ret->posn[ri][0])); // R(x) = 0.26e^(2.3x) + 1.05,
            EphA_max =  this->branches[i].EphA > EphA_max ? branches[i].EphA : EphA_max;
            EphA_min =  this->branches[i].EphA < EphA_min ? branches[i].EphA : EphA_min;
            // Set as in the authors' paper - starting at bottom in region x=(0,1), y=(-0.2,0)
            morph::Vector<T, 3> initpos = { rn_x[ri] + rn_p[2*i], rn_y[ri] + rn_p[2*i+1], 0 };
            morph::Vector<T, 2> initpos2 = { rn_x[ri] + rn_p[2*i], rn_y[ri] + rn_p[2*i+1] };
            this->ax_centroids.p[ri] += initpos / static_cast<T>(bpa);
            this->branches[i].path.clear();
            this->branches[i].path.push_back (initpos2);
            this->branches[i].id = i;
        }
        // The min/max of EphA is used below to set a morph::Scale in branchvisual
        std::cout << "EphA range: " << EphA_min << " to " << EphA_max << std::endl;

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
        retvisual<float>* retv = new retvisual<float>(v->shaderprog, v->tshaderprog, ret, offset);
        retv->cm.setType (morph::ColourMapType::Duochrome);
        retv->cm.setHueRG();
        retv->finalize();
        v->addVisualModel (retv);

        // Visualise the branches with a custom VisualModel
        offset[0] += 1.3f;
        this->bv = new BranchVisual<T> (v->shaderprog, offset, &this->branches);
        this->bv->EphA_scale.compute_autoscale (EphA_min, EphA_max);
        this->bv->finalize();
        v->addVisualModel (this->bv);

        // This one gives an 'axon view'
        offset[0] += 1.3f;
        this->av = new BranchVisual<T> (v->shaderprog, offset, &this->branches);
        this->av->axonview = true;
        this->av->bpa = this->bpa;
        this->av->seeaxons.insert(21);
        this->av->seeaxons.insert(38);
        this->av->seeaxons.insert(189);
        this->av->seeaxons.insert(378);
        this->av->seeaxons.insert(361);
        this->av->EphA_scale.compute_autoscale (EphA_min, EphA_max);
        this->av->finalize();
        v->addVisualModel (this->av);

        // Centroids of branches viewed with a NetVisual
        offset[0] += 1.3f;
        this->cv = new NetVisual<T> (v->shaderprog, offset, &this->ax_centroids);
        this->cv->finalize();
        v->addVisualModel (this->cv);
    }

    std::vector<T> ephcolourdata;
    std::vector<morph::Vector<float, 3>> rgcposcolourdata;
    std::vector<morph::Vector<float, 3>> coords;

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
    retina<T>* ret;
    // Parameters vecto (See Table 2 in the paper)
    morph::Vector<T, 4> m = { T{0.02}, T{0.2}, T{0.15}, T{0.1} };
    // The centre coordinate
    morph::Vector<T,2> centre = { T{0.5}, T{0.5} }; // FIXME get from CartGrid
    // (rgcside^2 * bpa) branches, as per the paper
    std::vector<branch<T>> branches;
    // Centroid of the branches for each axon
    net<T> ax_centroids;
    // A visual environment
    morph::Visual* v;
    // Specialised visualization of agents with a history
    BranchVisual<T>* bv;
    // Another visualization to show axon paths with only a few axons
    BranchVisual<T>* av;
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
        paramsfile = "./sg.json";
        morph::Tools::copyStringToFile ("{}\n", paramsfile);
    }

    morph::Config conf(paramsfile);
    if (!conf.ready) {
        std::cerr << "Failed to read config " << paramsfile << ". Exiting.\n";
        return 1;
    }
    Agent1<float> model (&conf);
    model.run();

    return 0;
}
