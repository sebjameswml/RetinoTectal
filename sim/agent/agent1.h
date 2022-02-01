/*
 * Retinotectal model resembling one presented by Hugh Simpson and Geoffrey
 * Goodhill in "A simple model can unify a broad range of phenomena in retinotectal map
 * development", Biol Cybern (2011) 104:9-29
 *
 * I'm bringing the idea of variable interaction with signalling gradients and
 * competition to try to get rid of the non-biological part of Simpson & Goodhill's
 * work.
 */

// If true, turns parallel execution OFF and branch interaction debugging ON.
static constexpr bool debug_compute_branch = false;
// If true, then compile code to collect minimum and maximum interactions from branches
static constexpr bool branch_min_maxes = false;
static constexpr bool branch_min_max_i = false; // make one false
static constexpr bool branch_min_max_j = true; // and one true

#include <ostream>
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

#include "branch.h"
#include "branch_stochastic.h"
#include "branch_geb.h"
#include "net.h"
#include "tissue.h"

#ifdef VISUALISE
# include <morph/Visual.h>
# include <morph/GraphVisual.h>
# include "branchvisual.h"
# include "netvisual.h"
# include "tissuevisual.h"
static constexpr bool visualise = true;
#else
static constexpr bool visualise = false;
#endif

// Permit one using directive, to shorten morph::unicode
#include <morph/unicode.h>
using morph::unicode;

// Declare class
template<typename T> struct AgentMetrics;
// Declare operator<< function
template<typename T> std::ostream& operator<< (std::ostream& os, const AgentMetrics<T>& am);
// Define class. To hold information about an Agent1 model. Single run or multiple runs?
template<typename T>
struct AgentMetrics
{
    // An identifier
    std::string id = std::string("");
    // The number of agents
    size_t n_agents = 0;
    // Sum of squared errors
    morph::vVector<T> sos;
    // RMS error = sqrt(sos/n_agents)
    morph::vVector<T> rms;
    // Number of crossings
    morph::vVector<T> crosscount;
    // Output a string
    std::string str() const
    {
        std::stringstream ss;
        if (sos.empty()) {
            ss << "No metrics to report.";
        } else if (sos.size() == 1) {
            ss << "SOS: " << this->sos.back() << " RMS: " << this->rms.back() << " #X: " << this->crosscount.back();
        } else {
            ss << "SOS mean(std): " << this->sos.mean() << "(" << this->sos.std() <<  ") RMS: "
               << this->rms.mean() << "(" << this->rms.std() << ") #X: "
               << this->crosscount.mean() << "(" << this->crosscount.std() << ")";
        }
        return ss.str();
    }
    // operator +=
    void operator+= (const AgentMetrics<T>& rhs)
    {
        // Adding means combinding sos, rms and crosscount
        this->sos.insert (this->sos.end(), rhs.sos.begin(), rhs.sos.end());
        this->rms.insert (this->rms.end(), rhs.rms.begin(), rhs.rms.end());
        this->crosscount.insert (this->crosscount.end(), rhs.crosscount.begin(), rhs.crosscount.end());
    }

    // Overload the stream output operator
    friend std::ostream& operator<< <> (std::ostream& os, const AgentMetrics<T>& am);
};
// Define operator<< friend function
template <typename T>
std::ostream& operator<< (std::ostream& os, const AgentMetrics<T>& am)
{
    os << am.str();
    return os;
}

// A selection of possible graph layouts to show when running the program
enum class graph_layout { a, b, c, d, e, f, g, h, i, j, k };

// Agent1 coordinates an agent based simulation of axon branches of type B. This class
// also incorporates the visualisation code that displays the state of the simulation as
// it progresses.
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
        if constexpr (visualise == true) {
            if (this->v != nullptr) { delete this->v; }
        }
    }

    unsigned int showevery = 1000;
    unsigned int visevery = 5;

    //! Just show the tissue. Don't use at same time as run()
    void showtissue()
    {
        if constexpr (visualise == true) {
            // Get layout from config
            this->layout = (graph_layout)this->conf->getUInt ("tissue_layout", 0);
            this->tvisinit();
            this->tvv->render();
            this->tvv->keepOpen();
        } // else do nothing
    }

    //! Run this model!
    void run()
    {
        if constexpr (visualise == true) {
            this->visevery = this->conf->getUInt ("visevery", 5);
            if (this->visinit_done == false) {
                // Get layout from config file. Default to 0 or 'a'
                this->layout = (graph_layout)this->conf->getUInt ("graph_layout", 0);

                // Set up (from config file if necssary) the times at which the various
                // graphs will be frozen (i.e. no longer updated from the sim)
                this->freeze_times[0] = 0;
                this->freeze_times[1] = this->conf->getUInt ("freeze_time1", 20);
                this->freeze_times[2] = this->conf->getUInt ("freeze_time2", 80);
                this->freeze_times[3] = this->conf->getUInt ("freeze_time3", 120);

                // How early to start showing the crossings metric?
                this->crosscount_from = this->conf->getUInt ("crosscount_from", 1000);
                this->crosscount_every = this->conf->getUInt ("crosscount_every", 50);

                this->visinit();
            }
        }

        std::chrono::steady_clock::time_point laststep = std::chrono::steady_clock::now();

        typename std::vector<B>::iterator pending_br_it = this->pending_branches.begin();
        typename morph::vVector<size_t>::iterator pb_sz_it = this->pb_sizes.begin();
        // How often to introduce groups of axons? 0 means 'all at once'
        unsigned int intro_every = this->mconf->getUInt ("intro_every", 0);
        if (intro_every == 0) {
            this->branches.resize (this->pending_branches.size());
            std::copy (this->pending_branches.begin(), this->pending_branches.end(), this->branches.begin());
            //this->branches.swap (this->pending_branches);
        }

        // Are we running the random 'model'?
        if (this->mconf->getString ("model", "axgrad") == "random") {
            //std::cout << "branches.size(): " << this->branches.size() << std::endl;
            this->steprandom();
            //std::cout << "RMS error of axon centroids: " << this->ax_centroids.rms() << std::endl;
            if constexpr (visualise == true) {
                this->vis(1);
                if (this->immediate_exit == false) { this->v->keepOpen(); }
            }
            return;
        }

        this->gradient_rng = new morph::RandNormal<T, std::mt19937>(1, this->conf->getDouble("gradient_rng_width", 0.0));
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

            if constexpr (visualise == true) { if (i%this->visevery == 0) { this->vis(i); } }

            if (i%showevery == 0) {
                std::chrono::steady_clock::duration since = std::chrono::steady_clock::now() - laststep;
                std::cout << "Step " << i << ".";
                if (i) {
                    std::cout << " Step duration: "
                              << std::chrono::duration_cast<std::chrono::milliseconds>(since).count()/showevery
                              << " ms\n";
                    //std::cout << "RMS error of axon centroids: " << this->ax_centroids.rms() << std::endl;
                }
                laststep = std::chrono::steady_clock::now();
            }
        }
        delete this->gradient_rng;

        if constexpr (visualise == true) {

            // Update retinal NT position vs tectal RC position graph for the 'd' layout:
            if (this->layout == graph_layout::d || this->layout == graph_layout::e) {
                morph::vVector<T>nt;
                morph::vVector<T>rc;
                morph::vVector<T>nt_m; // manipulated
                morph::vVector<T>rc_m;

                size_t ii = 0;
                for (ii = 0; ii < this->ret->posn.size(); ++ii) {
                    if (this->ret->rcpt_manipulated[ii][0]) {
                        nt_m.push_back(1-this->ret->posn[ii][0]);
                        rc_m.push_back(this->ax_centroids.p[ii][1]);
                    } else {
                        nt.push_back(1-this->ret->posn[ii][0]);
                        rc.push_back(this->ax_centroids.p[ii][1]);
                    }
                }

                this->gv->setdata (nt, rc, "wt");
                T ki_amount = this->conf->getDouble ("knockin", 1);
                T kd_amount = this->conf->getDouble ("knockdown", 0);
                std::stringstream ss;
                ss << "ski:" << ki_amount << ", gkd:" << kd_amount;
                this->gv->setdata (nt_m, rc_m, ss.str());
                this->gv->reinit();
                this->v->render();

            } else if (this->layout == graph_layout::h || this->layout == graph_layout::k) {
                // Graph layout h/k shows the nt vs rc (for the mid band of the retina/central column of the tectum?).
                morph::vVector<T>nt;
                morph::vVector<T>rc;
                morph::vVector<T>rc_expt;
                size_t ii = 0;
                for (ii = 0; ii < this->ret->posn.size(); ++ii) {
                    nt.push_back(1-this->ret->posn[ii][0]);
                    rc.push_back(this->ax_centroids.p[ii][1]);
                    rc_expt.push_back(this->ax_centroids.targ[ii][1]);
                }
                // Also show expt. version.

                this->gv->setdata (nt, rc_expt, "Expt");
                this->gv->setdata (nt, rc, "Sim");
                this->gv->reinit();
                this->v->render();
            }

            // Save final image based on config file names
            std::stringstream nss;
            nss << this->imagedir << "/" << this->title << ".png";
            this->v->saveImage (nss.str());
            std::cout << "Saved image: " << nss.str() << std::endl;
            this->vis(this->conf->getUInt ("steps", 1000));
            AgentMetrics<T> am = this->get_metrics();
            std::cout << "Final SOS: " << am.sos << " and crosscount = " << am.crosscount << std::endl;
            if (this->immediate_exit == false) { this->v->keepOpen(); }
        }
    }
    std::string imagedir = std::string("./paper/images");

    //! Compute SOS, crossings count etc and return info
    AgentMetrics<T> get_metrics()
    {
        AgentMetrics<T> am;
        am.n_agents = this->ax_centroids.targ.size();
        am.sos.emplace_back (this->ax_centroids.sos());
        am.rms.emplace_back (this->ax_centroids.rms());
        am.crosscount.emplace_back (this->ax_centroids.crosscount());
        return am;
    }

    //! Save any relevant results of the simulation to an HdfData object.
    void save (const std::string& outfile)
    {
        std::cout << "save results (sos=" << this->ax_centroids.sos() << ")\n";
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

        switch (this->layout) {
        case graph_layout::a:
        {
            this->bv->reinit();  // Branches
            this->cv->reinit();  // Centroids to end
            this->tcv->reinit(); // Experiment
            this->av->reinit();  // Selected axons
            this->gv->append ((float)stepnum, this->ax_centroids.rms(), 0); // RMS/Crossings
            if (stepnum > this->crosscount_from && stepnum%this->crosscount_every == 0) {
                this->gv->append ((float)stepnum, this->ax_centroids.crosscount(), 1);
            }
            break;
        }
        case graph_layout::b:
        {
            this->bv->reinit(); // Branches
            this->cv->reinit(); // Centroids
            this->av->reinit(); // Selected axons
            break;
        }
        case graph_layout::c:
        {
            this->bv->reinit();  // Branches
            if (stepnum <= this->freeze_times[1]) { this->cv1->reinit(); } // Time freeze branches
            this->cv->reinit();  // Centroids
            this->tcv->reinit(); // Experiment
            this->av->reinit();  // Selected axons
            this->gv->append ((float)stepnum, this->ax_centroids.rms(), 0); // RMS/Crossings
            if (stepnum > this->crosscount_from && stepnum%this->crosscount_every == 0) {
                this->gv->append ((float)stepnum, this->ax_centroids.crosscount(), 1);
            }
            break;
        }
        case graph_layout::d:
        {
            this->bv->reinit();  // Branches
            this->cv->reinit();  // Centroids
            this->tcv->reinit(); // Experiment
            this->av->reinit();  // Selected axons
            break;
        }
        case graph_layout::e:
        {
            this->bv->reinit(); // Branches
            this->cv->reinit(); // Centroids
            break;
        }
        case graph_layout::f:
        {
            this->cv->reinit(); // Centroids
            this->sim_time_txt->setupText (std::to_string(stepnum)); // Text
            AgentMetrics<T> am = this->get_metrics();
            this->emetric_txt->setupText (am.rms.str());
            this->crossings_txt->setupText (am.crosscount.str());
            break;
        }
        case graph_layout::g:
        {
            this->bv->reinit(); // Branches
            this->cv->reinit(); // Centroids
            this->gv->append ((float)stepnum, this->ax_centroids.rms(), 0); // SOS (whole)
            this->gv->append ((float)stepnum, this->ax_centroids.rms_outside({0.25,0.25,0},{0.75,0.75,0}), 1); // SOS (surround)
            this->gv->append ((float)stepnum, this->ax_centroids.rms_inside({0.25,0.25,0},{0.75,0.75,0}), 2); // SOS (patch)
            // Show RMS error
            AgentMetrics<T> am = this->get_metrics();
            this->emetric_txt->setupText (std::to_string(am.rms[0]));
            break;
        }
        case graph_layout::h:
        {
            this->bv->reinit(); // Branches
            this->cv->reinit(); // Centroids
            // Show RMS error
            AgentMetrics<T> am = this->get_metrics();
            this->emetric_txt->setupText (std::to_string(am.rms[0]));
            break;
        }
        case graph_layout::i:
        {
            if (stepnum <= this->freeze_times[1]) { this->cv1->reinit(); } // Time freeze branches
            this->cv->reinit();  // Centroids
            this->av->reinit();  // Selected axons
            this->gv->append ((float)stepnum, this->ax_centroids.rms(), 0); // RMS/Crossings
            if (stepnum > this->crosscount_from && stepnum%this->crosscount_every == 0) {
                this->gv->append ((float)stepnum, this->ax_centroids.crosscount(), 1);
            }
            break;
        }
        case graph_layout::j:
        {
            this->cv->reinit(); // Centroids
            this->gv->append ((float)stepnum, this->ax_centroids.rms(), 0); // SOS (whole)
            this->gv->append ((float)stepnum, this->ax_centroids.rms_outside({0.25,0.25,0},{0.75,0.75,0}), 1); // SOS (surround)
            this->gv->append ((float)stepnum, this->ax_centroids.rms_inside({0.25,0.25,0},{0.75,0.75,0}), 2); // SOS (patch)
            // Show RMS error
            AgentMetrics<T> am = this->get_metrics();
            this->emetric_txt->setupText (std::to_string(am.rms[0]));
            break;
        }
        case graph_layout::k:
        {
            this->cv->reinit(); // Centroids
            // Show RMS error
            AgentMetrics<T> am = this->get_metrics();
            this->emetric_txt->setupText (std::to_string(am.rms[0]));
            break;
        }
        default:
            std::cerr << "Warning: Unknown graph_layout\n";
            break;
        }

        this->v->render();
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
#endif // VISUALISE

    //! One step of the simulation in which branches positions are randomly set
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

    //! RNG for adding noise to gradient sampling
    morph::RandNormal<T, std::mt19937>* gradient_rng = nullptr;

    //! Perform one step of the simulation
    void step()
    {
        if constexpr (branch_min_maxes == true) {
            // Reset minses/maxes for branches to compute min/max in this step
            for (auto& b : this->branches) {
                b.minses.set_from (std::numeric_limits<T>::max());
                b.maxes.set_from(std::numeric_limits<T>::lowest());
            }
        }

        // Compute the next position for each branch:
#ifdef __OSX__
        // Mac compiler didn't like omp parallel for in front of a for(auto...
#pragma omp parallel for
        for (unsigned int i = 0; i < this->branches.size(); ++i) {
            morph::Vector<T, 2*N> rns;
            this->gradient_rng->get (rns); // Hmmn. Is rng thread safe?
            this->branches[i].compute_next (this->branches, this->ret, this->tectum, this->m, rns);
        }
#else
        if constexpr (debug_compute_branch == true) {
            for (auto& b : this->branches) {
                morph::Vector<T, 2*N> rns;
                this->gradient_rng->get (rns); // Hmmn. Is rng thread safe?
                b.compute_next (this->branches, this->ret, this->tectum, this->m, rns);
            }
        } else {
#pragma omp parallel for
            for (auto& b : this->branches) {
                morph::Vector<T, 2*N> rns;
                this->gradient_rng->get (rns);
                b.compute_next (this->branches, this->ret, this->tectum, this->m, rns);
            }
#endif
        }

        if constexpr (branch_min_maxes == true) {
            // Check mins/maxes
            morph::Vector<T,N> minses;
            minses.set_from(std::numeric_limits<T>::max());
            morph::Vector<T,N> maxes;
            maxes.set_from(std::numeric_limits<T>::lowest());
            for (auto b : this->branches) {
                for (size_t i = 0; i < N; ++i) {
                    minses[i] = b.minses[i] < minses[i] ? b.minses[i] : minses[i];
                    maxes[i] = b.maxes[i] > maxes[i] ? b.maxes[i] : maxes[i];
                }
            }
            std::cout << "Minimum axon-axon signal: " << minses << std::endl;
            std::cout << "Maximum axon-axon signal: " << maxes << std::endl;
        }
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
    //! Add a suitable set of orientation labels to a square-map visual
    void addOrientationLabels (morph::VisualModel* vm, const std::string& tag)
    {
        // Tag tells us which of the orientation labels to draw
        if (tag == "Tectal") { // Caudal, Rostral, Medial, Lateral
            vm->addLabel (std::string("C"), {0.5f, 1.05f, 0.0f});
            vm->addLabel (std::string("R"), {0.5f, -0.1f, 0.0f});
            vm->addLabel (std::string("M"), {-0.1f, 0.5f, 0.0f});
            vm->addLabel (std::string("L"), {1.04f, 0.5f, 0.0f});
        } else { // Retina: Dorsal, Ventral, Temporal, Nasal
            vm->addLabel (std::string("D"), {0.5f, 1.05f, 0.0f});
            vm->addLabel (std::string("V"), {0.5f, -0.1f, 0.0f});
            vm->addLabel (std::string("T"), {-0.08f, 0.5f, 0.0f});
            vm->addLabel (std::string("N"), {1.05f, 0.5f, 0.0f});
        }
    }

    //! Create a tissue visual, to reduce boilerplate code in init()
    tissuevisual<float, N>* createTissueVisual (morph::Vector<T,3>& offset, guidingtissue<T, N>* gtissue,
                                                const std::string& tag,
                                                expression_view exview, size_t pair_to_view, int alt_cmap=0)
    {
        return this->createTissueVisual (this->tvv->shaderprog, this->tvv->tshaderprog,
                                         offset, gtissue, tag, exview, pair_to_view, alt_cmap);
    }

    //! Show x/y gradients in an alternative colour map (Twilight)?
    static constexpr bool grads_in_jet = false;

    //! Create a tissue visual, to reduce boilerplate code in init(). Use either shader/tshader progs
    tissuevisual<float, N>* createTissueVisual (GLuint shad_prog, GLuint tshad_prog,
                                                morph::Vector<T,3>& offset, guidingtissue<T, N>* gtissue,
                                                const std::string& tag,
                                                expression_view exview, size_t pair_to_view, int alt_cmap=0)
    {
        tissuevisual<float, N>* tv = new tissuevisual<float, N>(shad_prog, tshad_prog, gtissue, offset);
        tv->view = exview;
        tv->pair_to_view = pair_to_view;
        if (exview == expression_view::ligand_grad_x_single || exview == expression_view::ligand_grad_y_single) {
            if constexpr (grads_in_jet == false) {
                tv->cm.setType (morph::ColourMapType::Monochrome);
                // Hue from pair_to_view.
                switch (pair_to_view) {
                case 0:
                    tv->cm.setHue (1.0f);   // red
                    break;
                case 1:
                    tv->cm.setHue (0.333f); // green
                    break;
                case 2:
                    tv->cm.setHue (0.833f); // magenta?
                    break;
                case 3:
                    tv->cm.setHue (0.5f);   // cyan?
                    break;
                default:
                    tv->cm.setHue (0.15f);  // something else
                    break;
                }
            } else {
                tv->cm.setType (morph::ColourMapType::Twilight);
            }
        } else {
            tv->cm.setType (morph::ColourMapType::Duochrome);
            if (alt_cmap == 1) {
                tv->cm.setHueCM();
            } else if (alt_cmap == 2) {
                tv->cm.setHueRG(); // Special - for Retinal positions (hence RG map)
            } else {
                tv->cm.setHueRG(); // WAS setHueRB, but when mixed you get magenta, so need RG for these maps still
            }
        }
        std::stringstream ss;
        if (exview == expression_view::receptor_exp) {
            ss << tag << " receptor expression " << (pair_to_view*2) << "/" << (1+pair_to_view*2);
        } else if (exview == expression_view::receptor_grad_x) {
            ss << tag << " receptor gradient "   << (pair_to_view*2) << "/" << (1+pair_to_view*2) << " x";
        } else if (exview == expression_view::receptor_grad_y) {
            ss << tag << " receptor gradient "   << (pair_to_view*2) << "/" << (1+pair_to_view*2) << " y";
        } else if (exview == expression_view::ligand_exp) {
            ss << tag << " ligand expression "   << (pair_to_view*2) << "/" << (1+pair_to_view*2);
        } else if (exview == expression_view::ligand_grad_x) {
            ss << tag << " ligand gradient "     << (pair_to_view*2) << "/" << (1+pair_to_view*2) << " x";
        } else if (exview == expression_view::ligand_grad_x_single) {
            ss << tag << " ligand gradient "     << (pair_to_view) << " x";
        } else if (exview == expression_view::ligand_grad_y) {
            ss << tag << " ligand gradient "     << (pair_to_view*2) << "/" << (1+pair_to_view*2) << " y";
        } else if (exview == expression_view::ligand_grad_y_single) {
            ss << tag << " ligand gradient "     << (pair_to_view) << " y";
        } else if (exview == expression_view::cell_positions) {
            ss << tag << " cell positions";
        }
        tv->addLabel (ss.str(), {0.0f, 1.15f, 0.0f});

        this->addOrientationLabels (tv, tag);

        tv->finalize();

        return tv;
    }
#endif // VISUALISE

    static constexpr size_t singleaxon_idx = 210;
    // Set true to use orthographic projection
    static constexpr bool use_ortho = false;
    static constexpr bool use_ortho_tvv = false;

    morph::Vector<expression_direction, N> get_directions (const std::string& dirns_tag)
    {
        morph::Vector<expression_direction, N> dirns;
        for (auto& dd : dirns) { dd = expression_direction::x_increasing; }
        nlohmann::json arr = this->mconf->get (dirns_tag);
        if (arr.size() > 0) {
            for (unsigned int i = 0; i < arr.size(); ++i) {
                std::string d = arr[i].get<std::string>();
                if (d == "+x" || d == "x") {
                    dirns[i] = expression_direction::x_increasing;
                } else if (d == "-x") {
                    dirns[i] = expression_direction::x_decreasing;
                } else if (d == "+y" || d == "y") {
                    dirns[i] = expression_direction::y_increasing;
                } else if (d == "-y") {
                    dirns[i] = expression_direction::y_decreasing;
                } else {
                    std::stringstream ee;
                    ee << "Invalid entry in _directions: '" << d << "' should be +x, +y, -x or -y.";
                    throw std::runtime_error (ee.str());
                }
            }
        } else {
            std::stringstream ee;
            ee << "Need to set content of JSON parameter " << dirns_tag;
            throw std::runtime_error (ee.str());
        }
        return dirns;
    }

    morph::Vector<interaction, N> get_interactions (const std::string& interactions_tag)
    {
        morph::Vector<interaction, N> interactions;
        for (auto& ii : interactions) { ii = interaction::null; }
        nlohmann::json arr = this->mconf->get (interactions_tag);
        if (arr.size() > 0) {
            for (unsigned int i = 0; i < arr.size(); ++i) {
                int ai = arr[i];
                if (ai < 0) {
                    interactions[i] = interaction::repulsion;
                } else if (ai > 0) {
                    interactions[i] = interaction::attraction;
                } // else if 0 then leave as interaction::null
            }
        } else {
            std::stringstream ee;
            ee << "Need to set content of JSON parameter " << interactions_tag;
            throw std::runtime_error (ee.str());
        }
        return interactions;
    }

    morph::Vector<expression_form, N> get_forms (const std::string& tag)
    {
        morph::Vector<expression_form, N> function_forms;
        for (auto& ff : function_forms) { ff = expression_form::lin; }
        nlohmann::json arr = this->mconf->get (tag);
        if (arr.size() > 0) {
            for (unsigned int i = 0; i < arr.size(); ++i) {
                function_forms[i] = (expression_form)arr[i].get<unsigned int>();
            }
        } else {
            std::stringstream ee;
            ee << "Need to set content of JSON parameter " << tag;
            throw std::runtime_error (ee.str());
        }
        return function_forms;
    }

    bool randomly_seeded = true;
    void setup_pending_branches()
    {
        std::vector<T> rn_x;
        std::vector<T> rn_y;
        std::vector<T> rn_p;
        std::vector<T> rn_p0;
        if (randomly_seeded) {
            // Axon initial positions x and y can be uniformly randomly selected...
            morph::RandUniform<T, std::mt19937> rng_x(T{0}, T{1.0});
            morph::RandUniform<T, std::mt19937> rng_y(T{-0.2}, T{0}); // S&G
            //morph::RandUniform<T, std::mt19937> rng_y(T{0.0001}, T{0.2}); // All within field
            // ...or set from the ideal position plus a random perturbation
            morph::RandNormal<T, std::mt19937> rng_p0(T{0}, T{0.1});
            // A normally distributed perturbation is added for each branch. SD=0.1.
            morph::RandNormal<T, std::mt19937> rng_p(T{0}, T{0.1});
            // Generate random number sequences all at once
            size_t axc_sz = this->ax_centroids.p.size();
            rn_x = rng_x.get (axc_sz); // ax_centroids size?
            rn_y = rng_y.get (axc_sz);
            rn_p = rng_p.get (axc_sz * 2 * this->bpa);
            rn_p0 = rng_p0.get (axc_sz * 2 * this->bpa);
        } else {
            // Axon initial positions x and y can be uniformly randomly selected...
            morph::RandUniform<T, std::mt19937> rng_x(T{0}, T{1.0}, 1000);
            morph::RandUniform<T, std::mt19937> rng_y(T{-0.2}, T{0}, 2000); // S&G
            //morph::RandUniform<T, std::mt19937> rng_y(T{0.0001}, T{0.2}); // All within field
            // ...or set from the ideal position plus a random perturbation
            morph::RandNormal<T, std::mt19937> rng_p0(T{0}, T{0.1}, 3000);
            // A normally distributed perturbation is added for each branch. SD=0.1.
            morph::RandNormal<T, std::mt19937> rng_p(T{0}, T{0.1}, 4000);
            // Generate random number sequences all at once
            size_t axc_sz = this->ax_centroids.p.size();
            rn_x = rng_x.get (axc_sz); // ax_centroids size?
            rn_y = rng_y.get (axc_sz);
            rn_p = rng_p.get (axc_sz * 2 * this->bpa);
            rn_p0 = rng_p0.get (axc_sz * 2 * this->bpa);
        }
        bool totally_random = this->mconf->getBool ("totally_random_init", true);
        std::string branch_model = this->mconf->getString ("branch_model", "james_agent");

        std::array<float, 3> red = { 1.0f, 0.0f, 0.0f };
        std::array<float, 3> blue = { 0.0f, 0.0f, 1.0f };

        float r_conf = this->mconf->getFloat ("r", 0.05f);
        float r_c_conf = this->mconf->getFloat ("r_c", 0.0f);
        float r_j_conf = this->mconf->getFloat ("r_j", 0.0f);
        float r_i_conf = this->mconf->getFloat ("r_i", 0.0f);
        T s = this->mconf->getFloat ("s", 1.1f);
        // A loop to set up each branch object in pending_branches.
        for (unsigned int i = 0; i < this->pending_branches.size(); ++i) {
            // Set the branch's termination zone
            unsigned int ri = i/bpa; // retina index
            this->pending_branches[i].init();
            this->pending_branches[i].s = s;
            if constexpr (branch_min_maxes == true) {
                this->pending_branches[i].maxes.set_from(std::numeric_limits<T>::lowest());
                this->pending_branches[i].minses.set_from(std::numeric_limits<T>::max());
            }
            this->pending_branches[i].setr (r_conf);
            this->pending_branches[i].setr_c (r_c_conf);
            this->pending_branches[i].setr_j (r_j_conf);
            this->pending_branches[i].setr_i (r_i_conf);
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

            if (this->genetic_manipulation == true) { // genetic manipulation of retinal receptor 0
                // Set colour red or blue depending on if receptor 0 in the retina was manipulated or not.

                this->ax_centroids.clr[ri] = this->ret->rcpt_manipulated[ri][0] == true ? red : blue;
            }

            // "experiment suggests": The target for axon centroids is defined by their
            // origin location on the retina. However, their target on the tectum is the
            // retinal position *transformed*. ALSO, if experimental manipulations have
            // been made, then the target positions will need to be modified
            // accordingly.

            // To convert from retinal position to tectal position, x'=y and y'=x.
            morph::Vector<T,2> tpos = this->ret->posn[ri];
            tpos.rotate();

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

            //for (auto xx : x_breaks) { std::cout << "x_break: " << xx << std::endl; }
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
    }

    /*
     * Setting up the 'experiment suggests' information
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
    void setup_expt_suggests()
    {
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
            throw std::runtime_error("WARNING: graft rotation applied to retina, but axon centroids net was not updated with a prediction");
        }

        if (this->conf->getBool ("compound_retina", false)) {
            // copy tissue::compound_tissue
            this->ax_centroids.targ_compound_tissue();
        }

        if (this->conf->getBool ("reber", false)) {
            this->ax_centroids.targ_reber();
        }

        // Don't have a brown-specific target yet
        if (this->conf->getBool ("brown", false)) {
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
    }

    //! Re-compute initial positions of branches; reset content of graphs
    void reset()
    {
        size_t num_branches = this->ret->num() * this->bpa;
        if (this->pending_branches.size() != num_branches) { throw std::runtime_error ("num_branches is wrong."); }
        this->setup_pending_branches();
        this->setup_expt_suggests();
#ifdef VISUALISE
        this->ax_history.clear();
        if (this->gv != nullptr) { this->gv->clear(); }
#endif
    }

    //! Update the 'm' parameters by reading mconf
    void update_m()
    {
        this->m[0] = this->mconf->getDouble ("m_g", 0.001);    // G rcpt-lgnd (axon-tectum)
        this->m[1] = this->mconf->getDouble ("m_j", 0.0);      // J rcpt-lgnd (axon-axon)
        this->m[2] = this->mconf->getDouble ("m_i", 0.0);      // I rcpt-rcpt (axon-axon)
        this->m[3] = this->mconf->getDouble ("m_c", 0.0);      // C (if used)
        this->m[4] = this->mconf->getDouble ("mborder", 0.5); // B
    }

    //! Graftswap "locn1" parameter from JSON, stored as a Vector. Member as it's used
    //! in setup_expt_suggests
    morph::Vector<size_t, 2> l1v;
    //! Graftswap "locn2" parameter from JSON, stored as a Vector.
    morph::Vector<size_t, 2> l2v;
    //! Graftswap "patchsize" parameter from Json.
    morph::Vector<size_t, 2> psv;

    //! Simulation init
    void init()
    {
        /*
         * Create tectum and retina tissue objects
         */

        this->rgcside = this->mconf->getUInt ("rgcside", this->rgcside);
        T gr_denom = rgcside-1;
        T gr = T{1}/gr_denom; // gr is grid element length

        // Get tissue parameters - expression directions, forms, interactions - from JSON
        morph::Vector<expression_form, N> ret_receptor_forms = this->get_forms ("ret_receptor_forms");
        morph::Vector<expression_form, N> ret_ligand_forms = this->get_forms ("ret_ligand_forms");
        morph::Vector<expression_direction, N> ret_receptor_directions = this->get_directions ("ret_receptor_directions");
        morph::Vector<expression_direction, N> ret_ligand_directions = this->get_directions ("ret_ligand_directions");
        morph::Vector<interaction, N> ret_forward_interactions = this->get_interactions ("ret_forward_interactions");
        // Set reverse interactions same as forward for now:
        morph::Vector<interaction, N> ret_reverse_interactions = this->get_interactions ("ret_forward_interactions");
        // Retinal receptor-receptor interactions. This models rcpt[i]-to-rcpt[i]
        // interactions. What about rcpt[i]-to-rcpt[j]? Could be a 4x4 matrix.
        morph::Vector<interaction, N> ret_rcptrcpt_interactions = this->get_interactions ("ret_rcptrcpt_interactions");

        morph::Vector<expression_form, N> tectum_receptor_forms = this->get_forms ("tectum_receptor_forms");
        morph::Vector<expression_form, N> tectum_ligand_forms = this->get_forms ("tectum_ligand_forms");
        morph::Vector<expression_direction, N> tectum_receptor_directions =  this->get_directions ("tectum_receptor_directions");
        morph::Vector<expression_direction, N> tectum_ligand_directions =  this->get_directions ("tectum_ligand_directions");

        // Tectum interactions don't matter - they don't have an effect, nor are they
        // visualised, but a value is needed for guidingtissue constructor.
        morph::Vector<interaction, N> tectum_forward_interactions;
        for (auto& ii : tectum_forward_interactions) { ii = interaction::repulsion; }
        morph::Vector<interaction, N> tectum_reverse_interactions;
        for (auto& ii : tectum_reverse_interactions) { ii = interaction::repulsion; }

        if constexpr (N==4 || N==2) {
            this->ret = new guidingtissue<T, N>(this->rgcside, this->rgcside, {gr, gr}, {0.0f, 0.0f},
                                                ret_receptor_forms,
                                                ret_ligand_forms,
                                                ret_receptor_directions,
                                                ret_ligand_directions,
                                                ret_forward_interactions,
                                                ret_reverse_interactions,
                                                ret_rcptrcpt_interactions);

            this->tectum = new guidingtissue<T, N>(this->rgcside, this->rgcside, {gr, gr}, {0.0f, 0.0f},
                                                   tectum_receptor_forms,
                                                   tectum_ligand_forms,
                                                   tectum_receptor_directions,
                                                   tectum_ligand_directions,
                                                   tectum_forward_interactions,
                                                   tectum_reverse_interactions,
                                                   tectum_reverse_interactions); // a dummy arg really
        } else {
            // C++-20 mechanism to trigger a compiler error for the else case. Not user friendly!
            []<bool flag = false>() { static_assert (flag, "N must be 2 or 4"); }();
        }

        /*
         * Apply any manipulations to retina or tectum
         */

        // Use a variable to prevent multiple manipulations from being applied. May
        // need to tweak this to allow selected combinations of manipulations.
        bool manipulated = false;

        // Graft swap manipulation
        nlohmann::json gs_coords = this->conf->root["graftswap_coords"];
        std::cout << gs_coords << std::endl;
        nlohmann::json l1;
        if (gs_coords.contains("locn1")) {
            l1 = gs_coords["locn1"]; // def [0,0]
        } else {
            l1.push_back (0); l1.push_back (0);
        }
        nlohmann::json l2 = gs_coords["locn2"]; //gs_coords.get ("locn2", "[0,0]");
        nlohmann::json ps = gs_coords["patchsize"];// gs_coords.get ("patchsize", "[0,0]");
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
            l1v = { l1[0].get<unsigned int>(), l1[1].get<unsigned int>() };
            l2v = { l2[0].get<unsigned int>(), l2[1].get<unsigned int>() };
            psv = { ps[0].get<unsigned int>(), ps[1].get<unsigned int>() };
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
            l1v = { l1[0].get<unsigned int>(), l1[1].get<unsigned int>() };
            psv = { ps[0].get<unsigned int>(), ps[1].get<unsigned int>() };
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
        this->genetic_manipulation = false;
        T affected = T{0.5};
        T ki_amount = this->conf->getDouble ("knockin", 1);
        T kd_amount = this->conf->getDouble ("knockdown", 0);
        if (this->conf->getBool ("reber", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            // Knockin/knockdown receptor 0:
            this->ret->receptor_knockin (0, affected, ki_amount);
            this->ret->receptor_knockdown (0, kd_amount);
            manipulated = true;
            this->genetic_manipulation = true;
        }

        if (this->conf->getBool ("brown", false)) {
            if (manipulated) { throw std::runtime_error ("Code is only tested for one manipulation at a time!"); }
            // Knockin receptor 0 in half of RGCs (but with no knockdown):
            this->ret->receptor_knockin (0, affected, ki_amount);
            manipulated = true;
            this->genetic_manipulation = true;
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

        this->setup_pending_branches();
        this->setup_expt_suggests();

        /*
         * Now initialise the Visualisation code
         */

        // The min/max of rcpt[0] is used below to set a morph::Scale in branchvisual
        std::cout << "Receptor expression range: " << rcpt_min << " to " << rcpt_max << std::endl;

        // Which axons to see?
        nlohmann::json seelist = this->conf->get ("seeaxons");
        if (seelist.size() > 0) {
            seeaxons.clear();
            for (unsigned int i = 0; i < seelist.size(); ++i) {
                seeaxons.insert (seelist[i].get<unsigned int>());
            }
        }

        // Parameters settable from json
        this->update_m();

        // Finally, set any additional parameters that will be needed with calling Agent1::run
        this->goslow = this->conf->getBool ("goslow", false);
        this->immediate_exit = this->conf->getBool ("exit", false);
    }

#ifdef VISUALISE
    static constexpr bool show_tectal_receptors = false;
    static constexpr float widthtoheight = 0.5625f;

    // Initialise tissue visualisation
    void tvisinit()
    {
        if (this->layout == graph_layout::a) {
            this->tvislayout_a();
        } else if (this->layout == graph_layout::g) {
            this->tvislayout_g();
        } else if (this->layout == graph_layout::h) {
            this->tvislayout_h();
        } else {
            this->tvislayout_b();
        }
    }

    // 2x3 graphs
    void tvislayout_a()
    {
        const unsigned int ww = this->conf->getUInt ("win_width", 1800);
        const unsigned int wh = this->conf->getUInt ("win_height", 1200);
        this->tvv = new morph::Visual (ww, wh, "Retinal and Tectal expression");
        this->tvv->setSceneTransXY(-0.405548,-0.139761);
        this->tvv->setSceneTransZ(-5.1);
        if constexpr (use_ortho_tvv) {
            float orthoside = 2.5;
            this->tvv->ptype = morph::perspective_type::orthographic;
            this->tvv->ortho_bl = {-orthoside,-orthoside};
            this->tvv->ortho_tr = {orthoside,orthoside};
        }

        this->tvv->setCurrent();

        morph::Vector<float> offset = { -1.5f, -2.5f, 0.0f };
        morph::Vector<float> offset2 = offset;
        size_t show_pair = 0; // 0 means show gradients for receptor/ligands 0 and 1.

        float sqside = 1.4f;

        // Retina
        offset2[1] += sqside;
        tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_exp, show_pair));
        offset2[0] += sqside;
        tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::ligand_exp, show_pair));
#ifdef SHOW_RET_GRADS
        offset2[0] += sqside;
        tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_grad_x, show_pair));
        offset2[1] += sqside;
        tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_grad_y, show_pair));
        offset2[1] -= sqside;
#endif
        // Tectum
        if constexpr (show_tectal_receptors == true) {
            offset2[0] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::receptor_exp, show_pair));
        }
        offset2[0] += sqside;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_exp, show_pair));
        // Tectal gradients for 0/1
        offset2[0] += sqside;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_x, show_pair));
        offset2[1] += sqside;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_y, show_pair));
        offset2[1] -= sqside;

        if constexpr (N>2) {
            show_pair = 1; // 1 means show for 2 and 3.
            offset2[0] = offset[0];
            offset2[1] = offset[1] + sqside;
            // Retina
            offset2[1] += sqside;
            // true means cyan-magenta
            tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_exp, show_pair, 1));
            offset2[0] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::ligand_exp, show_pair, 1));
#ifdef SHOW_RET_GRADS
            offset2[0] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_grad_x, show_pair));
            offset2[1] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, ret, "Retinal", expression_view::receptor_grad_y, show_pair));
            offset2[1] -= sqside;
#endif
            // Tectum
            if constexpr (show_tectal_receptors == true) {
                offset2[0] += sqside;
                tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::receptor_exp, show_pair, 1));
            }
            offset2[0] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_exp, show_pair, 1));
            offset2[0] += sqside;
            offset2[1] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_x, show_pair));
            offset2[1] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_y, show_pair));
            offset2[1] -= sqside;
        }
    }

    // Alt layout for ligand expression only. 2x5. expression and gradients
    void tvislayout_b()
    {
        const unsigned int ww = this->conf->getUInt ("win_width", 1200);
        const unsigned int wh = this->conf->getUInt ("win_height", 600);
        this->tvv = new morph::Visual (ww, wh, "Tissuevisb");
        this->tvv->setSceneTransXY(-0.405548,-0.139761);
        this->tvv->setSceneTransZ(-5.1);
        this->tvv->setCurrent();

        morph::Vector<float> offset = { -1.5f, 0.0f, 0.0f };
        morph::Vector<float> offset2 = offset;
        size_t show_pair = 0; // 0 means show gradients for ligands 0 and 1.
        float sqside = 1.4f;

        // Tectal expression for 0/1
        //offset2[1] -= sqside;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_exp, show_pair));
        // Tectal gradients for 0/1
        offset2[0] += sqside;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_x_single, 0));
        offset2[0] += sqside;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_x_single, 1));
        if constexpr (N>2) {
            offset2[0] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_x_single, 2));
            offset2[0] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_x_single, 3));
        }

        offset2 = offset;
        offset2[1] -= sqside;
        if constexpr (N>2) {
            show_pair = 1; // 1 means show for ligands 2 and 3.
            tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_exp, show_pair, 1));
        }
        offset2[0] += sqside;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_y_single, 0));
        offset2[0] += sqside;
        tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_y_single, 1));
        if constexpr (N>2) {
            offset2[0] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_y_single, 2));
            offset2[0] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_y_single, 3));
        }

        if constexpr (N>2) {
#if 0
            show_pair = 1; // 1 means show for ligands 2 and 3.
            offset2 = offset;
            tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_exp, show_pair, 1));
            offset2[0] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_x_single, show_pair));
            offset2[0] += sqside;
            tvv->addVisualModel (this->createTissueVisual (offset2, tectum, "Tectal", expression_view::ligand_grad_y, show_pair));
#endif
        }
    }

    // Alt layout for tectal ligand expression 0/1 only. 1x1
    void tvislayout_g()
    {
        const unsigned int ww = this->conf->getUInt ("win_width", 800);
        const unsigned int wh = this->conf->getUInt ("win_height", 900);
        this->tvv = new morph::Visual (ww, wh, "Tissuevisg");
        this->tvv->setSceneTrans(-0.485459f, -0.508987308f, -2.900002f);
        this->tvv->setCurrent();
        morph::Vector<float> offset = { 0.0f, 0.0f, 0.0f };
        tvv->addVisualModel (this->createTissueVisual (offset, tectum, "Tectal", expression_view::ligand_exp, 0));
    }

    void tvislayout_h()
    {
        const unsigned int ww = this->conf->getUInt ("win_width", 800);
        const unsigned int wh = this->conf->getUInt ("win_height", 900);
        this->tvv = new morph::Visual (ww, wh, "Tissuevish");
        this->tvv->setSceneTrans(-0.485459f, -0.508987308f, -2.900002f);
        this->tvv->setCurrent();
        morph::Vector<float> offset = { 0.0f, 0.0f, 0.0f };
        tvv->addVisualModel (this->createTissueVisual (offset, ret, "Retinal", expression_view::receptor_exp, 0));
    }

    // Updatable simulation time text
    morph::VisualTextModel* sim_time_txt = nullptr;
    morph::VisualTextModel* emetric_txt = nullptr;
    morph::VisualTextModel* crossings_txt = nullptr;

    // Set up the simulation visualisation scene. This depends on whether this->layout
    // is graph_layout::a, ::b or ::c, etc.
    //
    // 0 a: 2x3
    // 1 b: 1x3 branches, centroids, selected
    // 2 c: 2x4 layout with staggered time "snapshots"
    // 3 d: 2x3 like standard but with final DV-NT graph for e.g. Reber result
    // 4 e: 1x3 Branches, centroids and the DV-NT graph
    // 5 f: 1x1 with inset (just centroids and expt inset)
    // 6 g: 1x4 expt, branches, centroids, selected
    void visinit()
    {
        // morph::Visual init
        unsigned int wdefault = 1200;
        unsigned int hdefault = 800;
        switch (this->layout) {
        case graph_layout::a:
        case graph_layout::d:
        {
            // 2x3
            wdefault = 1920; hdefault = 1200;
            break;
        }
        case graph_layout::b:
        case graph_layout::e:
        case graph_layout::g:
        case graph_layout::h:
        case graph_layout::i:
        {
            // 1x4
            wdefault = 2480; hdefault = 574;
            break;
        }
        case graph_layout::c:
        {
            wdefault = 2200; hdefault = 1180;
            break;
        }
        case graph_layout::f:
        {
            // 1x1
            wdefault = 675; hdefault = 550;
            break;
        }
        case graph_layout::j:
        case graph_layout::k:
        {
            // 1x2
            wdefault = 1240; hdefault = 630;
            break;
        }
        default:
            break;
        }

        const unsigned int ww = this->conf->getUInt ("win_width", wdefault);
        const unsigned int wh = this->conf->getUInt ("win_height", hdefault);
        if (ww != wdefault || wh != hdefault) {
            std::cout << "visinit() with win height: " << wh << ", and win width: " << ww << " (from json)."<< std::endl;
        }

        std::string tt("Agent based retinotectal model: ");
        tt += this->title;
        this->v = new morph::Visual (ww, wh, tt);

        if (this->layout == graph_layout::a || this->layout == graph_layout::d) {
            this->v->setSceneTrans (-0.3762f, 0.6880f, -5.3f);
        } else if (this->layout == graph_layout::b) {
            this->v->setSceneTrans (-0.2375f, -0.03f, -2.8f);
        } else if (this->layout == graph_layout::c) {
            this->v->setSceneTrans (-0.91053f, 0.7895f, -6.1f);
        } else if (this->layout == graph_layout::e) {
            this->v->setSceneTrans (-0.890077f, -0.0236414f, -2.6f);
        } else if (this->layout == graph_layout::f) {
            this->v->setSceneTrans (1.00190961f,0.0175217576f,-2.70000315f);
        } else if (this->layout == graph_layout::g) {
            this->v->setSceneTrans (-1.22557163,-0.0104022622,-2.60000014);
        } else if (this->layout == graph_layout::h) {
            this->v->setSceneTrans (-1.22557163,-0.0104022622,-2.60000014);
        } else if (this->layout == graph_layout::i) {
            this->v->setSceneTrans (-1.22557163,-0.0104022622,-2.60000014);
        } else if (this->layout == graph_layout::j) {
            this->v->setSceneTrans (0.309693903,-0.021264188,-2.8000021);
        } else if (this->layout == graph_layout::k) {
            this->v->setSceneTrans (0.309693903,-0.021264188,-2.8000021);
        }

        if constexpr (use_ortho) {
            this->v->ptype = morph::perspective_type::orthographic;
            this->v->ortho_bl = {-2,-2*widthtoheight};
            this->v->ortho_tr = {2,2*widthtoheight};
        }

        if (this->conf->getBool ("lighting", false)) { this->v->lightingEffects(); }

        this->v->setCurrent(); // Don't have 2 windows, so this probably unnecessary
        this->set_graph_layout();
        this->visinit_done = true;
    }

    // This chooses which graph_layout function to call
    void set_graph_layout()
    {
        morph::Vector<float> offset = { -1.5f, -0.5f, 0.0f };
        if (this->layout == graph_layout::a           // 2 rows, 3 cols. Standard layout for investigations
            || this->layout == graph_layout::d) {     // Standard layout tweaked with graphs like in Brown et al
            this->graph_layout_ad (offset);
        } else if (this->layout == graph_layout::b) { // 1x3; branches, centroids, selected
            this->graph_layout_b (offset);
        } else if (this->layout == graph_layout::c) { // 2x4 layout with diff. time end points
            this->graph_layout_c (offset);
        } else if (this->layout == graph_layout::e) { // Branches, centroids and the position graph
            this->graph_layout_e (offset);
        } else if (this->layout == graph_layout::f) { // Just centroids. Used by agent1_eval.cpp
            this->graph_layout_f (offset);
        } else if (this->layout == graph_layout::g) { // 1x4; expt, branches, centroids, selected
            const std::string sl = this->conf->getString ("startletter", "A");
            this->graph_layout_g (offset, sl);
        } else if (this->layout == graph_layout::h) { // 1x4; expt, branches, centroids, ret NT vs tec RC
            const std::string sl = this->conf->getString ("startletter", "A");
            this->graph_layout_h (offset, sl);
        } else if (this->layout == graph_layout::i) {
            this->graph_layout_i (offset);
        } else if (this->layout == graph_layout::j) { // 1x2; centroids, selected
            const std::string sl = this->conf->getString ("startletter", "A");
            this->graph_layout_j (offset, sl);
        } else if (this->layout == graph_layout::k) { // 1x2; centroids, ret NT vs tec RC
            const std::string sl = this->conf->getString ("startletter", "A");
            this->graph_layout_k (offset, sl);
        } else {
            throw std::runtime_error ("Unknown layout");
        }
    }

    // 2 rows, 3 cols. Standard layout for investigations, or similar, tweaked with graphs for Reber/Brown
    void graph_layout_ad (const morph::Vector<float>& offset0)
    {
        morph::Vector<float> offset = offset0;
        // Top left
        v->addVisualModel (this->createTissueVisual (v->shaderprog, v->tshaderprog, offset, ret, "Retinal", expression_view::cell_positions, 0, 2));

        offset[1] -= 1.4f;
        // Branches: Visualise the branches with a custom VisualModel
        this->bv = new BranchVisual<T, N, B> (v->shaderprog, v->tshaderprog, offset, &this->branches, &this->ax_history);
        this->bv->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->bv->target_scale.compute_autoscale (0, 1);
        //this->bv->view = branchvisual_view::detail;
        this->bv->view = branchvisual_view::discs;
        //this->bv->view = branchvisual_view::discint;
        this->bv->finalize();
        this->bv->addLabel ("Branches", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->bv, std::string("Tectal"));
        v->addVisualModel (this->bv);

        // Axon centroids: Centroids of branches viewed with a NetVisual
        offset[0] += 1.3f;
        this->cv = new NetVisual<T> (v->shaderprog, v->tshaderprog, offset, &this->ax_centroids);
        this->cv->maxlen = this->conf->getDouble ("maxnetline", 1.0);
        this->cv->viewmode = netvisual_viewmode::actual;
        this->cv->finalize();
        this->cv->addLabel ("Axon centroids", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv, std::string("Tectal"));
        v->addVisualModel (this->cv);

        // Experiment: Another NetVisual view showing the target locations
        offset[1] += 1.4f;
        this->tcv = new NetVisual<T> (v->shaderprog, v->tshaderprog, offset, &this->ax_centroids);
        this->tcv->viewmode = netvisual_viewmode::targetplus;
        this->tcv->finalize();
        this->tcv->addLabel ("Experiment", {0.0f, 1.1f, 0.0f});
        v->addVisualModel (this->tcv);

        // Selected axons: This one gives an 'axon view'
        offset[0] += 1.5f;
        this->av = new BranchVisual<T, N, B> (v->shaderprog, v->tshaderprog, offset, &this->branches, &this->ax_history);
        this->av->view = branchvisual_view::axonview;
        for (auto sa : this->seeaxons) { this->av->seeaxons.insert(sa); }
        this->av->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->av->target_scale.compute_autoscale (0, 1);
        this->av->finalize();
        this->av->addLabel ("Selected axons", {0.0f, 1.1f, 0.0f});
        v->addVisualModel (this->av);

        offset[1] -= 1.4f;

        // Graph: A graph of the SOS diffs between axon position centroids and target positions from retina
        this->gv = new morph::GraphVisual<T> (v->shaderprog, v->tshaderprog, offset);
        this->gv->twodimensional = false;
        if (this->layout == graph_layout::a) {
            this->gv->axisstyle = morph::axisstyle::twinax;
            this->gv->setsize (0.9f, 1.0f);
            this->gv->setlimits (0, this->conf->getFloat ("steps", 1000),
                                 0, this->conf->getFloat("graph_ymax", 1.0f),
                                 0, this->conf->getFloat("graph_ymax2", 200.0f));
            this->gv->axislabelgap = 0.03f;
            this->gv->policy = morph::stylepolicy::lines;
            this->gv->ylabel = unicode::toUtf8 (unicode::epsilon);
            this->gv->ylabel2 = unicode::toUtf8 (unicode::eta);
            this->gv->xlabel = "t";
            this->gv->prepdata (unicode::toUtf8 (unicode::epsilon));
            this->gv->prepdata (unicode::toUtf8 (unicode::eta), morph::axisside::right);
        } else {
            this->gv->setlimits (0, 1, 0, 1);
            this->gv->policy = morph::stylepolicy::markers;
            this->gv->ylabel = "R " + unicode::toUtf8 (unicode::longrightarrow) + " tectum " + unicode::toUtf8 (unicode::longrightarrow) + " C";
            this->gv->xlabel = "N " + unicode::toUtf8 (unicode::longrightarrow) + " retina " + unicode::toUtf8 (unicode::longrightarrow) + " Y";
            //this->gv->prepdata ("0"); // without this, GraphVisual code crashes at first render.
        }
        this->gv->finalize();
        v->addVisualModel (this->gv);

        // A 'text' only visual model to list the parameters and incorporate figure letters
        offset[0] += 1.4f;
        morph::VisualModel* jtvm = new morph::VisualModel (v->shaderprog, v->tshaderprog, offset);
        float ty = 1.0f; // text y position
        float th = 0.1f; // text height
        float cw = 0.17f;
        jtvm->addLabel ("Parameters", {0.0f, ty, 0.0f});
        ty -= th;

        if (this->m[0] > 0) {
            jtvm->addLabel ("m_g", {0.0f, ty, 0.0f});
            jtvm->addLabel (std::to_string(this->m[0]), {cw, ty, 0.0f});
            ty -= th;
        }
        if (this->m[1] > 0) {
            jtvm->addLabel ("m_j", {0.0f, ty, 0.0f});
            jtvm->addLabel (std::to_string(this->m[1]), {cw, ty, 0.0f});
            ty -= th;
        }
        if (this->m[2] > 0) {
            jtvm->addLabel ("m_i", {0.0f, ty, 0.0f});
            jtvm->addLabel (std::to_string(this->m[2]), {cw, ty, 0.0f});
            ty -= th;
        }
        if (this->m[3] > 0) {
            jtvm->addLabel ("m_c", {0.0f, ty, 0.0f});
            jtvm->addLabel (std::to_string(this->m[3]), {cw, ty, 0.0f});
            ty -= th;
        }

        // r parameters
        if (this->m[1] > 0) {
            jtvm->addLabel ("r_j", {0.0f, ty, 0.0f});
            jtvm->addLabel (std::to_string(this->mconf->getFloat ("r_j", 0.0f)), {cw, ty, 0.0f});
            ty -= th;
        }
        if (this->m[2] > 0) {
            jtvm->addLabel ("r_i", {0.0f, ty, 0.0f});
            jtvm->addLabel (std::to_string(this->mconf->getFloat ("r_i", 0.0f)), {cw, ty, 0.0f});
            ty -= th;
        }
        if (this->m[3] > 0) {
            jtvm->addLabel ("r_c", {0.0f, ty, 0.0f});
            jtvm->addLabel (std::to_string(this->mconf->getFloat ("r_c", 0.0f)), {cw, ty, 0.0f});
            ty -= th;
        }

        v->addVisualModel (jtvm);
    }

    // 1x3 graphs (branches, centroids, selected)
    void graph_layout_b (const morph::Vector<float>& offset0)
    {
        morph::Vector<float> g_A = offset0 + morph::Vector<float>({0.0f, 0.0f, 0.0f});
        morph::Vector<float> g_B = offset0 + morph::Vector<float>({1.3f, 0.0f, 0.0f});
        morph::Vector<float> g_C = offset0 + morph::Vector<float>({2.6f, 0.0f, 0.0f});

        this->bv = new BranchVisual<T, N, B> (v->shaderprog, v->tshaderprog, g_A, &this->branches, &this->ax_history);
        this->bv->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->bv->target_scale.compute_autoscale (0, 1);
        this->bv->view = branchvisual_view::discs;
        this->bv->finalize();
        this->bv->addLabel ("Branches", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->bv, std::string("Tectal"));
        v->addVisualModel (this->bv);

        // Axon centroids: Centroids of branches viewed with a NetVisual
        this->cv = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_B, &this->ax_centroids);
        //this->cv->maxlen = this->conf->getDouble ("maxnetline", 1.0);
        //this->cv->viewmode = netvisual_viewmode::actual_nolines;
        this->cv->viewmode = netvisual_viewmode::actual;
        if (this->layout == graph_layout::b) {
            this->cv->radiusFixed = 0.02;
        }
        this->cv->finalize();
        this->cv->addLabel ("Axon centroids", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv, std::string("Tectal"));
        v->addVisualModel (this->cv);

        // Selected axons: This one gives an 'axon view'
        this->av = new BranchVisual<T, N, B> (v->shaderprog, v->tshaderprog, g_C, &this->branches, &this->ax_history);
        this->av->view = branchvisual_view::axonview;
        for (auto sa : this->seeaxons) { this->av->seeaxons.insert(sa); }
        this->av->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->av->target_scale.compute_autoscale (0, 1);
        this->av->finalize();
        this->av->addLabel ("Selected axons", {0.0f, 1.1f, 0.0f});
        v->addVisualModel (this->av);

        // Figure letters
        morph::Vector<float> ozero = {-0.2f, 1.1f, 0.0f};
        float lfs = 0.08f; // letter font size
        int lpts = 36; // letter point resolution
        morph::VisualModel* jtvm = new morph::VisualModel (v->shaderprog, v->tshaderprog, ozero);
        jtvm->addLabel ("A", g_A, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        jtvm->addLabel ("B", g_B, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        jtvm->addLabel ("C", g_C, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        this->v->addVisualModel (jtvm);
    }

    // Create a layout with diff. time end points. 2x4.
    void graph_layout_c (const morph::Vector<float>& offset0)
    {
        // This is a 2x4 layout
        morph::Vector<float> g_A = offset0 + morph::Vector<float>({0.0f, 0.0f, 0.0f});
        morph::Vector<float> g_B = offset0 + morph::Vector<float>({1.3f, 0.0f, 0.0f});
        morph::Vector<float> g_C = offset0 + morph::Vector<float>({2.6f, 0.0f, 0.0f});
        morph::Vector<float> g_D = offset0 + morph::Vector<float>({3.9f, 0.0f, 0.0f});

        morph::Vector<float> g_E = offset0 + morph::Vector<float>({0.0f, -1.6f, 0.0f});
        morph::Vector<float> g_F = offset0 + morph::Vector<float>({1.3f, -1.6f, 0.0f});
        morph::Vector<float> g_G = offset0 + morph::Vector<float>({2.6f, -1.6f, 0.0f});
        morph::Vector<float> g_H = offset0 + morph::Vector<float>({3.9f, -1.6f, 0.0f});

        // A Retinal cell positions
        this->v->addVisualModel (this->createTissueVisual (v->shaderprog, v->tshaderprog, g_A, ret, "Retinal", expression_view::cell_positions, 0, 2));

        // B Experiment: Another NetVisual view showing the target locations
        this->tcv = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_E, &this->ax_centroids);
        this->tcv->viewmode = netvisual_viewmode::targetplus;
        this->tcv->finalize();
        this->tcv->addLabel ("Experiment", {0.0f, 1.1f, 0.0f});
        this->v->addVisualModel (this->tcv);

        // t=0
        this->cv0 = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_B, &this->ax_centroids);
        this->cv0->maxlen = this->conf->getDouble ("maxnetline", 1.0);
        this->cv0->viewmode = netvisual_viewmode::actual;
        this->cv0->finalize();
        this->cv0->addLabel ("Axon centroids (t=0)", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv0, std::string("Tectal"));
        this->v->addVisualModel (this->cv0);

        // Axon centroids: Centroids of branches viewed with a NetVisual
        this->cv1 = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_C, &this->ax_centroids);
        this->cv1->maxlen = this->conf->getDouble ("maxnetline", 1.0);
        this->cv1->viewmode = netvisual_viewmode::actual;
        this->cv1->finalize();
        std::stringstream ss1;
        ss1 << "t=" << this->freeze_times[1];
        this->cv1->addLabel (ss1.str(), {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv1, std::string("Tectal"));
        this->v->addVisualModel (this->cv1);

        // Axon centroids, final positions
        this->cv = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_D, &this->ax_centroids);
        this->cv->maxlen = this->conf->getDouble ("maxnetline", 1.0);
        this->cv->viewmode = netvisual_viewmode::actual;
        this->cv->finalize();
        std::stringstream ssf;
        ssf << "t=" << this->conf->getInt("steps", -1);
        this->cv->addLabel (ssf.str(), {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv, std::string("Tectal"));
        this->v->addVisualModel (this->cv);

        // A graph of the SOS diffs between axon position centroids and target positions from retina
        // SEB working on this one
        this->gv = new morph::GraphVisual<T> (v->shaderprog, v->tshaderprog, g_F);
        this->gv->axisstyle = morph::axisstyle::twinax;
        this->gv->twodimensional = false;
        this->gv->setsize (0.9f, 1.0f);
        this->gv->setlimits (0, this->conf->getFloat ("steps", 1000),
                             0, this->conf->getFloat("graph_ymax", 1.0f),
                             0, this->conf->getFloat("graph_ymax2", 200.0f));
        this->gv->axislabelgap = 0.03f;
        this->gv->policy = morph::stylepolicy::lines;
        this->gv->ylabel = unicode::toUtf8 (unicode::epsilon);
        this->gv->ylabel2 = unicode::toUtf8 (unicode::eta);
        this->gv->xlabel = "t";
        this->gv->prepdata (unicode::toUtf8 (unicode::epsilon));
        this->gv->prepdata (unicode::toUtf8 (unicode::eta), morph::axisside::right);
        this->gv->finalize();
        this->v->addVisualModel (this->gv);

        // Branches: Visualise the branches with a custom VisualModel
        this->bv = new BranchVisual<T, N, B> (v->shaderprog, v->tshaderprog, g_G, &this->branches, &this->ax_history);
        this->bv->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->bv->target_scale.compute_autoscale (0, 1);
        this->bv->view = branchvisual_view::discs;
        this->bv->finalize();
        this->bv->addLabel ("Branches", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->bv, std::string("Tectal"));
        this->v->addVisualModel (this->bv);

        // C Selected axons: This one gives an 'axon view'
        this->av = new BranchVisual<T, N, B> (v->shaderprog, v->tshaderprog, g_H, &this->branches, &this->ax_history);
        this->av->view = branchvisual_view::axonview;
        for (auto sa : this->seeaxons) { this->av->seeaxons.insert(sa); }
        this->av->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->av->target_scale.compute_autoscale (0, 1);
        this->av->finalize();
        this->av->addLabel ("Selected axons", {0.0f, 1.1f, 0.0f});
        this->v->addVisualModel (this->av);

        // A 'text' only visual model to incorporate figure letters
        morph::Vector<float> ozero = {-0.2f, 1.1f, 0.0f};
        float lfs = 0.08f; // letter font size
        int lpts = 36; // letter point resolution
        morph::VisualModel* jtvm = new morph::VisualModel (v->shaderprog, v->tshaderprog, ozero);
        jtvm->addLabel ("A", g_A, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        jtvm->addLabel ("B", g_B, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        jtvm->addLabel ("C", g_C, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        jtvm->addLabel ("D", g_D, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        jtvm->addLabel ("E", g_E, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        jtvm->addLabel ("F", g_F, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        jtvm->addLabel ("G", g_G, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        jtvm->addLabel ("H", g_H, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        this->v->addVisualModel (jtvm);
    }

    void graph_layout_e (const morph::Vector<float>& offset0)
    {
        morph::Vector<float> offset = offset0;
        // Branches: Visualise the branches with a custom VisualModel
        this->bv = new BranchVisual<T, N, B> (v->shaderprog, v->tshaderprog, offset, &this->branches, &this->ax_history);
        this->bv->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->bv->target_scale.compute_autoscale (0, 1);
        this->bv->view = branchvisual_view::discs;
        this->bv->finalize();
        this->bv->addLabel ("Branches", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->bv, std::string("Tectal"));
        v->addVisualModel (this->bv);

        // Axon centroids: Centroids of branches viewed with a NetVisual
        offset[0] += 1.3f;
        this->cv = new NetVisual<T> (v->shaderprog, v->tshaderprog, offset, &this->ax_centroids);
        this->cv->viewmode = netvisual_viewmode::actual;
        if (this->layout == graph_layout::b) {
            this->cv->radiusFixed = 0.02;
        }
        this->cv->finalize();
        this->cv->addLabel ("Axon centroids", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv, std::string("Tectal"));
        v->addVisualModel (this->cv);

        offset[0] += 1.4f;

        // The position graph, like Brown/Reber and S&G papers
        this->gv = new morph::GraphVisual<T> (v->shaderprog, v->tshaderprog, offset);
        this->gv->twodimensional = false;
        this->gv->setlimits (0, 1, 0, 1);
        this->gv->policy = morph::stylepolicy::markers;
        this->gv->ylabel = "R " + unicode::toUtf8 (unicode::longrightarrow) + " tectum " + unicode::toUtf8 (unicode::longrightarrow) + " C";
        this->gv->xlabel = "N " + unicode::toUtf8 (unicode::longrightarrow) + " retina " + unicode::toUtf8 (unicode::longrightarrow) + " Y";
        v->addVisualModel (this->gv);
    }

    void graph_layout_f (const morph::Vector<float>& offset0)
    {
        morph::Vector<float> offset = offset0;
        // Axon centroids: Centroids of branches viewed with a NetVisual
        this->cv = new NetVisual<T> (v->shaderprog, v->tshaderprog, offset, &this->ax_centroids);
        this->cv->viewmode = netvisual_viewmode::actual;
        if (this->layout == graph_layout::b) {
            this->cv->radiusFixed = 0.02;
        }
        this->cv->finalize();
        this->cv->addLabel ("Axon centroids", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv, std::string("Tectal"));
        v->addVisualModel (this->cv);

        // Experiment: Another NetVisual view showing the target locations
        morph::Vector<float> expoff = {0.8f, -0.2f, 0.01f};
        offset += expoff;
        this->tcv = new NetVisual<T> (v->shaderprog, v->tshaderprog, offset, &this->ax_centroids);
        this->tcv->viewmode = netvisual_viewmode::targetplus;
        this->tcv->zoom = 0.5f;
        this->tcv->finalize();
        this->tcv->addLabel ("Experiment", {0.32f, 0.55f, 0.0f},
                             morph::colour::black, morph::VisualFont::Vera, 0.03, 24);
        v->addVisualModel (this->tcv);
        offset -= expoff;

        // A 'text' only visual model to display the sim time
        offset[0] += 1.05f;
        morph::VisualModel* jtvm = new morph::VisualModel (v->shaderprog, v->tshaderprog, offset);
        float ty = 1.1f; // text y position
        float th = 0.1f; // text height
        float l_x = 0.0f;
        float cw = 0.1f;
        jtvm->addLabel ("t", {l_x, ty, 0.0f});
        ty -= th;
        jtvm->addLabel ("t = ", {l_x, ty, 0.0f});
        jtvm->addLabel ("0", {l_x+cw, ty, 0.0f}, this->sim_time_txt);
        ty -= th;
        jtvm->addLabel ("sos: ", {l_x, ty, 0.0f});
        jtvm->addLabel ("0", {l_x, ty-0.8f*th, 0.0f}, this->emetric_txt);
        ty -= 2*th;
        jtvm->addLabel ("crossings: ", {l_x, ty, 0.0f});
        jtvm->addLabel ("0", {l_x, ty-0.8f*th, 0.0f}, this->crossings_txt);
        v->addVisualModel (jtvm);
    }

    // 1x4 graphs (expt, branches, centroids, selected)
    void graph_layout_g (const morph::Vector<float>& offset0, const std::string& startletter)
    {
        morph::Vector<float> g_A = offset0 + morph::Vector<float>({0.0f, 0.0f, 0.0f});
        morph::Vector<float> g_B = offset0 + morph::Vector<float>({1.5f, 0.0f, 0.0f});
        morph::Vector<float> g_C = offset0 + morph::Vector<float>({3.0f, 0.0f, 0.0f});
        morph::Vector<float> g_D = offset0 + morph::Vector<float>({4.5f, 0.0f, 0.0f});

        // Expt net
        this->tcv = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_A, &this->ax_centroids);
        this->tcv->viewmode = netvisual_viewmode::targetplus;
        this->tcv->finalize();
        this->tcv->addLabel ("Experiment", {0.0f, 1.1f, 0.0f});
        this->v->addVisualModel (this->tcv);

        // Branches
        this->bv = new BranchVisual<T, N, B> (v->shaderprog, v->tshaderprog, g_B, &this->branches, &this->ax_history);
        this->bv->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->bv->target_scale.compute_autoscale (0, 1);
        this->bv->view = branchvisual_view::discs;
        this->bv->finalize();
        this->bv->addLabel ("Branches", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->bv, std::string("Tectal"));
        this->v->addVisualModel (this->bv);

        // Axon centroids: Centroids of branches viewed with a NetVisual
        this->cv = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_C, &this->ax_centroids);
        this->cv->viewmode = netvisual_viewmode::actual;
        if (this->layout == graph_layout::b) {
            this->cv->radiusFixed = 0.02;
        }
        this->cv->finalize();
        this->cv->addLabel ("Axon centroids", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv, std::string("Tectal"));
        this->v->addVisualModel (this->cv);

        // A graph of the SOS diffs between axon position centroids and target positions from retina
        this->gv = new morph::GraphVisual<T> (v->shaderprog, v->tshaderprog, g_D);
        this->gv->twodimensional = false;
        this->gv->setlimits (0, this->conf->getFloat ("steps", 1000),
                             0, this->conf->getFloat("graph_ymax", 1.0f));
        this->gv->policy = morph::stylepolicy::lines;
        this->gv->ylabel = unicode::toUtf8(unicode::epsilon);
        this->gv->xlabel = "t";
        this->gv->prepdata ("All");
        this->gv->prepdata ("Surround");
        this->gv->prepdata ("Graft");
        this->gv->finalize();
        this->v->addVisualModel (this->gv);

        // A 'text' only visual model to display the RMS error
        morph::VisualModel* sosvm = new morph::VisualModel (v->shaderprog, v->tshaderprog, g_D);
        float ty = 0.87f; // text y position
        float l_x = 0.55f; // text x pos
        float th = 0.1f; // text height
        morph::TextGeometry lgeom = sosvm->addLabel (unicode::toUtf8(unicode::epsilon) + "[All]=", {l_x, ty, 0.0f});
        sosvm->addLabel ("0", {l_x + lgeom.width(), ty, 0.0f}, this->emetric_txt);
        ty -= 2*th;
        v->addVisualModel (sosvm);

        // Figure letters
        morph::Vector<float> ozero = {-0.2f, 1.1f, 0.0f};
        float lfs = 0.08f; // letter font size
        int lpts = 36; // letter point resolution
        char sl = 'A';
        if (!startletter.empty()) { sl = startletter[0]; }
        morph::VisualModel* jtvm = new morph::VisualModel (v->shaderprog, v->tshaderprog, ozero);
        jtvm->addLabel (std::string({sl}), g_A, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        this->v->addVisualModel (jtvm);
    }

    // 1x4 graphs (expt, branches, centroids, ret NT vs tec RC)
    void graph_layout_h (const morph::Vector<float>& offset0, const std::string& startletter)
    {
        morph::Vector<float> g_A = offset0 + morph::Vector<float>({0.0f, 0.0f, 0.0f});
        morph::Vector<float> g_B = offset0 + morph::Vector<float>({1.5f, 0.0f, 0.0f});
        morph::Vector<float> g_C = offset0 + morph::Vector<float>({3.0f, 0.0f, 0.0f});
        morph::Vector<float> g_D = offset0 + morph::Vector<float>({4.5f, 0.0f, 0.0f});

        this->tcv = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_A, &this->ax_centroids);
        this->tcv->viewmode = netvisual_viewmode::targetplus;
        this->tcv->finalize();
        this->tcv->addLabel ("Experiment", {0.0f, 1.1f, 0.0f});
        this->v->addVisualModel (this->tcv);

        this->bv = new BranchVisual<T, N, B> (v->shaderprog, v->tshaderprog, g_B, &this->branches, &this->ax_history);
        this->bv->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->bv->target_scale.compute_autoscale (0, 1);
        this->bv->view = branchvisual_view::discs;
        this->bv->finalize();
        this->bv->addLabel ("Branches", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->bv, std::string("Tectal"));
        this->v->addVisualModel (this->bv);

        // Axon centroids: Centroids of branches viewed with a NetVisual
        this->cv = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_C, &this->ax_centroids);
        this->cv->viewmode = netvisual_viewmode::actual;
        this->cv->finalize();
        this->cv->addLabel ("Axon centroids", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv, std::string("Tectal"));
        this->v->addVisualModel (this->cv);

        // NT vs RC graph here
        this->gv = new morph::GraphVisual<T> (v->shaderprog, v->tshaderprog, g_D);
        this->gv->twodimensional = false;
        this->gv->setlimits (0, 1, 0, 1);
        this->gv->policy = morph::stylepolicy::markers;
        this->gv->ylabel = "R " + unicode::toUtf8 (unicode::longrightarrow) + " tectum " + unicode::toUtf8 (unicode::longrightarrow) + " C";
        this->gv->xlabel = "N " + unicode::toUtf8 (unicode::longrightarrow) + " retina " + unicode::toUtf8 (unicode::longrightarrow) + " Y";
        this->gv->finalize();
        v->addVisualModel (this->gv);

        // A 'text' only visual model to display the RMS error
        morph::VisualModel* sosvm = new morph::VisualModel (v->shaderprog, v->tshaderprog, g_D);
        float ty = 0.87f; // text y position
        float l_x = 0.55f; // text x pos
        float th = 0.1f; // text height
        morph::TextGeometry lgeom = sosvm->addLabel (unicode::toUtf8(unicode::epsilon) + "=", {l_x, ty, 0.0f});
        sosvm->addLabel ("0", {l_x + lgeom.width(), ty, 0.0f}, this->emetric_txt);
        ty -= 2*th;
        v->addVisualModel (sosvm);

        // Figure letters
        morph::Vector<float> ozero = {-0.2f, 1.1f, 0.0f};
        float lfs = 0.08f; // letter font size
        int lpts = 36; // letter point resolution
        char sl = 'A';
        if (!startletter.empty()) { sl = startletter[0]; }
        morph::VisualModel* jtvm = new morph::VisualModel (v->shaderprog, v->tshaderprog, ozero);
        jtvm->addLabel (std::string({sl}), g_A, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        this->v->addVisualModel (jtvm);
    }

    void graph_layout_i (const morph::Vector<float>& offset0)
    {
        morph::Vector<float> g_A = offset0 + morph::Vector<float>({0.0f, 0.0f, 0.0f});
        morph::Vector<float> g_B = offset0 + morph::Vector<float>({1.5f, 0.0f, 0.0f});
        morph::Vector<float> g_C = offset0 + morph::Vector<float>({3.0f, 0.0f, 0.0f});
        morph::Vector<float> g_D = offset0 + morph::Vector<float>({4.5f, 0.0f, 0.0f});

        // Axon centroids: Centroids of branches viewed with a NetVisual
        this->cv1 = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_A, &this->ax_centroids);
        this->cv1->maxlen = this->conf->getDouble ("maxnetline", 1.0);
        this->cv1->viewmode = netvisual_viewmode::actual;
        this->cv1->finalize();
        std::stringstream ss1;
        ss1 << "t=" << this->freeze_times[1];
        this->cv1->addLabel (ss1.str(), {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv1, std::string("Tectal"));
        this->v->addVisualModel (this->cv1);

        // Axon centroids, final positions
        this->cv = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_B, &this->ax_centroids);
        this->cv->maxlen = this->conf->getDouble ("maxnetline", 1.0);
        this->cv->viewmode = netvisual_viewmode::actual;
        this->cv->finalize();
        std::stringstream ssf;
        ssf << "t=" << this->conf->getInt("steps", -1);
        this->cv->addLabel (ssf.str(), {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv, std::string("Tectal"));
        this->v->addVisualModel (this->cv);

        // Selected axons: This one gives an 'axon view'
        this->av = new BranchVisual<T, N, B> (v->shaderprog, v->tshaderprog, g_C, &this->branches, &this->ax_history);
        this->av->view = branchvisual_view::axonview;
        for (auto sa : this->seeaxons) { this->av->seeaxons.insert(sa); }
        this->av->rcpt_scale.compute_autoscale (rcpt_min, rcpt_max);
        this->av->target_scale.compute_autoscale (0, 1);
        this->av->finalize();
        this->av->addLabel ("Selected axons", {0.0f, 1.1f, 0.0f});
        this->v->addVisualModel (this->av);

        // Graph: A graph of the SOS diffs between axon position centroids and target positions from retina
        this->gv = new morph::GraphVisual<T> (v->shaderprog, v->tshaderprog, g_D);
        this->gv->axisstyle = morph::axisstyle::twinax;
        this->gv->twodimensional = false;
        this->gv->setsize (0.9f, 1.0f);
        this->gv->setlimits (0, this->conf->getFloat ("steps", 1000),
                             0, this->conf->getFloat("graph_ymax", 1.0f),
                             0, this->conf->getFloat("graph_ymax2", 200.0f));
        this->gv->axislabelgap = 0.03f;
        this->gv->policy = morph::stylepolicy::lines;
        this->gv->ylabel = unicode::toUtf8 (unicode::epsilon);
        this->gv->ylabel2 = unicode::toUtf8 (unicode::eta);
        this->gv->xlabel = "t";
        this->gv->prepdata (unicode::toUtf8 (unicode::epsilon));
        this->gv->prepdata (unicode::toUtf8 (unicode::eta), morph::axisside::right);
        this->gv->finalize();
        this->v->addVisualModel (this->gv);
    }

    // 1x2 graphs (centroids, selected)
    void graph_layout_j (const morph::Vector<float>& offset0, const std::string& startletter)
    {
        morph::Vector<float> g_A = offset0 + morph::Vector<float>({0.0f, 0.0f, 0.0f});
        morph::Vector<float> g_B = offset0 + morph::Vector<float>({1.5f, 0.0f, 0.0f});

        // Axon centroids: Centroids of branches viewed with a NetVisual
        this->cv = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_A, &this->ax_centroids);
        this->cv->viewmode = netvisual_viewmode::actual;
        if (this->layout == graph_layout::b) {
            this->cv->radiusFixed = 0.02;
        }
        this->cv->finalize();
        this->cv->addLabel ("Axon centroids", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv, std::string("Tectal"));
        this->v->addVisualModel (this->cv);

        // A graph of the SOS diffs between axon position centroids and target positions from retina
        this->gv = new morph::GraphVisual<T> (v->shaderprog, v->tshaderprog, g_B);
        this->gv->twodimensional = false;
        this->gv->setlimits (0, this->conf->getFloat ("steps", 1000),
                             0, this->conf->getFloat("graph_ymax", 1.0f));
        this->gv->policy = morph::stylepolicy::lines;
        this->gv->ylabel = unicode::toUtf8(unicode::epsilon);
        this->gv->xlabel = "t";
        this->gv->prepdata ("All");
        this->gv->prepdata ("Surround");
        this->gv->prepdata ("Graft");
        this->gv->finalize();
        this->v->addVisualModel (this->gv);

        // A 'text' only visual model to display the RMS error
        morph::VisualModel* sosvm = new morph::VisualModel (v->shaderprog, v->tshaderprog, g_B);
        float ty = 0.87f; // text y position
        float l_x = 0.55f; // text x pos
        float th = 0.1f; // text height
        morph::TextGeometry lgeom = sosvm->addLabel (unicode::toUtf8(unicode::epsilon) + "[All]=", {l_x, ty, 0.0f});
        sosvm->addLabel ("0", {l_x + lgeom.width(), ty, 0.0f}, this->emetric_txt);
        ty -= 2*th;
        v->addVisualModel (sosvm);

        // Figure letters
        morph::Vector<float> ozero = {-0.2f, 1.1f, 0.0f};
        float lfs = 0.08f; // letter font size
        int lpts = 36; // letter point resolution
        char sl = 'A';
        if (!startletter.empty()) { sl = startletter[0]; }
        morph::VisualModel* jtvm = new morph::VisualModel (v->shaderprog, v->tshaderprog, ozero);
        jtvm->addLabel (std::string({sl}), g_A, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        this->v->addVisualModel (jtvm);
    }

    // 1x4 graphs (expt, branches, centroids, ret NT vs tec RC)
    void graph_layout_k (const morph::Vector<float>& offset0, const std::string& startletter)
    {
        morph::Vector<float> g_A = offset0 + morph::Vector<float>({0.0f, 0.0f, 0.0f});
        morph::Vector<float> g_B = offset0 + morph::Vector<float>({1.5f, 0.0f, 0.0f});

        // Axon centroids: Centroids of branches viewed with a NetVisual
        this->cv = new NetVisual<T> (v->shaderprog, v->tshaderprog, g_A, &this->ax_centroids);
        this->cv->viewmode = netvisual_viewmode::actual;
        this->cv->finalize();
        this->cv->addLabel ("Axon centroids", {0.0f, 1.1f, 0.0f});
        this->addOrientationLabels (this->cv, std::string("Tectal"));
        this->v->addVisualModel (this->cv);

        // NT vs RC graph here
        this->gv = new morph::GraphVisual<T> (v->shaderprog, v->tshaderprog, g_B);
        this->gv->twodimensional = false;
        this->gv->setlimits (0, 1, 0, 1);
        this->gv->policy = morph::stylepolicy::markers;
        this->gv->ylabel = "R " + unicode::toUtf8 (unicode::longrightarrow) + " tectum " + unicode::toUtf8 (unicode::longrightarrow) + " C";
        this->gv->xlabel = "N " + unicode::toUtf8 (unicode::longrightarrow) + " retina " + unicode::toUtf8 (unicode::longrightarrow) + " Y";
        this->gv->finalize();
        v->addVisualModel (this->gv);

        // A 'text' only visual model to display the RMS error
        morph::VisualModel* sosvm = new morph::VisualModel (v->shaderprog, v->tshaderprog, g_B);
        float ty = 0.87f; // text y position
        float l_x = 0.55f; // text x pos
        float th = 0.1f; // text height
        morph::TextGeometry lgeom = sosvm->addLabel (unicode::toUtf8(unicode::epsilon) + "=", {l_x, ty, 0.0f});
        sosvm->addLabel ("0", {l_x + lgeom.width(), ty, 0.0f}, this->emetric_txt);
        ty -= 2*th;
        v->addVisualModel (sosvm);

        // Figure letters
        morph::Vector<float> ozero = {-0.2f, 1.1f, 0.0f};
        float lfs = 0.08f; // letter font size
        int lpts = 36; // letter point resolution
        char sl = 'A';
        if (!startletter.empty()) { sl = startletter[0]; }
        morph::VisualModel* jtvm = new morph::VisualModel (v->shaderprog, v->tshaderprog, ozero);
        jtvm->addLabel (std::string({sl}), g_A, morph::colour::black, morph::VisualFont::VeraBold, lfs, lpts);
        this->v->addVisualModel (jtvm);
    }

#endif // VISUALISE

    // The axons to see - these will have their path information stored
    std::set<size_t> seeaxons = {21, 38, 189, 378, 361};
    // Branches per axon
    unsigned int bpa = 8;
    // Number of RGCs on a side
    unsigned int rgcside = 20;
    // If true, then slow things down a bit in the visualization
    bool goslow = false;
    // Exit or keep showing graphics?
    bool immediate_exit = false;
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
    // Parameters vector (See Table 2 in the paper)
    //                        G        C       I        J        B
    morph::Vector<T, 5> m = { T{0.02}, T{0.2}, T{0.15}, T{0.15}, T{0.1}};
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
    // Has a genetic manipulation been applied?
    bool genetic_manipulation = false;
    // Path history is a map indexed by axon id. 3D as it's used for vis.
    std::map<size_t, morph::vVector<morph::Vector<T, 3>>> ax_history;
    // Receptor max and min - used across init() and visinit()
    T rcpt_max = -1e9;
    T rcpt_min = 1e9;
    // a title for this simulation
    std::string title = "";
#ifdef VISUALISE
    // A visual environment for the sim
    morph::Visual* v = nullptr;
    // A visual environment specifically for the tissue visualisation
    morph::Visual* tvv = nullptr;
    // Has visualisation been done already?
    bool visinit_done = false;
    // Specialised visualization of agents as spheres with a little extra colour patch on top
    BranchVisual<T, N, B>* bv = nullptr;
    // Another visualization to show axon paths with only a few axons
    BranchVisual<T, N, B>* av = nullptr;
    // Centroid visual
    NetVisual<T>* cv = nullptr;
    // Centroid visuals at different time points
    NetVisual<T>* cv0 = nullptr;
    NetVisual<T>* cv1 = nullptr;
    NetVisual<T>* cv2 = nullptr;
    NetVisual<T>* cv3 = nullptr;
    // Simulation times to stop updating graphs (see graph_layout::c)
    morph::Vector<size_t, 4> freeze_times;
    // How early to start showing the crossings count metric.
    unsigned int crosscount_from = 1000;
    // How often to show the crossings count metric.
    unsigned int crosscount_every = 50;
    // Centroid visual for targets
    NetVisual<T>* tcv = nullptr;
    // A graph for the SOS metric
    morph::GraphVisual<T>* gv = nullptr;
    // Make a couple of options for the graph layout for different figure types
    graph_layout layout = graph_layout::c;
#endif // VISUALISE
};
