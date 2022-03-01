/*
 * Parameter searching with the RetinoTectal model
 *
 * Launch with a base json config for the model which can incorporate initial values for
 * parameters and with a second config which will specify the parameters to
 * optimise. Any experiment configs will be either generated by this code or loaded from
 * known locations.
 */

// Uniform random number generator, for branches
#include "branch_rng_float.h"

#include <string>
#include <vector>
#include <morph/tools.h>
#include <morph/Config.h>
using nlohmann::json;
#include <morph/Anneal.h>
#include "agent1.h"
#include <csignal>

// Do we visualise the optimization itself with graphs?
#ifdef OPTVIS
# include <sstream>
# include <morph/Visual.h>
# include <morph/GraphVisual.h>
# include <morph/ScatterVisual.h>
# include <morph/TriaxesVisual.h>
// Global access to the visual object
morph::Visual* pv = nullptr;
#endif

// A count of the number of sims used in objfn() for text output
unsigned int model_sim_count = 0;
// To allow the signal handler to close the optimiser, make it a global pointer.
morph::Anneal<double>* optimiser = nullptr;
// The name of the model and search files are added to "anneal1" for saving data out.
// The 'model id' and 'search id' are derived from their JSON files' names
std::string s_id("");
std::string m_id("");
// If true, output a json config representing every single parameter set computed
static constexpr bool debug_json_config_on_every_step = false;

// Objective function involves running the agent based model
template <typename T=float, size_t N=4>
T objfn (Agent1<T, N, branch<T, N>>& model1,
         morph::Config* mconf,
         const std::vector<std::string>& params,
         const morph::vVector<T>& param_values)
{
    // Set params in model(s)
    for (size_t i = 0; i < params.size(); ++i) {
        mconf->set (params[i], param_values[i]);
    }
    model1.update_m();
    model1.reset(); // Updates initial positions based on r_i, r_j etc

    // Run model and then get metrics
    model1.run();
    ++model_sim_count;
    AgentMetrics<T> m1m = model1.get_metrics();

    // Here's a combination of the sos differences between the expected map and the actual map, plus a cross count.
    std::cout << "wt expt metrics: " << m1m.str() << " for parameters: ";
    for (size_t i = 0; i < params.size(); ++i) {
        std::cout << params[i] << " = " << param_values[i] << ", ";
    }
    std::cout << " (Sim count: " << model_sim_count << ")\n";

    // Compute the metric.
    T scaled_nc = T{1} + T{0.01} * m1m.crosscount.back();
    T rtn = m1m.sos.back() * scaled_nc;

    return rtn;
}

// A signal handler for optimiser.
void signalHandler (int signum)
{
    std::cout << "Optimisation interrupted. Saving data...\n";
    if (optimiser != nullptr) {
        morph::Tools::createDirIf ("log/anneal1");
        std::string fname = "log/anneal1/anneal1_" + m_id + std::string("_")
        + morph::Tools::filenameTimestamp() + std::string(".h5");
        std::cout << " to file "<< fname << std::endl;
        optimiser->save (fname);
        std::cout << "...saved. Best objective so far was " << optimiser->f_x_best
                  << " for params " << optimiser->x_best << std::endl;
    }
#ifdef OPTVIS
    std::cout << "Saved data. Leaving optimisation visualisation window open...\n";
    if (pv != nullptr) { pv->keepOpen(); }
#endif
    exit (signum);
}

#ifdef OPTVIS
void tav_setup (morph::TriaxesVisual<float>* tav, const size_t start_idx, const size_t dimensions,
                const morph::vVector<float>& param_range_min, const morph::vVector<float>& param_range_max,
                const std::vector<std::string>& params)
{
    size_t _0 = start_idx % dimensions;
    size_t _1 = (start_idx+1) % dimensions;
    size_t _2 = (start_idx+2) % dimensions;
    std::cout << "index _0:" << _0 << ", _1:" << _1 << ", _2:" << _2 << std::endl;
    tav->input_min = {param_range_min[_0], param_range_min[_1], param_range_min[_2]};
    tav->input_max = {param_range_max[_0], param_range_max[_1], param_range_max[_2]};
    std::cout << "tav->input_min: " << tav->input_min << "; tav->input_max: " << tav->input_max << std::endl;
    tav->xlabel = params[_0];
    tav->ylabel = params[_1];
    tav->zlabel = params[_2];
    tav->finalize();
}

// Specialise AnnealVisual to have an extra keybind and a 'start index' attribute
class AnnealVisual : public morph::Visual
{
public:
    AnnealVisual (int width, int height, const std::string& title) : morph::Visual (width, height, title) {}
    //! Dimension view start index
    size_t start_idx = 0;
    size_t dimensions = 1;
protected:
    //! Act on keys and toggle 'start index' for the relevant VisualModels
    virtual void key_callback_extra (GLFWwindow* window, int key, int scancode, int action, int mods)
    {
        if (key == GLFW_KEY_D && action == GLFW_PRESS) {
            std::cout << "Changing start index from " << start_idx << " to " << (start_idx + 1) << std::endl;
            start_idx = (start_idx+1) % dimensions;
        }

        if (key == GLFW_KEY_H && action == GLFW_PRESS) {
            std::cout << "AnnealVisual help:\n";
            std::cout << "d: cycle start index\n";
        }
    }
};
#endif

// E.g.: pbm && ./build/sim/agent/search1c configs/a1/m_eE_GCI.json configs/a1/s_GCI.json
int main (int argc, char **argv)
{
    signal(SIGINT, signalHandler);
    signal(SIGTERM, signalHandler);

    // Set up config objects
    std::string paramsfile_mdl("");
    std::string paramsfile_srch("");

    static constexpr size_t expected_args = 2;
    if (argc == (expected_args+1)) {
        // Two arguments. The first is the model config and the second is the search config
        paramsfile_mdl = std::string(argv[1]);
        std::string fname = paramsfile_mdl; // e.g.: m_el.json
        morph::Tools::stripUnixPath (fname);
        // Find the thing between "m_" and ".json"
        std::vector<std::string> one = morph::Tools::stringToVector (fname, ".json");
        if (!one.empty()) {
            std::vector<std::string> two = morph::Tools::stringToVector (one[0], "m_");
            // two should have size 2, with the first entry an empty string
            if (two.size() > 1) {
                m_id = two[1];
            } else {
                std::cout << "unknown model - was probably not in form 'm_ID.json'.\n";
                m_id = "unknown";
            }
            one.clear();
        }

        paramsfile_srch = std::string(argv[2]);
        fname = paramsfile_srch; // e.g.: s_GJ.json or s_GI.json
        morph::Tools::stripUnixPath (fname);
        // Find the thing between "s_" and ".json"
        one = morph::Tools::stringToVector (fname, ".json");
        if (!one.empty()) {
            std::vector<std::string> two = morph::Tools::stringToVector (one[0], "s_");
            if (two.size() > 1) { s_id = two[1]; }
        }
    } else {
        std::cerr << "Require " << expected_args << " args\n";
        return 1;
    }

    morph::Config* sconf = new morph::Config (paramsfile_srch);
    morph::Config* mconf = new morph::Config (paramsfile_mdl);

    // Check configs
    bool need_exit = false;
    if (!sconf->ready) {
        std::cerr << "Failed to read search(i.e. optimization) config " << paramsfile_srch << ". Exiting.\n";
        need_exit = true;
    }
    if (!mconf->ready) {
        std::cerr << "Failed to read model config " << paramsfile_mdl << ". Exiting.\n";
        need_exit = true;
    }
    if (need_exit) {
        delete sconf;
        delete mconf;
        return 1;
    }

    morph::Tools::createDirIf ("./log/agent");
    std::string base_path = std::string("./log/agent/") + m_id + std::string("_") + s_id;
    std::string outfile = base_path + std::string(".h5");
    std::string outfile_js = base_path + std::string(".json");
    std::string branch_model = mconf->getString ("branch_model", "james_agent");
    size_t num_guiders = mconf->getInt("num_guiders", 4);

    // Open s_*.json and get the array "params".
    std::vector<std::string> params;
    morph::vVector<double> param_values;
    morph::vVector<morph::Vector<double,2>> param_ranges;
    json params_j = sconf->get("params");
    for (auto p : params_j) {
        params.push_back (p.get<std::string>());
        param_values.push_back (mconf->getFloat(p.get<std::string>(), 0));
        // Set ranges for the params
        if (params.back() == "r_c") {
            json r_c_range = sconf->get("r_c_range");
            // Fixme: Check r_c_range contains 2 entries
            param_ranges.push_back (morph::Vector<double, 2>{r_c_range[0].get<double>(), r_c_range[1].get<double>()});
        } else if (params.back() == "r_i") {
            json r_i_range = sconf->get("r_i_range");
            param_ranges.push_back (morph::Vector<double, 2>{r_i_range[0].get<double>(), r_i_range[1].get<double>()});
        } else if (params.back() == "r_j") {
            json r_j_range = sconf->get("r_j_range");
            param_ranges.push_back (morph::Vector<double, 2>{r_j_range[0].get<double>(), r_j_range[1].get<double>()});
        } else if (params.back() == "s") {
            json s_range = sconf->get("s_range");
            param_ranges.push_back (morph::Vector<double, 2>{s_range[0].get<double>(), s_range[1].get<double>()});
        } else if (params.back() == "m_g") {
            json m_g_range = sconf->get("m_g_range");
            param_ranges.push_back (morph::Vector<double, 2>{m_g_range[0].get<double>(), m_g_range[1].get<double>()});
        } else if (params.back() == "m_c") {
            json m_c_range = sconf->get("m_c_range");
            param_ranges.push_back (morph::Vector<double, 2>{m_c_range[0].get<double>(), m_c_range[1].get<double>()});
        } else if (params.back() == "m_i") {
            json m_i_range = sconf->get("m_i_range");
            param_ranges.push_back (morph::Vector<double, 2>{m_i_range[0].get<double>(), m_i_range[1].get<double>()});
        } else if (params.back() == "m_j") {
            json m_j_range = sconf->get("m_j_range");
            param_ranges.push_back (morph::Vector<double, 2>{m_j_range[0].get<double>(), m_j_range[1].get<double>()});
        }
    }
    morph::vVector<float> param_range_max(param_ranges.size(), float{0});
    morph::vVector<float> param_range_min(param_ranges.size(), float{0});
    morph::vVector<float> one_over_param_maxes(param_ranges.size(), float{0});
    for (size_t i = 0; i < param_ranges.size(); ++i) {
        one_over_param_maxes[i] = static_cast<float>(1.0/param_ranges[i][1]);
        param_range_min[i] = param_ranges[i][0];
        param_range_max[i] = param_ranges[i][1];
    }
    morph::vVector<float> param_range_diff = param_range_max - param_range_min;
    morph::vVector<float> param_range_offs = param_range_min / param_range_diff;

    optimiser = new morph::Anneal<double>(param_values, param_ranges);
    // Anneal ASA params from sconf:
    optimiser->temperature_ratio_scale = sconf->getDouble ("temperature_ratio_scale", 1e-2);
    optimiser->temperature_anneal_scale = sconf->getDouble ("temperature_anneal_scale", 100.0);
    optimiser->cost_parameter_scale_ratio = sconf->getDouble ("cost_parameter_scale_ratio", 1.5);
    optimiser->acc_gen_reanneal_ratio = sconf->getDouble ("acc_gen_reanneal_ratio", 1e-6);
    optimiser->objective_repeat_precision = sconf->getDouble ("objective_repeat_precision", 1e-1);
    optimiser->f_x_best_repeat_max = sconf->getUInt ("f_x_best_repeat_max", 15);
    optimiser->reanneal_after_steps = sconf->getUInt ("reanneal_after_steps", 100);
    for (auto pn : params) {
        optimiser->param_names.push_back (pn);
    }
    optimiser->init();
    std::cout << "optimiser's objective_repeat_precision is " << optimiser->objective_repeat_precision << std::endl;

#ifdef OPTVIS
    // Set up the visualisation
    AnnealVisual v (1920, 1080, "Optimization progress");
    v.dimensions = param_range_max.size();
    v.zNear = 0.001;
    v.setSceneTransZ (-3.0f);
    v.lightingEffects (true);
    pv = &v;
    morph::Vector<float, 3> offset = { -0.7, 0.0, 0.0 };

    // First a scatter plot that can be updated. Just using a ScatterVisual for this.
    size_t sv_start_idx_last = v.start_idx;
    morph::ScatterVisual<double>* sv = new morph::ScatterVisual<double> (v.shaderprog, offset);
    sv->radiusFixed = 0.002f;
    sv->colourScale.compute_autoscale (0, 30);
    sv->cm.setType (morph::ColourMapType::Plasma);
    sv->finalize();
    v.addVisualModel (sv);

    morph::TriaxesVisual<float>* tav = new morph::TriaxesVisual<float> (v.shaderprog, v.tshaderprog, offset);
    tav_setup (tav, v.start_idx, v.dimensions, param_range_min, param_range_max, params);
    v.addVisualModel (tav);

    offset[0] += 2.0f;
    // Add a graph to track T_i and T_cost
    morph::GraphVisual<double>* graph1 = new morph::GraphVisual<double> (v.shaderprog, v.tshaderprog, offset);
    graph1->twodimensional = false;
    graph1->setlimits (0, 1000, -10, 1);
    graph1->policy = morph::stylepolicy::lines;
    graph1->ylabel = "log(T)";
    graph1->xlabel = "Anneal time";
    graph1->prepdata ("Tparam");
    graph1->prepdata ("Tcost");
    graph1->finalize();
    v.addVisualModel (graph1);

    offset[0] += 1.4f;
    morph::GraphVisual<double>* graph2 = new morph::GraphVisual<double> (v.shaderprog, v.tshaderprog, offset);
    graph2->twodimensional = false;
    graph2->setlimits (0, 1000, 0.0f, 2.0f);
    graph2->policy = morph::stylepolicy::lines;
    graph2->ylabel = "obj value";
    graph2->xlabel = "Anneal time";
    graph2->prepdata ("f_x");
    graph2->prepdata ("f_x_best");
    graph2->finalize();
    v.addVisualModel (graph2);

    offset[0] += 1.4f;
    morph::GraphVisual<double>* graph3 = new morph::GraphVisual<double> (v.shaderprog, v.tshaderprog, offset);
    graph3->twodimensional = false;
    graph3->setlimits (0, 1000, -1.0f, 100.0f);
    graph3->policy = morph::stylepolicy::lines;
    graph3->ylabel = "obj value";
    graph3->xlabel = "Anneal time";
    graph3->prepdata ("f_x_cand");
    graph3->finalize();
    v.addVisualModel (graph3);

    // Text labels to show additional information that might update
    morph::Vector<float> lpos = {-0.08f, 0.03f, 0.0f};
    morph::VisualTextModel* fps_tm;
    v.addLabel ("Unset", lpos, fps_tm);

    // Fixed text labels
    std::stringstream ss;
    lpos[1] -= 0.02f;
    ss << "reanneal_after_steps = " << optimiser->reanneal_after_steps;
    v.addLabel (ss.str(), lpos);

    lpos[1] -= 0.02f; ss.str("");
    ss << "temperature_ratio_scale = " << optimiser->temperature_ratio_scale;
    v.addLabel (ss.str(), lpos);

    lpos[1] -= 0.02f; ss.str("");
    ss << "temperature_anneal_scale = " << optimiser->temperature_anneal_scale;
    v.addLabel (ss.str(), lpos);

    lpos[1] -= 0.02f; ss.str("");
    ss << "cost_parameter_scale_ratio = " << optimiser->cost_parameter_scale_ratio;
    v.addLabel (ss.str(), lpos);

    lpos[1] -= 0.02f; ss.str("");
    ss << "acc_gen_reanneal_ratio = " << optimiser->acc_gen_reanneal_ratio;
    v.addLabel (ss.str(), lpos);

    lpos[1] -= 0.02f; ss.str("");
    ss << "objective_repeat_precision = " << optimiser->objective_repeat_precision;
    v.addLabel (ss.str(), lpos);

    lpos[1] -= 0.02f; ss.str("");
    ss << "f_x_best_repeat_max = " << optimiser->f_x_best_repeat_max;
    v.addLabel (ss.str(), lpos);

#endif // OPTVIS

    if (num_guiders == 4) {

        // for (each expt) {

        // Create/use an expt config
        std::string paramsfile_expt = "./optimization/e_wt.json";
        morph::Config* econf1 = new morph::Config (paramsfile_expt);
        Agent1<float, 4, branch<float, 4>> model1 (econf1, mconf);
        model1.title = std::string("j4_") + m_id + std::string("_1_s_") + s_id;
        model1.immediate_exit = true;
        model1.randomly_seeded = false;

        // ASA optimisation
#ifdef OPTVIS
        unsigned int max_repeats_so_far = 0;
#endif
        while (optimiser->state != morph::Anneal_State::ReadyToStop) {
            if (optimiser->state == morph::Anneal_State::NeedToCompute) {
                morph::vVector<float> xc = optimiser->x_cand.as_float();
                optimiser->f_x_cand = objfn (model1, mconf, params, xc);
                if constexpr (debug_json_config_on_every_step ==  true) {
                    // Save model params and objective into json for analysis
                    mconf->set ("f_x", optimiser->f_x_cand);
                    std::string jspath = "log/anneal1/" + m_id + std::string("_")
                    + std::to_string(optimiser->steps) + std::string(".json");
                    mconf->write (jspath);
                }

            } else if (optimiser->state == morph::Anneal_State::NeedToComputeSet) {
                optimiser->f_x_plusdelta = objfn (model1, mconf, params, optimiser->x_plusdelta.as_float());
            } else {
                throw std::runtime_error ("Unexpected state for anneal object.");
            }

#ifdef OPTVIS
            // Add to optimization visualisation, which could be a scatter plot.
            // Append to the 2D graph of sums:
            double tkmean = optimiser->T_k.mean();
            double tcostmean = optimiser->T_cost.mean();
            graph1->append ((float)optimiser->steps, std::log(tkmean), 0);
            graph1->append ((float)optimiser->steps, std::log(tcostmean), 1);
            graph2->append ((float)optimiser->steps, optimiser->f_x, 0);
            graph2->append ((float)optimiser->steps, optimiser->f_x_best, 1);
            graph3->append ((float)optimiser->steps, optimiser->f_x_cand, 0);

            // Add parameter set to the scattervisual. Scale coords by param_maxes
            morph::vVector<float> x_cand_dbl (optimiser->x_cand.size());
            std::copy (optimiser->x_cand.begin(), optimiser->x_cand.end(), x_cand_dbl.begin());
            morph::vVector<float> coord = (x_cand_dbl / param_range_diff) - param_range_offs;
            // Deal with < 3 coords
            size_t sz = coord.size();
            if (sz < 3) {
                coord.resize (3);
                for (size_t i = sz; i < 3; ++i) { coord[i] = float{0}; }
            }
            static constexpr float scatter_max_sz = 0.05f; // How big should the bad blobs be allowed to get?
            static constexpr float blob_divisor = 600.0f; // For scaling the size of the scatter blobs

            // If there are >3 coords, then have a 'start index' and show 3 of
            // them. Allow user to swich the 'start index' so they can show dimensions
            // 0,1,2 or 1,2,3, or 2,3,4 etc, looping back to 0 as necessary. If start
            // changes, then Scatter needs to be re-constructed.
            if (v.start_idx != sv_start_idx_last) {
                std::cout << "Rebuilding the scatter plot!\n";
                // rebuild the scatter, start by clearing it
                sv->clear();
                // Make a single container of params from the Anneal object
                morph::vVector<morph::vVector<double>> param_hist = optimiser->param_hist_rejected;
                morph::vVector<float> f_param_hist = optimiser->f_param_hist_rejected.as_float();
                std::cout << "param_hist from rejected size is " << param_hist.size() << std::endl;
                param_hist.insert (param_hist.end(), optimiser->param_hist_accepted.begin(), optimiser->param_hist_accepted.end());
                f_param_hist.insert (f_param_hist.end(), optimiser->f_param_hist_accepted.begin(), optimiser->f_param_hist_accepted.end());
                std::cout << "param_hist with accepted size is " << param_hist.size() << std::endl;
                // Now add all the points
                for (size_t i = 0; i < param_hist.size(); ++i) {
                    morph::vVector<float> c = (param_hist[i].as_float() / param_range_diff) - param_range_offs;
                    std::cout << "coordinate c = " << c << " has obj f value " << f_param_hist[i] << std::endl;
                    sv->add ({c[v.start_idx%v.dimensions], c[(v.start_idx+1)%v.dimensions], c[(v.start_idx+2)%v.dimensions]},
                             f_param_hist[i], (f_param_hist[i]/blob_divisor > scatter_max_sz ? scatter_max_sz : f_param_hist[i]/blob_divisor));
                }
                // Change the Triaxes visual
                tav->clear();
                tav->idx = 0;
                tav_setup (tav, v.start_idx, v.dimensions, param_range_min, param_range_max, params);
                sv_start_idx_last = v.start_idx;
            }
            sv->add ({coord[v.start_idx%v.dimensions], coord[(v.start_idx+1)%v.dimensions], coord[(v.start_idx+2)%v.dimensions]},
                     optimiser->f_x_cand, (optimiser->f_x_cand/blob_divisor > scatter_max_sz ? scatter_max_sz : optimiser->f_x_cand/blob_divisor));

            max_repeats_so_far = max_repeats_so_far < optimiser->f_x_best_repeats ? optimiser->f_x_best_repeats : max_repeats_so_far;
            std::stringstream ss;
            ss << "T_k = " << tkmean << ", T_cost = " << tcostmean
               << ", current f_x_best = " << optimiser->f_x_best
               << ", f_x_best_repeats = " << optimiser->f_x_best_repeats
               << " (highest num repeats so far: " << max_repeats_so_far << ")";
            fps_tm->setupText (ss.str());

            glfwWaitEventsTimeout (0.0166);
            v.render();
#endif
            // Every 100 steps save out data from optimiser? Or do it at and and catch
            // TERM signal and save data before exit? Probably that.
            std::cout << "For params (in order): ";
            bool _first = true;
            for (auto p : params) {
                if (!_first) {
                    std::cout << ", ";
                } else { _first = false; }
                std::cout << p;
            }
            std::cout << "\nCurrent x_best: " << optimiser->x_best
                      << " (objective: " << optimiser->f_x_best << ")\n";
            optimiser->step();
        }

        std::cout << "After optimization (simulated " << model_sim_count << " times):\n";

        std::cout << "Optimisation finished. Saving data...\n";
        morph::Tools::createDirIf ("log/anneal1");
        std::string ts = morph::Tools::filenameTimestamp();
        std::string fname = "log/anneal1/anneal1_" + m_id + std::string("_") + ts + std::string(".h5");
        std::cout << " to file "<< fname << std::endl;
        optimiser->save (fname);
        std::cout << "...saved. Best objective was " << optimiser->f_x_best
                  << " for params " << optimiser->x_best << std::endl;
        // Store origin json configs used for this anneal run
        {
            std::string fname_conf = "log/anneal1/anneal1_" + m_id + std::string("_") + ts + "_mconf.json";
            mconf->write (fname_conf);
            fname_conf = "log/anneal1/anneal1_" + m_id + std::string("_") + ts + "_sconf.json";
            sconf->write (fname_conf);
        }

        morph::vVector<double> final_params = optimiser->x_best;
        if (params.size() != final_params.size()) { throw std::runtime_error ("Uh oh"); }
        for (size_t i = 0; i < params.size(); ++i) {
            std::cout << params[i] << " = " << final_params[i] << std::endl;
        }

        std::string mdl_conf_out = "log/anneal1/anneal1_" + m_id + std::string("_") + ts + "_mopt.json";
        mconf->write (mdl_conf_out);
        if (mconf->ready == true) {
            std::cout << "Wrote optimised model to file " << mdl_conf_out << std::endl;
        } else {
            std::cerr << "Failed to write optimised model to file " << mdl_conf_out << std::endl;
        }
        delete econf1;

    } else if (num_guiders == 2) {
        std::cerr << "Not implemented for num_guiders == 2\n";
    }

    delete sconf;
    delete mconf;
    delete optimiser;
    return 0;
}
