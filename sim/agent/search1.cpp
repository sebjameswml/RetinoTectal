/*
 * Parameter searching with the RetinoTectal model
 *
 * Launch with a base json config for the model which can incorporate initial values for
 * parameters and with a second config which will specify the parameters to
 * optimise. Any experiment configs will be either generated by this code or loaded from
 * known locations.
 */

#include <string>
#include <vector>
#include <morph/tools.h>
#include <morph/Config.h>
#include <morph/NM_Simplex.h>
#define RANDSINGLE 1
#define NO_RANDDOUBLE 1
#include <morph/rng.h>
#include "agent1.h"

// Objective function involves running the agent based model
template <typename T=float, size_t N=4>
T objfn (Agent1<T, N, branch<T, N>>& model,
         morph::Config* mconf,
         std::vector<std::string>& params,
         std::vector<T>& param_values)
{
    // Set params in model
    model.reset();
    for (size_t i = 0; i < params.size(); ++i) {
        mconf->set (params[i], param_values[i]);
    }
    model.update_m();
    // Run it then get metrics
    model.run();
    AgentMetrics<T> m = model.get_metrics();
    T rtn = m.sos + m.crosscount;
    return rtn;
}

int main (int argc, char **argv)
{
    // Set up config objects
    std::string paramsfile_mdl("");
    std::string paramsfile_srch("");
    // The 'model id' and 'search id' are derived from their JSON files' names
    std::string m_id("");
    std::string s_id("");

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
            if (two.size() > 1) { m_id = two[1]; }
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
    std::vector<float> param_values;
    Json::Value params_j = sconf->getArray ("params");
    for (auto p : params_j) {
        params.push_back (p.asString());
        param_values.push_back (mconf->getFloat(p.asString(), 0));
    }

    std::cout << "Optimisable parameters:\n";
    size_t np = params.size(); // num params
    size_t nv = np+1;          // num vertices
    for (size_t i = 0; i < np; ++i) {
        std::cout << params[i] << " = " << param_values[i] << std::endl;
    }

    // NM Search parameters
    std::vector<std::vector<float>> i_vertices(nv);
    // First vertex set from config file
    i_vertices[0].resize(np);
    for (size_t i = 0; i < np; ++i) {
        i_vertices[0][i] = param_values[i];
    }
    // Other vertices obtained by adding/subtracting a bit (this could be a NM_Simplex method)
    for (size_t j = 1; j < nv; ++j) {
        i_vertices[j].resize(np);
        for (size_t i = 0; i < np; ++i) {
            float pv = param_values[i];
            pv = pv + (pv/2)*(morph::randSingle()-0.5f);
            i_vertices[j][i] = pv;
        }
    }

    morph::NM_Simplex<float> simp(i_vertices);
    simp.termination_threshold = sconf->getFloat ("nm_threshold", 0.001f);
    simp.too_many_operations = sconf->getUInt ("nm_toomany", 0); // 0 means go on for ever

    if (num_guiders == 4) {

        // for (each expt) {

        // Create/use an expt config
        //std::string paramsfile_expt = "./configs/a1/e_wt.json";
        std::string paramsfile_expt = "./configs/a1/e_retablate.json";
        morph::Config* econf = new morph::Config (paramsfile_expt);
        Agent1<float, 4, branch<float, 4>> model (econf, mconf);
        model.title = std::string("j4_") + m_id + std::string("_s_") + s_id;
        model.immediate_exit = true;
        model.randomly_seeded = false;

        AgentMetrics<float> metrics1;
        AgentMetrics<float> metrics2;
        // For each expt (or group of expts? many models?) run, then adjust params, etc:
        model.run();
        metrics1 = model.get_metrics();
        model.reset();

        // Condense metrics into a single value.
        float val = metrics1.sos + metrics1.crosscount;

        while (simp.state != morph::NM_Simplex_State::ReadyToStop) {

            std::cout << "During optimization:\n";
            std::vector<float> final_params = simp.best_vertex();
            if (params.size() != final_params.size()) { throw std::runtime_error ("Uh oh"); }
            for (size_t i = 0; i < np; ++i) {
                std::cout << params[i] << " = " << final_params[i] << std::endl;
            }

            if (simp.state == morph::NM_Simplex_State::NeedToComputeThenOrder) {
                // 1. apply objective to each vertex
                for (unsigned int i = 0; i <= simp.n; ++i) {
                    simp.values[i] = objfn (model, mconf, params, simp.vertices[i]);
                    std::cout << "Vertex " << i << " returns objective " << simp.values[i] << std::endl;
                }
                simp.order();

            } else if (simp.state == morph::NM_Simplex_State::NeedToOrder) {
                simp.order();

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeReflection) {
                val = objfn (model, mconf, params, simp.xr);
                std::cout << "Before apply_reflection, objective is " << val << std::endl;
                simp.apply_reflection (val);

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeExpansion) {
                val = objfn (model, mconf, params, simp.xe);
                std::cout << "Before apply_expansion, objective is " << val << std::endl;
                simp.apply_expansion (val);

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeContraction) {
                val = objfn (model, mconf, params, simp.xc);
                std::cout << "Before apply_contraction, objective is " << val << std::endl;
                simp.apply_contraction (val);
            }
        }

        std::cout << "After optimization:\n";
        std::vector<float> final_params = simp.best_vertex();
        if (params.size() != final_params.size()) { throw std::runtime_error ("Uh oh"); }
        for (size_t i = 0; i < np; ++i) {
            std::cout << params[i] << " = " << final_params[i] << std::endl;
        }

        delete econf;

    } else if (num_guiders == 2) {
        std::cerr << "Not implemented for num_guiders == 2\n";
    }

    delete sconf;
    delete mconf;
    return 0;
}
