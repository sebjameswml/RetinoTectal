/*
 * Evaluate the agent across several manipulations.
 */

// Uniform random number generator, for branches
#include "branch_rng_float.h"

#include <string>
#include <vector>
#include <morph/tools.h>
#include <morph/Config.h>
#include "agent1.h"

int main (int argc, char **argv)
{
    // Set up config objects
    std::string paramsfile_mdl("");
    // The 'model id' is derived from the JSON file name
    std::string m_id("");

    int runtimesteps = -1;
    if (argc >= 2) {
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
        // A 2nd arg is converted into a number to be simulation steps
        if (argc > 2) { runtimesteps = std::stoi (argv[2]); }
    } else {
        std::cerr << "Usage: " << argv[0] << " <model.json>\n";
        return -1;
    }

    morph::Config* mconf = new morph::Config(paramsfile_mdl);
    if (!mconf->ready) {
        std::cerr << "Failed to read model config " << paramsfile_mdl << ". Exiting.\n";
        delete mconf;
        return 1;
    }

    morph::Tools::createDirIf ("./log/agent");
    size_t num_guiders = mconf->getInt("num_guiders", 4);
    if (num_guiders != 4) { delete mconf; return -4; }

    morph::Config* econf = (morph::Config*)0;
    std::vector<std::string> e_manips;
    std::vector<AgentMetrics<float>> e_manips_results;
    // Add all the manipulation experiments to evaluate here
    e_manips.push_back ("configs/a1/e_retablate.json");
    e_manips.push_back ("configs/a1/e_tecrot90.json");
    e_manips.push_back ("configs/a1/e_tecrot180.json");
    e_manips.push_back ("configs/a1/e_tecablate.json");
    e_manips.push_back ("configs/a1/e_tecswap.json");

    for (unsigned int i = 0; i < e_manips.size(); ++i) {

        econf = new morph::Config(e_manips[i]);
        std::string e_id("");
        std::string fname = e_manips[i];
        morph::Tools::stripUnixPath (fname);
        std::vector<std::string> one = morph::Tools::stringToVector (fname, ".json");
        if (!one.empty()) {
            std::vector<std::string> two = morph::Tools::stringToVector (one[0], "e_");
            if (two.size() > 1) { e_id = two[1]; }
        }
        econf->set ("exit", true); // Exit at end of sim
        econf->set ("graph_layout", 5); // Choose a layout to view
        // Update simulation steps if user provided a number
        if (runtimesteps > 0) { econf->set ("steps", runtimesteps); }
        std::string outfile = std::string("./log/agent/") + m_id + std::string("_") + e_id + std::string(".h5");
        Agent1<float, 4, branch<float, 4>> model (econf, mconf);
        model.imagedir = std::string ("./log/agent1_eval");
        model.title = std::string("j4_") + m_id + std::string("_") + e_id;
        model.run();
        //model.save (outfile);
        AgentMetrics<float> am = model.get_metrics();
        am.id = e_id;
        std::cout << "For manipulation " << am.id << ", SOS: " << am.sos << " and crosscount = " << am.crosscount << std::endl;
        e_manips_results.push_back (am);
        delete econf;
    }

    std::cout << "------------------------------------------\n";
    std::cout << "Results for model " << m_id << std::endl;
    std::cout << "------------------------------------------\n";
    for (auto am : e_manips_results) {
        std::cout << "For manipulation " << am.id << ", SOS: " << am.sos << " and crosscount = " << am.crosscount << std::endl;
    }
    std::cout << "------------------------------------------\n";

    delete mconf;
    return 0;
}
