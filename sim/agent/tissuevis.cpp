/*
 * Just visualise the tissue of agent1.cpp (without running the model)
 */

#include "branch_rng_float.h"

#include <string>
#include <vector>
#include <morph/tools.h>
#include <morph/Config.h>
#include "agent1.h"

int main (int argc, char **argv)
{
    // Set up config objects
    std::string paramsfile("");
    std::string paramsfile_mdl("");
    // The 'model id' and 'expt id' are derived from their JSON files' *names*
    std::string m_id("");
    std::string e_id("");

    if (argc >= 3) {
        // If given two arguments, then the first is the model config and the second is the expt/sim config
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

        paramsfile = std::string(argv[2]);
        fname = paramsfile; // e.g.: e_wt.json
        morph::Tools::stripUnixPath (fname);
        // Find the thing between "m_" and ".json"
        one = morph::Tools::stringToVector (fname, ".json");
        if (!one.empty()) {
            std::vector<std::string> two = morph::Tools::stringToVector (one[0], "e_");
            if (two.size() > 1) { e_id = two[1]; }
        }

    } else if (argc == 2) {
        // With one argument, we use the same file for both model and expt/sim config
        paramsfile = std::string(argv[1]);
    } else {
        // Create an empty/default json file to run the sim with default values
        paramsfile = "./a1.json";
        morph::Tools::copyStringToFile ("{}\n", paramsfile);
    }

    morph::Config* conf = (morph::Config*)0;
    morph::Config* mconf = (morph::Config*)0;

    bool need_exit = false;
    std::cout << "Opening params config files " << paramsfile << " and model config " << paramsfile_mdl << std::endl;
    conf = new morph::Config(paramsfile);
    if (argc < 3) {
        mconf = conf; // with < two arguments, we use one file for both configs
    } else {
        mconf = new morph::Config(paramsfile_mdl);
    }

    if (!conf->ready) {
        std::cerr << "Failed to read sim/expt config " << paramsfile << ". Exiting.\n";
        need_exit = true;
    }
    if (argc > 2 && !mconf->ready) {
        std::cerr << "Failed to read model config " << paramsfile_mdl << ". Exiting.\n";
        need_exit = true;
    }

    if (need_exit) {
        delete conf;
        if (argc > 2) { delete mconf; }
        return 1;
    }

    mconf->process_args (argc, argv);
    conf->process_args (argc, argv);

    morph::Tools::createDirIf ("./log/agent");

    std::string outfile = std::string("./log/agent/") + m_id
    + std::string("_") + e_id + std::string(".h5");

    std::string branch_model = mconf->getString ("branch_model", "james_agent");

    size_t num_guiders = mconf->getInt("num_guiders", 4);
    if (num_guiders == 4) {
        if (branch_model == "gebhardt") {
            Agent1<float, 4, branch_geb<float, 4>> model (conf, mconf);
            model.showtissue();
        } else {
            std::cout << "branch_model: " << branch_model << std::endl;
            Agent1<float, 4, branch<float, 4>> model (conf, mconf);
            model.showtissue();
        }
    } else if (num_guiders == 2) {
        Agent1<float, 2, branch<float, 2>> model (conf, mconf);
        model.showtissue();
    }

    //std::cout << "conf:\n" << conf->str() << std::endl;

    delete conf;
    if (argc > 2) { delete mconf; }
    return 0;
}
