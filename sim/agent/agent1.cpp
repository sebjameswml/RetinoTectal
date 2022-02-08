/*
 * Retinotectal model resembling one presented by Hugh Simpson and Geoffrey
 * Goodhill in "A simple model can unify a broad range of phenomena in retinotectal map
 * development", Biol Cybern (2011) 104:9-29
 *
 * I'm bringing the idea of variable interaction with signalling gradients and
 * competition to try to get rid of the non-biological part of Simpson & Goodhill's
 * work.
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
    std::string paramsfile("");
    std::string paramsfile_mdl("");
    // The 'model id' and 'expt id' are derived from their JSON files' *names*
    std::string m_id("");
    std::string e_id("");

    if (argc >= 2) {
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

    // Cool idea! Use morph::Config to process cmd line args in standard way to override parameters. Nice!
    mconf->process_args (argc, argv);
    conf->process_args (argc, argv);

    morph::Tools::createDirIf ("./log/agent");

    // Create an additional file identifier
    std::stringstream coss;
    for (auto co : conf->config_overrides) {
        coss <<  "_" << co.first << "_" << co.second;
    }

    std::string outfile = std::string("./log/agent/") + m_id
    + std::string("_") + e_id + coss.str() + std::string(".h5");

    std::string branch_model = mconf->getString ("branch_model", "james_agent");

    size_t num_guiders = mconf->getInt("num_guiders", 4);
    if (num_guiders == 4) {
        if (branch_model == "gebhardt") {
            Agent1<float, 4, branch_geb<float, 4>> model (conf, mconf);
            model.title = std::string("geb_") + m_id + std::string("_") + e_id;
            model.run();
            model.am.save (outfile);
        } else if (branch_model == "stochastic") {
            Agent1<float, 4, branch_stochastic<float, 4, 11>> model (conf, mconf);
            model.title = std::string("stoc_") + m_id + std::string("_") + e_id;
            model.run();
            model.am.save (outfile);
        } else {
            Agent1<float, 4, branch<float, 4>> model (conf, mconf);
            model.title = std::string("j4_") + m_id + std::string("_") + e_id;
            model.run();
            model.am.save (outfile);
        }
    } else if (num_guiders == 2) {
        Agent1<float, 2, branch<float, 2>> model (conf, mconf);
        model.title = std::string("j2_") + m_id + std::string("_") + e_id;
        model.run();
        model.am.save (outfile);
    }

    //std::cout << "conf:\n" << conf->str() << std::endl;

    delete conf;
    if (argc > 2) { /* mconf->write ("mconf.json"); */ delete mconf; }
    return 0;
}
