#include "ktt.h"

// Visualization includes
#include <morph/vec.h>
#include <morph/Visual.h>
#include <morph/GraphVisual.h>
#include <morph/Grid.h>
#include <morph/GridVisual.h>
#include <morph/unicode.h>
#include <morph/Config.h>
#include <morph/ConfigVisual.h>

int main (int argc, char** argv)
{
    // constexpr experiment expt = experiment::wildtype;
    constexpr experiment expt = experiment::knockin_hetero;
    // constexpr experiment expt = experiment::knockin_homo;
    constexpr int tissue_n = 100;
    // Do you want a 3D graph world?
    constexpr bool threeD = false;
    // How many models? (Vary a param across these)
    constexpr int num_models = 10;

    std::map<int, std::unique_ptr<ktt1d<expt, tissue_n>>> models;

    morph::Config conf;
    if (argc > 1) {
        std::string conf_path(argv[1]);
        std::cout << "Setting parameters from " << conf_path << std::endl;
        conf.init (conf_path);
    }

    morph::vvec<float> gammas;
    gammas.linspace (0.03, 0.1, num_models);

    for (int mi = 0; mi < num_models; ++mi) {
        models[mi] = std::make_unique<ktt1d<expt, tissue_n>>();

        if (conf.ready) {
            models[mi]->comp_param_A = conf.get<float> ("A", models[mi]->comp_param_A);
            models[mi]->comp_param_B = conf.get<float> ("B", models[mi]->comp_param_B);
            models[mi]->comp_param_D = conf.get<float> ("D", models[mi]->comp_param_D);

            models[mi]->chem_param_alpha = conf.get<float> ("alpha", models[mi]->chem_param_alpha);
            models[mi]->set_act_param_a (conf.get<float> ("a", models[mi]->get_act_param_a()));
            models[mi]->set_act_param_b (conf.get<float> ("b", models[mi]->get_act_param_b()));
            models[mi]->set_act_param_gamma (conf.get<float> ("gamma", models[mi]->get_act_param_gamma()));
        }

        // Hard code a range of gammas. Replace with any other param you want to examine
        models[mi]->set_act_param_gamma (gammas[mi]);
    }

    morph::Visual v(940, 420, "Tsigankov, Koulakov, Trippett");
    v.setSceneTrans (morph::vec<float,3>{-0.413069f, -0.461113f, -2.6f});

    // Dataset style for bargraphs
    morph::DatasetStyle dsb(morph::stylepolicy::bar); // Draw a bar graph by creating a bar policy DatasetStyle
    dsb.markercolour = morph::colour::blue; // bar colour
    dsb.markersize = 0.01f;
    dsb.showlines = false;

    // A bar graph for ligand expression
    auto gv = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>{ -0.8f, 0.65f, 0.0f });
    v.bindmodel (gv);
    morph::vvec<float> ret_nt_axis; // N-T index from 0 to 100
    ret_nt_axis.arange (0.0f, static_cast<float>(tissue_n), 1.0f);
    morph::vvec<float> sc_cr_axis = ret_nt_axis; // CR index from 0 to 100
    // Plots in the paper have and R-C axis for termination site
    morph::vvec<float> sc_rc_axis;
    sc_rc_axis.arange (static_cast<float>(tissue_n-1), -1.0f, -1.0f);

    gv->setsize(1.0f, 0.35f);
    gv->scalingpolicy_y = morph::scalingpolicy::manual_min;
    gv->datamin_y = 0;
    gv->setdataaxisdist (0.01f + dsb.markersize/2.0f);
    gv->setlimits (0.0f, static_cast<float>(tissue_n-1), 0.0f, 1.5f);
    gv->setdata (sc_cr_axis, models[0]->col_la, dsb);
    gv->xlabel = "SC posn (Caudal -> Rostral)";
    gv->ylabel = "Ligand";
    gv->finalize();
    v.addVisualModel (gv);

    // A bar graph for rgc receptor expression
    auto gv2 = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({-0.8,0,0}));
    v.bindmodel (gv2);
    gv2->setsize(1.0f, 0.35f);
    gv2->scalingpolicy_y = morph::scalingpolicy::manual_min;
    gv2->datamin_y = 0;
    gv2->setdataaxisdist (0.01f + dsb.markersize/2.0f);
    gv2->setlimits (0.0f, static_cast<float>(tissue_n-1), 0.0f, 1.5f);
    dsb.markercolour = morph::colour::lime;
    gv2->setdata (ret_nt_axis, models[0]->ret_ra, dsb);
    gv2->xlabel = "Ret. posn (Nasal -> Temporal)";
    gv2->ylabel = "Receptor";
    gv2->finalize();
    v.addVisualModel (gv2);

    morph::DatasetStyle ds(morph::stylepolicy::markers);
    ds.markercolour = morph::colour::crimson;

    // Visualize stuff from Config
    std::vector<std::string> confkeys = { "alpha", "A", "B", "D", "a", "b", "gamma" };
    auto cvis = std::make_unique<morph::ConfigVisual<>> (&conf, confkeys,
                                                         morph::vec<float>({0.3,1.0,0}),
                                                         morph::TextFeatures(0.05));
    v.bindmodel (cvis);
    cvis->twodimensional = !threeD;
    cvis->finalize();
    v.addVisualModel (cvis);

    // Next job is to visualize the synapses
    morph::vec<float, 2> dx = { 0.01f, 0.01f };
    morph::Grid<int, float> grid1 (tissue_n, tissue_n, dx, {0,0},
                                   morph::GridDomainWrap::None,
                                   morph::GridOrder::bottomleft_to_topright_colmaj); // orients output

    // num_models graphs - one per model
    std::map<int, morph::GridVisual<float, int, float>*> gvptrs;
    for (int mi = 0; mi < num_models; ++mi) {
        morph::vec<float> gvloc = {0.8f + 1.2f * static_cast<float>(mi), 0.0f, 0.0f};
        auto gridv1 = std::make_unique<morph::GridVisual<float, int, float>>(&grid1, gvloc);
        v.bindmodel (gridv1);
        gridv1->setScalarData (&models[mi]->ret_synapse_density);
        if constexpr (threeD == true) {
            gridv1->gridVisMode = morph::GridVisMode::Columns;
            gridv1->interpolate_colour_sides = true;
            gridv1->twodimensional = false;
            gridv1->zScale.compute_autoscale (0.0f, 1.0f);
        } else {
            gridv1->twodimensional = true;
            gridv1->gridVisMode = morph::GridVisMode::RectInterp;
            gridv1->zScale.setParams(0,0);
        }
        gridv1->cm.setType (morph::ColourMapType::Twilight);
        gridv1->colourScale.compute_autoscale (0.0f, 1.0f);
        gridv1->addLabel ("Ret", morph::vec<float>({-0.18f, 0.49f, 0.0f}), morph::TextFeatures(0.09f));
        gridv1->addLabel ("SC index", morph::vec<float>({0.3f, -0.16f, 0.0f}), morph::TextFeatures(0.09f));
        std::string gammalbl = morph::unicode::toUtf8(morph::unicode::gamma) + "=" + std::to_string (models[mi]->get_act_param_gamma());
        gridv1->addLabel (gammalbl, morph::vec<float>({0.3f, -0.3f, 0.0f}), morph::TextFeatures(0.06f));
        gridv1->finalize();
        gvptrs[mi] = v.addVisualModel (gridv1);
    }

    constexpr int visevery = 1000;
    int loop = 0;
    morph::vvec<unsigned long long int> n_syns (num_models, 0);
    while (!v.readyToFinish /*&& loop < 20*/) {

#pragma omp parallel for
        for (int mi = 0; mi < num_models; ++mi) {
            for (int li = 0; li < visevery; ++li) { models[mi]->step(); }
        }
        loop += visevery;

        // Update graphs every 'visevery' model steps
        if (loop % 10000 == 0) {
            for (int mi = 0; mi < num_models; ++mi) { n_syns[mi] = models[mi]->n_syn; }
            std::cout << "Loop " << loop << ". n_syns: " << n_syns << "\n";
        }
        v.waitevents (0.001);
        for (int mi = 0; mi < num_models; ++mi) {
            models[mi]->compute_ret_synapse_density();
            gvptrs[mi]->reinit();
        }
        v.render();
        if (loop > 1000000) { break; } // 10^6 iterations to stationary soln

    }

    std::cout << "Stationary? Finished after " << (loop-1) << " iterations\n";
    v.keepOpen();

    return 0;
}
