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

// Fig 4 displays wildtype results
//static constexpr experiment expt = experiment::wildtype;
static constexpr experiment expt = experiment::knockin_hetero;
//static constexpr experiment expt = experiment::knockin_homo;

static constexpr int tissue_n = 100;

// Do you want a 3D graph?
static constexpr bool threeD = false;

int main (int argc, char** argv)
{
    ktt1d<expt, tissue_n> model;
    morph::Config conf;
    if (argc > 1) {
        std::string conf_path(argv[1]);
        std::cout << "Setting parameters from " << conf_path << std::endl;
        conf.init (conf_path);
        if (conf.ready) {
            model.comp_param_A = conf.get<float> ("A", model.comp_param_A);
            model.comp_param_B = conf.get<float> ("B", model.comp_param_B);
            model.comp_param_D = conf.get<float> ("D", model.comp_param_D);

            model.chem_param_alpha = conf.get<float> ("alpha", model.chem_param_alpha);
            model.set_act_param_a (conf.get<float> ("a", model.get_act_param_a()));
            model.set_act_param_b (conf.get<float> ("b", model.get_act_param_b()));
            model.set_act_param_gamma (conf.get<float> ("gamma", model.get_act_param_gamma()));
        }
    }

    morph::Visual v(940, 420, "Tsigankov, Koulakov, Trippett");
    v.setSceneTrans (morph::vec<float,3>({-0.413069f, -0.461113f, -2.6f}));

    // Dataset style for bargraphs
    morph::DatasetStyle dsb(morph::stylepolicy::bar); // Draw a bar graph by creating a bar policy DatasetStyle
    dsb.markercolour = morph::colour::blue; // bar colour
    dsb.markersize = 0.01f;
    dsb.showlines = false;

    // A bar graph for ligand expression
    auto gv = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({-0.8,0.65,0}));
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
    gv->setdata (sc_cr_axis, model.col_la, dsb);
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
    gv2->setdata (ret_nt_axis, model.ret_ra, dsb);
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
    cvis->finalize();
    v.addVisualModel (cvis);

    // Next job is to visualize the synapses
    morph::vec<float, 2> dx = { 0.01f, 0.01f };
    morph::Grid<int, float> grid1 (tissue_n, tissue_n, dx, {0,0},
                                   morph::GridDomainWrap::None,
                                   morph::GridOrder::bottomleft_to_topright_colmaj); // orients output
    auto gridv1 = std::make_unique<morph::GridVisual<float, int, float>>(&grid1, morph::vec<float>({0.8,0,0}));
    v.bindmodel (gridv1);
    gridv1->setScalarData (&model.ret_synapse_density);
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
    gridv1->finalize();
    auto gridv1p = v.addVisualModel (gridv1);

    constexpr int visevery = 1000;
    int loop = 0;
    while (!v.readyToFinish /*&& loop < 20*/) {
        model.step();
        // Update graphs every 'visevery' model steps
        if (loop % visevery == 0) {
            if (loop % 10000 == 0) {
                std::cout << "Loop " << loop << ". n_syn is now " << model.n_syn << "\n";
            }
            v.waitevents (0.001);
            model.compute_ret_synapse_density();
            gridv1p->reinit();
            v.render();
            if (loop > 1000000) { break; } // 10^6 iterations to stationary soln
        }
        loop += 1;
    }

    std::cout << "Stationary? Finished after " << (loop-1) << " iterations\n";
    v.keepOpen();

    return 0;
}
