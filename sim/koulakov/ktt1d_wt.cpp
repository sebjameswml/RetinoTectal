#include "ktt.h"

// Visualization includes
#include <morph/vec.h>
#include <morph/Visual.h>
#include <morph/GraphVisual.h>
#include <morph/Grid.h>
#include <morph/GridVisual.h>
#include <morph/unicode.h>

// Fig 4 displays wildtype results
static constexpr experiment expt = experiment::wildtype;

static constexpr int tissue_n = 100;

int main()
{
    ktt1d<expt, tissue_n> model; // Keep model default alpha of 30
    ktt1d<expt, tissue_n> model_bigalph;
    model_bigalph.alpha = 100000.0f;

    morph::Visual v(1024, 768, "Tsigankov and Koulakov, A unifying model 2006");
    v.setSceneTrans (morph::vec<float,3>({-0.492555f, 0.367545f, -5.6f}));

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

    // Next job is to visualize the synapses
    morph::vec<float, 2> dx = { 0.01f, 0.01f };
    morph::Grid<int, float> grid1 (tissue_n, tissue_n, dx);
    auto gridv1 = std::make_unique<morph::GridVisual<float, int, float>>(&grid1, morph::vec<float>({0.8,0,0}));
    v.bindmodel (gridv1);
    gridv1->setScalarData (&model.ret_synapse_density);
    gridv1->twodimensional = true;
    gridv1->cm.setType (morph::ColourMapType::Jet);
    gridv1->colourScale.compute_autoscale (0.0f, 1.0f);
    gridv1->zScale.compute_autoscale (0.0f, 0.0f);
    gridv1->finalize();
    auto gridv1p = v.addVisualModel (gridv1);

    int loop = 0;
    while (!v.readyToFinish /*&& loop < 20*/) {
        model.step();
        model_bigalph.step();
        // Update graphs every 1000 model steps
        if (loop++ % 1000 == 0) {
            std::cout << "1000 loops. n_syn is now " << model.n_syn << "\n";
            v.waitevents (0.018);
            model.compute_ret_synapse_density();
            gridv1p->reinit();
            v.render();
            if (loop > 1000000) { break; } // 10^6 iterations to stationary soln
        }
    }

    std::cout << "Stationary? Finished with loop=" << loop << std::endl;
    v.keepOpen();

    return 0;
}
