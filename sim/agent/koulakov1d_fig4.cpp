#include "koulakov1d.h"

// Visualization includes
#include <morph/vec.h>
#include <morph/Visual.h>
#include <morph/GraphVisual.h>
#include <morph/unicode.h>

// Fig 4 displays wildtype results
static constexpr experiment expt = experiment::wildtype;

int main()
{
    k1d<expt> model; // Keep model default alpha of 30
    k1d<expt> model_bigalph;
    model_bigalph.alpha = 100000.0f;

    morph::Visual v(1024, 768, "Koulakov and Tsigankov Fig. 4");
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
    ret_nt_axis.arange (0.0f, static_cast<float>(N), 1.0f);
    morph::vvec<float> sc_cr_axis = ret_nt_axis; // CR index from 0 to 100
    // Plots in the paper have and R-C axis for termination site
    morph::vvec<float> sc_rc_axis;
    sc_rc_axis.arange (static_cast<float>(N-1), -1.0f, -1.0f);

    gv->setsize(1.0f, 0.35f);
    gv->scalingpolicy_y = morph::scalingpolicy::manual_min;
    gv->datamin_y = 0;
    gv->setdataaxisdist (0.01f + dsb.markersize/2.0f);
    gv->setlimits (0.0f, static_cast<float>(N-1), 0.0f, 1.5f);
    gv->setdata (sc_cr_axis, model.la, dsb);
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
    gv2->setlimits (0.0f, static_cast<float>(N-1), 0.0f, 1.5f);
    dsb.markercolour = morph::colour::lime;
    gv2->setdata (ret_nt_axis, model.ra, dsb);
    gv2->xlabel = "Ret. posn (Nasal -> Temporal)";
    gv2->ylabel = "Receptor";
    gv2->finalize();
    v.addVisualModel (gv2);

    morph::DatasetStyle ds(morph::stylepolicy::markers);
    ds.markercolour = morph::colour::crimson;

    // Initial ordering (random)
    auto gv3 = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({0.8,0,0}));
    v.bindmodel (gv3);
    gv3->setlimits (0.0f, static_cast<float>(N), 0.0f, static_cast<float>(N));
    ds.datalabel = "Initial state";
    gv3->setdata (model.rgc_for_sc_idx.as_float(), sc_rc_axis, ds);
    gv3->xlabel = "Ret. posn  (Nasal -> Temporal)";
    gv3->ylabel = "SC term. (Rostral -> Caudal)";
    gv3->finalize();
    v.addVisualModel (gv3);

    // The main model (alpha 30)
    auto gv4 = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({-0.8,-1.6,0}));
    v.bindmodel (gv4);
    gv4->twodimensional = false;
    gv4->setlimits (0.0f, static_cast<float>(N), 0.0f, static_cast<float>(N));
    ds.datalabel = morph::unicode::toUtf8 (morph::unicode::alpha) + std::string(" = ") + std::to_string (model.alpha);
    gv4->setdata (model.rgc_for_sc_idx.as_float(), sc_rc_axis, ds);
    gv4->xlabel = "Ret. posn (Nasal -> Temporal)";
    gv4->ylabel = "SC term. (Rostral -> Caudal)";
    gv4->finalize();
    auto gv4p = v.addVisualModel (gv4);

    // The high alpha model (alpha 100000)
    auto gv5 = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({0.8,-1.6,0}));
    v.bindmodel (gv5);
    gv5->setlimits (0.0f, static_cast<float>(N), 0.0f, static_cast<float>(N));
    ds.datalabel = morph::unicode::toUtf8 (morph::unicode::alpha) + std::string(" = ") + std::to_string (model_bigalph.alpha);
    gv5->setdata (model_bigalph.rgc_for_sc_idx.as_float(), sc_cr_axis, ds);
    gv5->xlabel = "Ret. posn (Nasal -> Temporal)";
    gv5->ylabel = "SC term. (Rostral -> Caudal)";
    gv5->finalize();
    auto gv5p = v.addVisualModel (gv5);

    int loop = 0;
    while (!v.readyToFinish) {
        model.step();
        model_bigalph.step();
        // Update graphs every 1000 model steps
        if (loop++ % 1000 == 0) {
            v.waitevents (0.018);
            gv4p->update (model.rgc_for_sc_idx.as_float(), sc_rc_axis, 0);
            gv5p->update (model_bigalph.rgc_for_sc_idx.as_float(), sc_rc_axis, 0);
            v.render();
            if (loop > 100000) { break; }
        }
    }

    v.keepOpen();

    return 0;
}
