// Visualization includes
#include <morph/vec.h>
#include <morph/Visual.h>
#include <morph/GraphVisual.h>
#include <morph/CartGrid.h>
#include <morph/CartGridVisual.h>
#include <morph/unicode.h>

// For pointers to manipulate graphs
struct g_ptrs {
    morph::CartGridVisual<float>* g_prob_dist = nullptr;
    morph::GraphVisual<float>* g_prob_dist_maxima = nullptr;
    morph::GraphVisual<float>* g_big_alpha = nullptr;
};

// Nice for debug to be able to swing the graphs around, but otherwise prefer 2d
static constexpr bool two_dee = false;

// This figure has 3 similar columns
template<experiment E>
g_ptrs plot_col (morph::Visual<>& v, morph::vec<float> offset,
                 const k1d<E>& model, const k1d<E>& model_bigalpha,
                 morph::CartGrid* cg, morph::vvec<float>* prob_data)
{
    g_ptrs ptrs;
    float graph_step = 1.2f;

    // Dataset style for bargraphs
    morph::DatasetStyle dsb(morph::stylepolicy::bar);
    dsb.markercolour = morph::colour::blue;
    dsb.markersize = 0.01f;
    dsb.showlines = false;

    // A bar graph for ligand expression
    auto gv = std::make_unique<morph::GraphVisual<float>> (offset);
    v.bindmodel (gv);
    morph::vvec<float> ret_nt_axis; // N-T index from 0 to 100
    ret_nt_axis.arange (0.0f, static_cast<float>(N), 1.0f);
    morph::vvec<float> sc_cr_axis = ret_nt_axis; // CR index from 0 to 100
    // Plots in the paper have and R-C axis for termination site
    morph::vvec<float> sc_rc_axis;
    sc_rc_axis.arange (static_cast<float>(N-1), -1.0f, -1.0f);

    gv->twodimensional = two_dee;
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

    offset[1] -= graph_step/2.0f;

    // A bar graph for rgc receptor expression
    auto gv2 = std::make_unique<morph::GraphVisual<float>> (offset);
    v.bindmodel (gv2);
    gv2->twodimensional = two_dee;
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

    offset[1] -= graph_step;

    // Prob density here
    auto cgv = std::make_unique<morph::CartGridVisual<float>> (cg, offset);
    v.bindmodel (cgv);
    cgv->twodimensional = two_dee;
    cgv->cartVisMode = morph::CartVisMode::RectInterp;
    cgv->setScalarData (prob_data);
    cgv->cm.setType (morph::ColourMapType::Magma); // MonovalRed closest to paper
    if constexpr (two_dee == true) { cgv->zScale.setParams (0,0); }
    cgv->finalize();
    ptrs.g_prob_dist = v.addVisualModel (cgv);

    offset[1] -= graph_step;

    // Prob density maxima
    auto gv4 = std::make_unique<morph::GraphVisual<float>> (offset);
    v.bindmodel (gv4);
    gv4->twodimensional = two_dee;
    gv4->setlimits (0.0f, static_cast<float>(N), 0.0f, static_cast<float>(N));
    gv4->xlabel = "Ret. posn (Nasal -> Temporal)";
    gv4->ylabel = "SC term. (Rostral -> Caudal)";
    gv4->finalize();
    ptrs.g_prob_dist_maxima = v.addVisualModel (gv4);

    offset[1] -= 1.1f * graph_step;

    // Big alpha
    // The high alpha model (alpha 100000)
    auto gv5 = std::make_unique<morph::GraphVisual<float>> (offset);
    v.bindmodel (gv5);
    gv5->twodimensional = two_dee;
    gv5->setlimits (0.0f, static_cast<float>(N), 0.0f, static_cast<float>(N));
    morph::DatasetStyle ds(morph::stylepolicy::markers);
    ds.markercolour = morph::colour::crimson;
    ds.datalabel = morph::unicode::toUtf8 (morph::unicode::alpha) + std::string(" = ") + std::to_string (model_bigalpha.alpha);
    gv5->setdata (model_bigalpha.rgc_for_sc_idx.as_float(), sc_cr_axis, ds);
    gv5->xlabel = "Ret. posn (Nasal -> Temporal)";
    gv5->ylabel = "SC term. (Rostral -> Caudal)";
    gv5->finalize();
    ptrs.g_big_alpha = v.addVisualModel (gv5);

    return ptrs;
}
