#include "koulakov1d.h"

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

// Nice for debug to be able to swing the graphs around
bool two_dee = false;

// This figure has 3 similar columns
template<experiment E>
g_ptrs plot_col (morph::Visual& v, morph::vec<float> offset,
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

// Fig 5 displays wildtype and knockin results
int main()
{
    k1d<experiment::wildtype> model_wt;
    k1d<experiment::knockin_hetero> model_het;
    k1d<experiment::knockin_homo> model_hom;

    // Big alpha versions
    k1d<experiment::wildtype> model_wt_biga;
    k1d<experiment::knockin_hetero> model_het_biga;
    k1d<experiment::knockin_homo> model_hom_biga;
    model_wt_biga.alpha = 100000.0f;
    model_het_biga.alpha = model_wt_biga.alpha;
    model_hom_biga.alpha = model_wt_biga.alpha;

    morph::Visual v(1024, 1024, "Koulakov and Tsigankov Fig. 5");

    // Make a CartGrid
    float pix = 0.01f;
    auto cg = std::make_unique<morph::CartGrid> (pix, pix, 0.0f, 0.0f,
                                                 N*pix-pix,
                                                 N*pix-pix, 0.0f,
                                                 morph::CartDomainShape::Rectangle,
                                                 morph::CartDomainWrap::Horizontal);
    cg->setBoundaryOnOuterEdge();

    morph::vvec<float> pd_wt (N*N, 0.0f);
    morph::vvec<float> pd_het (N*N, 0.0f);
    morph::vvec<float> pd_hom (N*N, 0.0f);
    g_ptrs graphs_wt = plot_col (v, morph::vec<float>({0,0,0}), model_wt, model_wt_biga, cg.get(), &pd_wt);
    g_ptrs graphs_het = plot_col (v, morph::vec<float>({1.5,0,0}), model_het, model_het_biga, cg.get(), &pd_het);
    g_ptrs graphs_hom = plot_col (v, morph::vec<float>({3.0,0,0}), model_hom, model_hom_biga, cg.get(), &pd_hom);

    // The SC R-C axis for plot updates
    morph::vvec<float> sc_rc_axis;
    sc_rc_axis.arange (static_cast<float>(N-1), -1.0f, -1.0f);


    // Prepare biga graphs, which doesn't take long
    int loop = 0;
    while (!v.readyToFinish) {
        model_wt.step();
        model_wt_biga.step();
        model_het.step();
        model_het_biga.step();
        model_hom.step();
        model_hom_biga.step();
        // Update graphs every 1000 model steps
        if (loop++ % 1000 == 0) {
            v.waitevents (0.018);
            graphs_wt.g_big_alpha->update (model_wt_biga.rgc_for_sc_idx.as_float(), sc_rc_axis, 0);
            graphs_het.g_big_alpha->update (model_het_biga.rgc_for_sc_idx.as_float(), sc_rc_axis, 0);
            graphs_hom.g_big_alpha->update (model_hom_biga.rgc_for_sc_idx.as_float(), sc_rc_axis, 0);
            v.render();
            if (loop > 50000) { break; }
        }
    }

    std::cout << "Acquire stationary solution...\n";
    for (; loop < 1000000; ++loop) {
        model_wt.step();
        model_het.step();
        model_hom.step();
    }

    v.waitevents (0.018);
    v.render();

    std::cout << "Accumulate probability distributions...\n";
    // Now we have stationary solutions, prepare the distributions
    for (int i = 0; i < 50000; ++i) {
        // Run 1000 steps 50000 times...
        for (int j = 0; j < 1000; ++j) {
            model_wt.step();
            model_het.step();
            model_hom.step();
        }
        // Accumulate positions from models
        for (int sci = 0; sci < N; ++sci) {
            // retinal termination is model_wt.rgc_for_sc_idx[sci]
            int dataidx = (N-1-sci) * N + model_wt.rgc_for_sc_idx[sci];
            pd_wt[dataidx] += 1.0f;
            dataidx = (N-1-sci) * N + model_het.rgc_for_sc_idx[sci];
            pd_het[dataidx] += 1.0f;
            dataidx = (N-1-sci) * N + model_hom.rgc_for_sc_idx[sci];
            pd_hom[dataidx] += 1.0f;
        }
    }
    pd_wt /= 50000.0f;
    pd_het /= 50000.0f;
    pd_hom /= 50000.0f;

    graphs_wt.g_prob_dist->colourScale.reset();
    graphs_het.g_prob_dist->colourScale.reset();
    graphs_hom.g_prob_dist->colourScale.reset();

    graphs_wt.g_prob_dist->updateData (&pd_wt); // etc
    graphs_het.g_prob_dist->updateData (&pd_het);
    graphs_hom.g_prob_dist->updateData (&pd_hom);

    // Determine maxima
    morph::vvec<float> wt_x;
    wt_x.arange (0.0f, static_cast<float>(N-1), 2.0f);
    morph::vvec<float> ki_x;
    ki_x.arange (1.0f, static_cast<float>(N), 2.0f);

    morph::vvec<float> wt_max_wt(N/2, 0.0f);
    morph::vvec<float> wt_max_ki(N/2, 0.0f);
    morph::vvec<float> het_max_wt(N/2, 0.0f);
    morph::vvec<float> het_max_ki(N/2, 0.0f);
    morph::vvec<float> hom_max_wt(N/2, 0.0f);
    morph::vvec<float> hom_max_ki(N/2, 0.0f);
    // even
    int i = 0;
    morph::vvec<float> tmp(N, 0.0f);
    for (int x = 0; x < N; x+=2, ++i) {
        tmp.zero();
        for (int y = 0; y < N; ++y) { tmp[y] = pd_wt[y*N + x]; }
        wt_max_wt[i] = tmp.argmax();
        tmp.zero();
        for (int y = 0; y < N; ++y) { tmp[y] = pd_het[y*N + x]; }
        het_max_wt[i] = tmp.argmax();
        tmp.zero();
        for (int y = 0; y < N; ++y) { tmp[y] = pd_hom[y*N + x]; }
        hom_max_wt[i] = tmp.argmax();
    }
    // Odd (knockin)
    i = 0;
    tmp.zero();
    for (int x = 1; x < N; x+=2, ++i) { // retinal origin
        tmp.zero();
        for (int y = 0; y < N; ++y) { tmp[y] = pd_wt[y*N + x]; }
        wt_max_ki[i] = tmp.argmax();
        tmp.zero();
        for (int y = 0; y < N; ++y) { tmp[y] = pd_het[y*N + x]; }
        het_max_ki[i] = tmp.argmax();
        tmp.zero();
        for (int y = 0; y < N; ++y) { tmp[y] = pd_hom[y*N + x]; }
        hom_max_ki[i] = tmp.argmax();
    }

    morph::DatasetStyle ds(morph::stylepolicy::markers);
    ds.markercolour = morph::colour::crimson;
    ds.datalabel = "normal";
    morph::DatasetStyle ds2(morph::stylepolicy::markers);
    ds2.markercolour = morph::colour::black;
    ds2.datalabel = "knock-in";


    graphs_wt.g_prob_dist_maxima->setdata (wt_x, wt_max_wt, ds);
    graphs_wt.g_prob_dist_maxima->setdata (ki_x, wt_max_ki, ds2);
    graphs_wt.g_prob_dist_maxima->reinit();

    graphs_het.g_prob_dist_maxima->setdata (wt_x, het_max_wt, ds);
    graphs_het.g_prob_dist_maxima->setdata (ki_x, het_max_ki, ds2);
    graphs_het.g_prob_dist_maxima->reinit();

    graphs_hom.g_prob_dist_maxima->setdata (wt_x, hom_max_wt, ds);
    graphs_hom.g_prob_dist_maxima->setdata (ki_x, hom_max_ki, ds2);
    graphs_hom.g_prob_dist_maxima->reinit();

    std::cout << "Done accumulating probability distributions\n";

    v.keepOpen();

    return 0;
}
