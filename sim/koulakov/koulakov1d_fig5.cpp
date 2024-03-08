// Model
#include "koulakov1d.h"

// Some visualization helpers
#include "koulakov1d_figs_5_6.h"

#include <iostream>
#include <morph/vec.h>
#include <morph/vvec.h>
#include <morph/Visual.h>
#include <morph/GraphVisual.h>
#include <morph/CartGrid.h>

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
                                                 morph::GridDomainShape::Rectangle,
                                                 morph::GridDomainWrap::Horizontal);
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
