/*!
 * \file
 *
 * Visualise a bunch of agents, each of which has a history of locations
 * that it has visited previously, shown as lines. A coloured cap is used to indicate
 * EphA expression level. Saturn rings can be used to show interaction radius
 *
 * \author Seb James
 * \date 2021
 */
#pragma once

#ifdef __OSX__
# include <OpenGL/gl3.h>
#else
# include <GL3/gl3.h>
#endif
#include <morph/VisualModel.h>
#include <morph/Vector.h>
#include <morph/Scale.h>
#include <vector>
#include <array>
#include <set>

#include "branch.h"

enum class branchvisual_view
{
    detailed,  // Show receptor expression info as coloured pies and interaction radii
    discs,     // Just show a disc, of size the largest interaction radius
    axonview,  // discs with a grey history path
    discint    // disc and interaction hoop
};

template <typename Flt, size_t N, typename B=branch<Flt, N>>
class BranchVisual : public morph::VisualModel
{
public:
    BranchVisual(GLuint sp, GLuint tsp, const morph::Vector<float, 3> _offset,
                 std::vector<B>* _branches,
                 std::map<size_t, morph::vVector<morph::Vector<Flt, 3>>>* _ax_history,
                 branchvisual_view _view = branchvisual_view::detailed)
    {
        this->branches = _branches;
        this->ax_history = _ax_history;
        this->view = _view;
        this->shaderprog = sp;
        this->tshaderprog = tsp;
        this->mv_offset = _offset;
        this->viewmatrix.translate (this->mv_offset);
    }

    // Receptor and ligand need to be scaled
    morph::Scale<Flt, Flt> rcpt_scale;
    // Target location will have a different scale
    morph::Scale<Flt, Flt> target_scale;

    void drawBoundary (VBOint& idx)
    {
        std::array<float, 3> gry = { 0.2, 0.2, 0.2 };
        float w = 1, h = 1;
        this->computeLine (idx, morph::Vector<Flt, 3>({0, 0, 0}), morph::Vector<Flt, 3>({w, 0, 0}),
                           this->uz, gry, gry, this->blinewidth, this->blinewidth/4);
        this->computeLine (idx, morph::Vector<Flt, 3>({w, 0, 0}), morph::Vector<Flt, 3>({w, h, 0}),
                           this->uz, gry, gry, this->blinewidth, this->blinewidth/4);
        this->computeLine (idx, morph::Vector<Flt, 3>({w, h, 0}), morph::Vector<Flt, 3>({0, h, 0}),
                           this->uz, gry, gry, this->blinewidth, this->blinewidth/4);
        this->computeLine (idx, morph::Vector<Flt, 3>({0, h, 0}), morph::Vector<Flt, 3>({0, 0, 0}),
                           this->uz, gry, gry, this->blinewidth, this->blinewidth/4);
    }

    void initializeVertices()
    {
        if (this->view == branchvisual_view::axonview) { return this->initializeAxonview(); }

        // Set radius from the first branch
        if (!this->branches->empty()) {
            this->radiusFixed = (*this->branches)[0].getr();

            morph::Vector<Flt, 3> rads = { (*this->branches)[0].getr_c(), (*this->branches)[0].getr_j(), (*this->branches)[0].getr_i() };
            this->rad_interaction = rads.max() > this->radiusFixed ? rads.max() : Flt{0.0};
            //std::cout << "radiusFixed: " << radiusFixed << " rad_interaction: " << rad_interaction << std::endl;
        }

        VBOint idx = 0;
        // For each branch,simply draw a sphere for the current location, with a second
        // colour for the rcpt expression.
        for (auto b : *this->branches) {
            // Colour should come from original target location, rather than receptor value, to emphasise swaps in location.
            std::array<float, 3> clr = { this->target_scale.transform_one(b.target[0]),
                                         this->target_scale.transform_one(b.target[1]), 0 }; // position colour

            // A sphere at the last location. Tune number of rings (second last arg) in
            // sphere to change size of clr2 disc at top
            morph::Vector<float, 3> cur = { b.current[0], b.current[1], 0 };
            if (this->view == branchvisual_view::detailed) {
                std::array<float, 3> clr_r0 = { this->rcpt_scale.transform_one(b.rcpt[0]), 0, 0 }; // rcpt0 expression colour is red
                std::array<float, 3> clr_r1 = { 0, this->rcpt_scale.transform_one(b.rcpt[1]), 0 };
                std::array<float, 3> clr_r2 = { 0, 0, this->rcpt_scale.transform_one(b.rcpt[2]) };
                std::array<float, 3> clr_r3 = { this->rcpt_scale.transform_one(b.rcpt[3]),
                                                0,
                                                this->rcpt_scale.transform_one(b.rcpt[3]) };

                std::array<float, 3> clr_l0 = { this->rcpt_scale.transform_one(b.lgnd[0]), 0, 0 };
                std::array<float, 3> clr_l1 = { 0, this->rcpt_scale.transform_one(b.lgnd[1]), 0 };
                std::array<float, 3> clr_l2 = { 0, 0, this->rcpt_scale.transform_one(b.lgnd[2]) };
                std::array<float, 3> clr_l3 = { this->rcpt_scale.transform_one(b.lgnd[3]),
                                                0,
                                                this->rcpt_scale.transform_one(b.lgnd[3]) };
                this->computeDiscA (idx, cur, clr,
                                    clr_r0, clr_r1, clr_r2, clr_r3,
                                    clr_l0, clr_l1, clr_l2, clr_l3,
                                    this->radiusFixed);
                // Draw ring for the interaction, if it's larger than radiusFixed
                if (rad_interaction > 0.0f) {
                    this->computeRing (idx, cur, clr, this->rad_interaction, 0.1f*this->rad_interaction, 18);
                }

            } else if (this->view == branchvisual_view::discs) {
                this->computeTube (idx, cur, cur-morph::Vector<float,3>({0,0,this->radiusFixed*0.1f}), this->ux, this->uy,
                                   clr, clr,
                                   this->radiusFixed, 20);
            } else if (this->view == branchvisual_view::discint) {
                this->computeTube (idx, cur, cur-morph::Vector<float,3>({0,0,this->radiusFixed*0.1f}), this->ux, this->uy,
                                   clr, clr,
                                   this->radiusFixed, 20);
                if (this->rad_interaction > 0.0f) {
                    this->computeRing (idx, cur, clr, this->rad_interaction, 0.1f*this->rad_interaction, 20);
                }
            }
        }
    }

    // A special puck/disc with colours to show receptor and ligand expression
    void computeDiscA (VBOint& idx, morph::Vector<float> so,
                       std::array<float, 3> sc, // Main colour
                       std::array<float, 3> r0, // Rcpt 0 colour
                       std::array<float, 3> r1,
                       std::array<float, 3> r2,
                       std::array<float, 3> r3,
                       std::array<float, 3> l0, // ligand 0 colour
                       std::array<float, 3> l1,
                       std::array<float, 3> l2,
                       std::array<float, 3> l3,
                       float r = 1.0f, int segments = 16)
    {
        float r_1 = 2.6f * r/10.0f; // centre circle - rcpt
        float r_2 = 3.0f * r/10.0f; // white ring
        float r_3 = 5.0f * r/10.0f; // ligand ring
        float r_4 = 5.4f * r/10.0f; // outer white ring

        std::array<VBOint, 4> capMiddle;
        std::array<std::array<float, 3>, 4> rclrs = { r0, r1, r2, r3 };
        std::array<std::array<float, 3>, 4> lclrs = { l0, l1, l2, l3 };
        std::array<float, 3> white = {1,1,1};
        // Push the central point
        this->vertex_push (so[0], so[1], so[2], this->vertexPositions);
        this->vertex_push (this->uz, this->vertexNormals);
        this->vertex_push (r0, this->vertexColors);
        capMiddle[0] = idx++;
        this->vertex_push (so[0], so[1], so[2], this->vertexPositions);
        this->vertex_push (this->uz, this->vertexNormals);
        this->vertex_push (r1, this->vertexColors);
        capMiddle[1] = idx++;
        this->vertex_push (so[0], so[1], so[2], this->vertexPositions);
        this->vertex_push (this->uz, this->vertexNormals);
        this->vertex_push (r2, this->vertexColors);
        capMiddle[2] = idx++;
        this->vertex_push (so[0], so[1], so[2], this->vertexPositions);
        this->vertex_push (this->uz, this->vertexNormals);
        this->vertex_push (r3, this->vertexColors);
        capMiddle[3] = idx++;

        // Centre circle
        for (int jj = 0; jj < 4; ++jj) { // jj indexes quarters
            int quarlen = segments/4;
            int quarstart = jj*quarlen;
            int quarend = quarstart + quarlen;
            for (int j = quarstart; j < quarend; ++j) {

                float segment = 2 * M_PI * (float) (j) / segments;
                float x1 = r_1*cos(segment);
                float y1 = r_1*sin(segment);

                int jnext = (j+1)%segments;
                float segnext = 2 * M_PI * (float) (jnext) / segments;
                float x1n = r_1*cos(segnext);
                float y1n = r_1*sin(segnext);

                this->indices.push_back (capMiddle[jj]);

                this->vertex_push (so[0]+x1, so[1]+y1, so[2], this->vertexPositions);
                this->vertex_push (this->uz, this->vertexNormals);
                this->vertex_push (rclrs[jj], this->vertexColors);
                this->indices.push_back (idx++);

                this->vertex_push (so[0]+x1n, so[1]+y1n, so[2], this->vertexPositions);
                this->vertex_push (this->uz, this->vertexNormals);
                this->vertex_push (rclrs[jj], this->vertexColors);
                this->indices.push_back (idx++);
            }
        }

        // White ring
        for (int j = 0; j < segments; j++) {
            float segment = 2 * M_PI * (float) (j) / segments;
            float xin = r_1 * cos(segment);
            float yin = r_1 * sin(segment);
            float xout = r_2 * cos(segment);
            float yout = r_2 * sin(segment);
            int segjnext = (j+1) % segments;
            float segnext = 2 * M_PI * (float) (segjnext) / segments;
            float xin_n = r_1 * cos(segnext);
            float yin_n = r_1 * sin(segnext);
            float xout_n = r_2 * cos(segnext);
            float yout_n = r_2 * sin(segnext);
            // Quad corners
            morph::Vector<float> c1 = {xin, yin, 0};
            morph::Vector<float> c2 = {xout, yout, 0};
            morph::Vector<float> c3 = {xout_n, yout_n, 0};
            morph::Vector<float> c4 = {xin_n, yin_n, 0};
            this->computeFlatQuad (idx, so+c1, so+c2, so+c3, so+c4, white);
        }

        // Ligands ring
        for (int j = 0; j < segments; j++) {
            int jj = j < segments/2 ? (j < segments/4 ? 0 : 1) : (j < 3*segments/4 ? 2 : 3);
            float segment = 2 * M_PI * (float) (j) / segments;
            float xin = r_2 * cos(segment);
            float yin = r_2 * sin(segment);
            float xout = r_3 * cos(segment);
            float yout = r_3 * sin(segment);
            int segjnext = (j+1) % segments;
            float segnext = 2 * M_PI * (float) (segjnext) / segments;
            float xin_n = r_2 * cos(segnext);
            float yin_n = r_2 * sin(segnext);
            float xout_n = r_3 * cos(segnext);
            float yout_n = r_3 * sin(segnext);
            // Quad corners
            morph::Vector<float> c1 = {xin, yin, 0};
            morph::Vector<float> c2 = {xout, yout, 0};
            morph::Vector<float> c3 = {xout_n, yout_n, 0};
            morph::Vector<float> c4 = {xin_n, yin_n, 0};
            this->computeFlatQuad (idx, so+c1, so+c2, so+c3, so+c4, lclrs[jj]);
        }

        // White spacer ring
        for (int j = 0; j < segments; j++) {
            float segment = 2 * M_PI * (float) (j) / segments;
            float xin = r_3 * cos(segment);
            float yin = r_3 * sin(segment);
            float xout = r_4 * cos(segment);
            float yout = r_4 * sin(segment);
            int segjnext = (j+1) % segments;
            float segnext = 2 * M_PI * (float) (segjnext) / segments;
            float xin_n = r_3 * cos(segnext);
            float yin_n = r_3 * sin(segnext);
            float xout_n = r_4 * cos(segnext);
            float yout_n = r_4 * sin(segnext);
            // Quad corners
            morph::Vector<float> c1 = {xin, yin, 0};
            morph::Vector<float> c2 = {xout, yout, 0};
            morph::Vector<float> c3 = {xout_n, yout_n, 0};
            morph::Vector<float> c4 = {xin_n, yin_n, 0};
            this->computeFlatQuad (idx, so+c1, so+c2, so+c3, so+c4, white);
        }

        // Outer ring
        for (int j = 0; j < segments; j++) {
            float segment = 2 * M_PI * (float) (j) / segments;
            float xin = r_4 * cos(segment);
            float yin = r_4 * sin(segment);
            float xout = r * cos(segment);
            float yout = r * sin(segment);
            int segjnext = (j+1) % segments;
            float segnext = 2 * M_PI * (float) (segjnext) / segments;
            float xin_n = r_4 * cos(segnext);
            float yin_n = r_4 * sin(segnext);
            float xout_n = r * cos(segnext);
            float yout_n = r * sin(segnext);
            // Quad corners
            morph::Vector<float> c1 = {xin, yin, 0};
            morph::Vector<float> c2 = {xout, yout, 0};
            morph::Vector<float> c3 = {xout_n, yout_n, 0};
            morph::Vector<float> c4 = {xin_n, yin_n, 0};
            this->computeFlatQuad (idx, so+c1, so+c2, so+c3, so+c4, sc);
        }
    }

    void initializeAxonview()
    {
        // Set radius from the first branch
        if (!this->branches->empty()) { this->radiusFixed = (*this->branches)[0].getr(); }

        VBOint idx = 0;

        if (this->ax_history->empty()) { return; }

        // Now vis meanpaths
        for (auto mp : *this->ax_history) {
            std::array<float, 3> clr = {0.5f, 0.5f, 0.5f};
            for (unsigned int i = 1; i < mp.second.size(); ++i) {
                // draw line from mp[i-1] to mp[i]
                this->computeFlatLineRnd (idx, mp.second[i-1], mp.second[i], this->uz, clr, this->linewidth, 0.0f, true, false);
            }
        }

        // For each branch, draw spheres for the current location, with a second colour
        // for the rcpt expression. Also lines from end of common path to each sphere.
        for (auto b : *this->branches) {
            if (!this->seeaxons.count(b.aid)) { continue; }
            std::array<float, 3> clr = { this->target_scale.transform_one(b.target[0]),
                                         this->target_scale.transform_one(b.target[1]), 0 }; // position colour
            std::array<float, 3> clr2 = { this->rcpt_scale.transform_one(b.rcpt[0]), 0,
                                          this->rcpt_scale.transform_one(b.rcpt[1]) }; // rcpt expr colour
            morph::Vector<float, 3> cur = { 0, 0, 0 };
            cur[0] = b.next[0]; // or path.back()?
            cur[1] = b.next[1];
            this->computeSphere (idx, cur, clr, clr2, this->radiusFixed, 14, 12);
            this->computeFlatLineRnd (idx, (*this->ax_history)[b.aid].back(), cur, this->uz, clr, this->linewidth/2.0f, 0.0f, true, false);
        }

        this->drawBoundary (idx);
    }

    //! Set this->radiusFixed, then re-compute vertices.
    void setRadius (float fr)
    {
        this->radiusFixed = fr;
        this->reinit();
    }

    //! Pointer to a vector of branches to visualise
    std::vector<B>* branches = (std::vector<B>*)0;
    //! Pointer to axon history, if required.
    std::map<size_t, morph::vVector<morph::Vector<Flt, 3>>>* ax_history = (std::map<size_t, morph::vVector<morph::Vector<Flt, 3>>>*)0;
    //! Container for axon centroids. Compute here or only vis here?
    //! Change this to get larger or smaller spheres.
    Flt radiusFixed = 0.01;
    Flt rad_interaction = 0.0; // For a ring showing the largest interaction radius
    Flt linewidth = 0.008;
    Flt blinewidth = 0.004;
    //! A normal vector, fixed as pointing up
    morph::Vector<float, 3> ux = {1,0,0};
    morph::Vector<float, 3> uy = {0,1,0};
    morph::Vector<float, 3> uz = {0,0,1};
    //! In axon view, show just a small selection of axons
    //bool axonview = false;
    branchvisual_view view = branchvisual_view::detailed;
    //! How much history to show at max in the default view?
    unsigned int histlen = 20;
    //! Identities of axons to show in the axon view
    std::set<int> seeaxons;
};
