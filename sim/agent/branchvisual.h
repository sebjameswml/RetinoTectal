/*!
 * \file
 *
 * Visualise a bunch of agents (as spheres), each of which has a history of locations
 * that it has visited previously, shown as lines. A coloured cap is used to indicate
 * EphA expression level.
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

template <typename Flt, size_t N, typename B=branch<Flt, N>>
class BranchVisual : public morph::VisualModel
{
public:
    BranchVisual(GLuint sp, GLuint tsp, const morph::Vector<float, 3> _offset,
                 std::vector<B>* _branches,
                 std::map<size_t, morph::vVector<morph::Vector<Flt, 3>>>* _ax_history,
                 bool _axonview = false)
    {
        this->branches = _branches;
        this->ax_history = _ax_history;
        this->axonview = _axonview;
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
        if (this->axonview) { return this->initializeAxonview(); }

        // Set radius from the first branch
        if (this->setRadiusFromBranches && !this->branches->empty()) {
            this->radiusFixed = (*this->branches)[0].getr();
            this->rad_competition = (*this->branches)[0].getrc() > this->radiusFixed ? (*this->branches)[0].getrc() : 0.0f;
            this->rad_interaction = (*this->branches)[0].getrrl() > this->radiusFixed ? (*this->branches)[0].getrrl() : 0.0f;
        }

        VBOint idx = 0;
        // For each branch,simply draw a sphere for the current location, with a second
        // colour for the rcpt expression.
        for (auto b : *this->branches) {
            // Colour should come from original target location, rather than receptor value, to emphasise swaps in location.
            std::array<float, 3> clr = { this->target_scale.transform_one(b.target[0]),
                                         this->target_scale.transform_one(b.target[1]), 0 }; // position colour

            std::array<float, 3> clr_r0 = { this->rcpt_scale.transform_one(b.rcpt[0]), 0, 0 }; // rcpt0 expression colour is red
            std::array<float, 3> clr_r1 = { 0, this->rcpt_scale.transform_one(b.rcpt[1]), 0 };
            std::array<float, 3> clr_r2 = { 0, 0, this->rcpt_scale.transform_one(b.rcpt[2]) };
            std::array<float, 3> clr_r3 = { this->rcpt_scale.transform_one(b.rcpt[3]),
                                            this->rcpt_scale.transform_one(b.rcpt[3]),
                                            this->rcpt_scale.transform_one(b.rcpt[3]) };

            std::array<float, 3> clr_l0 = { this->rcpt_scale.transform_one(b.lgnd[0]), 0, 0 };
            std::array<float, 3> clr_l1 = { 0, this->rcpt_scale.transform_one(b.lgnd[1]), 0 };
            std::array<float, 3> clr_l2 = { 0, 0, this->rcpt_scale.transform_one(b.lgnd[2]) };
            std::array<float, 3> clr_l3 = { this->rcpt_scale.transform_one(b.lgnd[3]),
                                            this->rcpt_scale.transform_one(b.lgnd[3]),
                                            this->rcpt_scale.transform_one(b.lgnd[3]) };

            // A sphere at the last location. Tune number of rings (second last arg) in
            // sphere to change size of clr2 disc at top
            morph::Vector<float, 3> cur = { b.current[0], b.current[1], 0 };
            //this->computeSphere (idx, cur, clr, clr2, this->radiusFixed, 14, 12);
            this->computeSphereA (idx, cur, clr,
                                  clr_r0, clr_r1, clr_r2, clr_r3,
                                  clr_l0, clr_l1, clr_l2, clr_l3,
                                  this->radiusFixed, 14, 12);

            // Draw ring for the competition interaction
            if (rad_competition > 0.0f) {
                this->computeRing (idx, cur, clr, this->rad_competition, 0.1f*this->rad_competition, 12);
            }
            if (rad_interaction > 0.0f) {
                this->computeRing (idx, cur, clr, this->rad_interaction, 0.1f*this->rad_interaction, 12);
            }
        }
    }

    // A specialised computeSphere with a funky colour cap to visualise the 4 receptors and 4 ligands
    void computeSphereA (VBOint& idx, morph::Vector<float> so,
                         std::array<float, 3> sc, // Main sphere colour
                         std::array<float, 3> r0, // Rcpt 0 colour
                         std::array<float, 3> r1,
                         std::array<float, 3> r2,
                         std::array<float, 3> r3,
                         std::array<float, 3> l0, // ligands 0
                         std::array<float, 3> l1,
                         std::array<float, 3> l2,
                         std::array<float, 3> l3,
                         float r = 1.0f, int rings = 10, int segments = 12)
    {
        // First cap, draw as a triangle fan, but record indices so that
        // we only need a single call to glDrawElements.
        float rings0 = M_PI * -0.5;
        float _z0  = sin(rings0);
        float z0  = r * _z0;
        float r_0 =  cos(rings0);
        float rings1 = M_PI * (-0.5 + 1.0f / rings);
        float _z1 = sin(rings1);
        float z1 = r * _z1;
        float r_1 = cos(rings1);

        // Push the central point
        this->vertex_push (so[0]+0.0f, so[1]+0.0f, so[2]+z0, this->vertexPositions);
        this->vertex_push (0.0f, 0.0f, -1.0f, this->vertexNormals);
        this->vertex_push (sc, this->vertexColors);
        VBOint capMiddle = idx++;
        VBOint ringStartIdx = idx;
        VBOint lastRingStartIdx = idx;

        // Top cap actually appears as bottom cap in the vis.
        bool firstseg = true;
        for (int j = 0; j < segments; j++) {

            float segment = 2 * M_PI * (float) (j) / segments;
            float x = cos(segment);
            float y = sin(segment);

            float _x1 = x*r_1;
            float x1 = _x1*r;
            float _y1 = y*r_1;
            float y1 = _y1*r;

            this->vertex_push (so[0]+x1, so[1]+y1, so[2]+z1, this->vertexPositions);
            this->vertex_push (_x1, _y1, _z1, this->vertexNormals);
            this->vertex_push (sc, this->vertexColors);

            if (!firstseg) {
                this->indices.push_back (capMiddle);
                this->indices.push_back (idx-1);
                this->indices.push_back (idx++);
            } else {
                idx++;
                firstseg = false;
            }
        }
        this->indices.push_back (capMiddle);
        this->indices.push_back (idx-1);
        this->indices.push_back (capMiddle+1);

        std::array<float, 3> fillclr = { 1.0f, 1.0f, 1.0f };

        // Now add the triangles around the rings, with appropriate colours
        for (int i = 2; i < rings; i++) {

            rings0 = M_PI * (-0.5 + (float) (i) / rings);
            _z0  = sin(rings0);
            z0  = r * _z0;
            r_0 =  cos(rings0);

            for (int j = 0; j < segments; j++) {

                // "current" segment
                float segment = 2 * M_PI * (float)j / segments;
                float x = cos(segment);
                float y = sin(segment);

                // One vertex per segment
                float _x0 = x*r_0;
                float x0 = _x0*r;
                float _y0 = y*r_0;
                float y0 = _y0*r;

                std::array<float, 3> jrcol = (j < segments/2 ?
                                              (j < segments/4 ? r0 : r1)
                                              :
                                              (j < segments*3/4 ? r2 : r3));

                std::array<float, 3> jlcol = (j < segments/2 ?
                                              (j < segments/4 ? l0 : l1)
                                              :
                                              (j < segments*3/4 ? l2 : l3));

                // NB: Only add ONE vertex per segment. ALREADY have the first ring!
                this->vertex_push (so[0]+x0, so[1]+y0, so[2]+z0, this->vertexPositions);
                // The vertex normal of a vertex that makes up a sphere is
                // just a normal vector in the direction of the vertex.
                this->vertex_push (_x0, _y0, _z0, this->vertexNormals);
                if (i== (rings-5)) {
                    this->vertex_push (fillclr, this->vertexColors);
                } else if (i==(rings-4) || i == (rings-3)) {
                    this->vertex_push (jlcol, this->vertexColors);
                } else if (i==(rings-2)) {
                    this->vertex_push (fillclr, this->vertexColors);
                } else if (i==(rings-1)) {
                    this->vertex_push (jrcol, this->vertexColors);
                } else {
                    this->vertex_push (sc, this->vertexColors);
                }
                if (j == segments - 1) {
                    // Last vertex is back to the start
                    this->indices.push_back (ringStartIdx++);
                    this->indices.push_back (idx);
                    this->indices.push_back (lastRingStartIdx);
                    this->indices.push_back (lastRingStartIdx);
                    this->indices.push_back (idx++);
                    this->indices.push_back (lastRingStartIdx+segments);
                } else {
                    this->indices.push_back (ringStartIdx++);
                    this->indices.push_back (idx);
                    this->indices.push_back (ringStartIdx);
                    this->indices.push_back (ringStartIdx);
                    this->indices.push_back (idx++);
                    this->indices.push_back (idx);
                }
            }
            lastRingStartIdx += segments;
        }

        // bottom cap
        rings0 = M_PI * 0.5;
        _z0  = sin(rings0);
        z0  = r * _z0;
        r_0 =  cos(rings0);
        // Push the central point of the bottom cap

        this->vertex_push (so[0]+0.0f, so[1]+0.0f, so[2]+z0, this->vertexPositions);
        this->vertex_push (0.0f, 0.0f, 1.0f, this->vertexNormals);
        this->vertex_push (r0, this->vertexColors);
        capMiddle = idx++;
        this->vertex_push (so[0]+0.0f, so[1]+0.0f, so[2]+z0, this->vertexPositions);
        this->vertex_push (0.0f, 0.0f, 1.0f, this->vertexNormals);
        this->vertex_push (r1, this->vertexColors);
        VBOint capMiddle2 = idx++;
        this->vertex_push (so[0]+0.0f, so[1]+0.0f, so[2]+z0, this->vertexPositions);
        this->vertex_push (0.0f, 0.0f, 1.0f, this->vertexNormals);
        this->vertex_push (r2, this->vertexColors);
        VBOint capMiddle3 = idx++;
        this->vertex_push (so[0]+0.0f, so[1]+0.0f, so[2]+z0, this->vertexPositions);
        this->vertex_push (0.0f, 0.0f, 1.0f, this->vertexNormals);
        this->vertex_push (r3, this->vertexColors);
        VBOint capMiddle4 = idx++;
        firstseg = true;
        // No more vertices to push, just do the indices for the bottom cap
        ringStartIdx = lastRingStartIdx;
        for (int j = 0; j < segments; j++) {
            if (j != segments - 1) {
                this->indices.push_back (
                              j < segments/2 ?
                              (j < segments/4 ? capMiddle : capMiddle2)
                              :
                               (j < 3*segments/4 ? capMiddle3 : capMiddle4)
                                  );
                this->indices.push_back (ringStartIdx++);
                this->indices.push_back (ringStartIdx);
            } else {
                // Last segment
                this->indices.push_back (capMiddle4);
                this->indices.push_back (ringStartIdx);
                this->indices.push_back (lastRingStartIdx);
            }
        }
    }

    void initializeAxonview()
    {
        // Set radius from the first branch
        if (this->setRadiusFromBranches && !this->branches->empty()) {
            this->radiusFixed = (*this->branches)[0].getr();
        }

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
    Flt rad_competition = 0.0; // For a ring showing competitive interaction radius
    Flt rad_interaction = 0.0; // For a ring showing recpt-lgnd interaction radius
    // If true, then use the r which is encoded in the first branch in branches to set radiusFixed. Else, don't change it.
    bool setRadiusFromBranches = true;
    Flt linewidth = 0.008;
    Flt blinewidth = 0.004;
    //! A normal vector, fixed as pointing up
    morph::Vector<float, 3> uz = {0,0,1};
    //! In axon view, show just a small selection of axons
    bool axonview = false;
    //! How much history to show at max in the default view?
    unsigned int histlen = 20;
    //! Identities of axons to show in the axon view
    std::set<int> seeaxons;
};
