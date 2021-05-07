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

template <typename Flt>
class BranchVisual : public morph::VisualModel
{
public:
    BranchVisual(GLuint sp, const morph::Vector<float, 3> _offset,
                 std::vector<branch<Flt>>* _branches, bool _axonview = false)
    {
        this->branches = _branches;
        this->axonview = _axonview;
        this->shaderprog = sp;
        this->mv_offset = _offset;
        this->viewmatrix.translate (this->mv_offset);
    }

    morph::Scale<Flt, Flt> EphA_scale;

    void initializeVertices()
    {
        if (this->axonview) { return this->initializeAxonview(); }

        VBOint idx = 0;
        // For each branch, draw lines for the path history and a sphere for the current
        // location, with a second colour for the EphA expression.
        for (auto b : *this->branches) {
            // Colour comes from target location.
            std::array<float, 3> clr = { b.tz[0], b.tz[1], 0 };
            std::array<float, 3> clr2 = { 0, 0, this->EphA_scale.transform_one(b.EphA) };
            morph::Vector<float, 3> last = { 0, 0, 0 };
            morph::Vector<float, 3> cur = { 0, 0, 0 };
#if 0
            // First draw the path. Problem here is where to *start* in the history
            for (unsigned int i = 1; i < b.path.size() && i < this->histlen; ++i) {
                last[0] = b.path[i-1][0];
                last[1] = b.path[i-1][1];
                cur[0] = b.path[i][0];
                cur[1] = b.path[i][1];
                this->computeFlatLineRnd (idx, last, cur, this->uz, clr, this->linewidth, 0.0f, true, false);
            }
#endif
            // Finally, a sphere at the last location. Tune number of rings (second last
            // arg) in sphere to change size of clr2 disc at top
            morph::Vector<float, 2> bk = b.path.back();
            cur[0] = bk[0];
            cur[1] = bk[1];
            this->computeSphere (idx, cur, clr, clr2, this->radiusFixed, 14, 12);
        }
    }

    void initializeAxonview()
    {
        VBOint idx = 0;
        // For each branch, draw lines for the path history and a sphere for the current
        // location, with a second colour for the EphA expression.
        for (auto b : *this->branches) {

            if (!this->seeaxons.count(b.aid)) { continue; }

            // Colour comes from target location.
            std::array<float, 3> clr = { b.tz[0], b.tz[1], 0 };
            std::array<float, 3> clr2 = { 0, 0, this->EphA_scale.transform_one(b.EphA) };
            morph::Vector<float, 3> last = { 0, 0, 0 };
            morph::Vector<float, 3> cur = { 0, 0, 0 };
            // First draw the path
            for (unsigned int i = 1; i < b.path.size(); ++i) {
                last[0] = b.path[i-1][0];
                last[1] = b.path[i-1][1];
                cur[0] = b.path[i][0];
                cur[1] = b.path[i][1];
                this->computeFlatLineRnd (idx, last, cur, this->uz, clr, this->linewidth, 0.0f, true, false);
            }
            // Finally, a sphere at the last location. Tune number of rings (second last
            // arg) in sphere to change size of clr2 disc at top
            this->computeSphere (idx, cur, clr, clr2, this->radiusFixed, 14, 12);
        }
    }

    //! Set this->radiusFixed, then re-compute vertices.
    void setRadius (float fr)
    {
        this->radiusFixed = fr;
        this->reinit();
    }

    //! Pointer to a vector of branches to visualise
    std::vector<branch<Flt>>* branches = (std::vector<branch<Flt>>*)0;
    //! Container for axon centroids. Compute here or only vis here?
    //! Change this to get larger or smaller spheres.
    Flt radiusFixed = 0.01;
    Flt linewidth = 0.008;
    //! A normal vector, fixed as pointing up
    morph::Vector<float, 3> uz = {0,0,1};
    //! In axon view, show just a small selection of axons
    bool axonview = false;
    //! How much history to show at max in the default view?
    unsigned int histlen = 20;
    //! Identities of axons to show in the axon view
    std::set<int> seeaxons;
};
