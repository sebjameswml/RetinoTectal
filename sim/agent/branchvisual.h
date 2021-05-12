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

template <typename Flt, size_t N>
class BranchVisual : public morph::VisualModel
{
public:
    BranchVisual(GLuint sp, GLuint tsp, const morph::Vector<float, 3> _offset,
                 std::vector<branch<Flt, N>>* _branches, bool _axonview = false)
    {
        this->branches = _branches;
        this->axonview = _axonview;
        this->shaderprog = sp;
        this->tshaderprog = tsp;
        this->mv_offset = _offset;
        this->viewmatrix.translate (this->mv_offset);
    }

    // Receptor and ligand need to be scaled to 0/1 range
    morph::Scale<Flt, Flt> rcpt_scale;

    void initializeVertices()
    {
        if (this->axonview) { return this->initializeAxonview(); }

        VBOint idx = 0;
        // For each branch, draw lines for the path history and a sphere for the current
        // location, with a second colour for the rcpt expression.
        for (auto b : *this->branches) {
            // Colour comes from target location.
            // Prolly need receptor scaling here:
            std::array<float, 3> clr = { this->rcpt_scale.transform_one(b.rcpt[0]),
                                         this->rcpt_scale.transform_one(b.rcpt[1]), 0 };
            std::array<float, 3> clr2 = { 0, 0, this->rcpt_scale.transform_one(b.rcpt[0]) };
            morph::Vector<float, 3> cur = { 0, 0, 0 };
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

        std::map<int, morph::vVector<morph::Vector<Flt, 2>>> meanpaths;
        std::map<int, std::array<float, 3>> mpcolours;

        // Add up the path to get the mean path for the
        // axon, to replicate the visualization in Simpson & Goodhill.
        for (auto b : *this->branches) {

            if (!this->seeaxons.count(b.aid)) { continue; }
            std::array<float, 3> clr = { this->rcpt_scale.transform_one(b.rcpt[0]),
                                         this->rcpt_scale.transform_one(b.rcpt[1]), 0 };
            mpcolours[b.aid] = clr;
            if (meanpaths.count(b.aid)== 0) {
                // b.path/b.bpa divides each Vector element of the vVector path by bpa (a scalar)
                meanpaths[b.aid] = (b.path/static_cast<Flt>(this->bpa));
            } else {
                meanpaths[b.aid] += (b.path/static_cast<Flt>(this->bpa));
            }
        }

        // Now vis meanpaths
        for (auto mp : meanpaths) {
            std::array<float, 3> clr = mpcolours[mp.first];
            morph::Vector<float, 3> last = { 0, 0, 0 };
            morph::Vector<float, 3> cur = { 0, 0, 0 };
            for (unsigned int i = 1; i < mp.second.size(); ++i) {
                // draw line from mp[i-1] to mp[i]
                last[0] = mp.second[i-1][0];
                last[1] = mp.second[i-1][1];
                cur[0] = mp.second[i][0];
                cur[1] = mp.second[i][1];
                this->computeFlatLineRnd (idx, last, cur, this->uz, clr, this->linewidth, 0.0f, true, false);
            }
        }

        // For each branch, draw spheres for the current location, with a second colour
        // for the rcpt expression. Also lines from end of common path to each sphere.
        for (auto b : *this->branches) {
            if (!this->seeaxons.count(b.aid)) { continue; }
            std::array<float, 3> clr = { this->rcpt_scale.transform_one(b.rcpt[0]),
                                         this->rcpt_scale.transform_one(b.rcpt[1]), 0 };
            std::array<float, 3> clr2 = { 0, 0, this->rcpt_scale.transform_one(b.rcpt[0]) };
            morph::Vector<float, 3> cur = { 0, 0, 0 };
            cur[0] = b.path.back()[0];
            cur[1] = b.path.back()[1];
            this->computeSphere (idx, cur, clr, clr2, this->radiusFixed, 14, 12);
            morph::Vector<float, 3> last = { 0, 0, 0 };
            last[0] = meanpaths[b.aid].back()[0];
            last[1] = meanpaths[b.aid].back()[1];
            this->computeFlatLineRnd (idx, last, cur, this->uz, clr, this->linewidth/2.0f, 0.0f, true, false);
        }
    }

    //! Set this->radiusFixed, then re-compute vertices.
    void setRadius (float fr)
    {
        this->radiusFixed = fr;
        this->reinit();
    }

    //! Pointer to a vector of branches to visualise
    std::vector<branch<Flt, N>>* branches = (std::vector<branch<Flt, N>>*)0;
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
    //! Branches per axon. Required for the 'single axon view' mode
    unsigned int bpa = 8;
};
