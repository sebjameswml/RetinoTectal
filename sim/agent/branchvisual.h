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

    void initializeVertices()
    {
        if (this->axonview) { return this->initializeAxonview(); }

        VBOint idx = 0;
        // For each branch,simply draw a sphere for the current location, with a second
        // colour for the rcpt expression.
        for (auto b : *this->branches) {
            // Colour should come from original target location, rather than receptor value, to emphasise swaps in location.
            std::array<float, 3> clr = { this->target_scale.transform_one(b.target[0]),
                                         this->target_scale.transform_one(b.target[1]), 0 }; // position colour
            std::array<float, 3> clr2 = { this->rcpt_scale.transform_one(b.rcpt[0]), 0,
                                          this->rcpt_scale.transform_one(b.rcpt[1]) }; // rcpt expr colour
            // A sphere at the last location. Tune number of rings (second last arg) in
            // sphere to change size of clr2 disc at top
            morph::Vector<float, 3> cur = { b.current[0], b.current[1], 0 };
            this->computeSphere (idx, cur, clr, clr2, this->radiusFixed, 14, 12);
        }
    }

    void initializeAxonview()
    {
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
            std::array<float, 3> clr = { this->rcpt_scale.transform_one(b.rcpt[0]),
                                         this->rcpt_scale.transform_one(b.rcpt[1]), 0 }; // position colour
            std::array<float, 3> clr2 = { this->rcpt_scale.transform_one(b.rcpt[0]), 0,
                                          this->rcpt_scale.transform_one(b.rcpt[1]) }; // rcpt expr colour
            morph::Vector<float, 3> cur = { 0, 0, 0 };
            cur[0] = b.next[0]; // or path.back()?
            cur[1] = b.next[1];
            this->computeSphere (idx, cur, clr, clr2, this->radiusFixed, 14, 12);
            this->computeFlatLineRnd (idx, (*this->ax_history)[b.aid].back(), cur, this->uz, clr, this->linewidth/2.0f, 0.0f, true, false);
        }
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
