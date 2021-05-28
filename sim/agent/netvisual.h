/*!
 * \file
 *
 * Visualise a network of locations on a rectangular grid, with lines to their intended
 * neighbours to the north, south, east and west. Intended to reproduce the
 * visualisations in the Simpson and Goodhill paper.
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

#include "net.h"

enum class netvisual_viewmode
{
    actual, // Show a visualisation of actual locations and relations
    target // Show a visualisation of target locations and relations
};

template <typename Flt>
class NetVisual : public morph::VisualModel
{
public:
    NetVisual(GLuint sp, GLuint tsp, const morph::Vector<float, 3> _offset, net<Flt>* _locations)
    {
        this->locations = _locations;
        this->shaderprog = sp;
        this->tshaderprog = tsp;
        this->mv_offset = _offset;
        this->viewmatrix.translate (this->mv_offset);
    }

    void initializeVertices()
    {
        if (this->viewmode == netvisual_viewmode::target) {
            this->initv_target();
        } else {
            this->initv_actual();
        }
    }

    //! Draw vertices for the net's actual locations p
    void initv_actual()
    {
        VBOint idx = 0;

        // Discs at the net vertices
        morph::Vector<float,3> puckthick = { 0, 0, 0.002 };
        for (unsigned int i = 0; i < this->locations->p.size(); ++i) {
            this->computeTube (idx,
                               this->locations->p[i]+puckthick,
                               this->locations->p[i]-puckthick,
                               morph::Vector<float,3>({1,0,0}), morph::Vector<float,3>({0,1,0}),
                               this->locations->clr[i], this->locations->clr[i],
                               this->radiusFixed, 16);
        }
        // Connection lines
        for (auto c : this->locations->c) {
            morph::Vector<Flt, 3> c1 = this->locations->p[c[0]];
            morph::Vector<Flt, 3> c2 = this->locations->p[c[1]];
            std::array<float, 3> clr1 = this->locations->clr[c[0]];
            std::array<float, 3> clr2 = this->locations->clr[c[1]];
            if ((c1-c2).length() < Flt{0.2}) {
                // omit long lines for a bit more clarity
                this->computeLine (idx, c1, c2, this->uz, clr1, clr2, this->linewidth, this->linewidth/4);
            }
        }
    }

    //! Draw where the net vertices are EXPECTED (according to their targ attribute)
    void initv_target()
    {
        VBOint idx = 0;


        for (unsigned int i = 0; i < this->locations->p.size(); ++i) {
            this->computeTube (idx,
                               this->locations->targ[i]+puckthick,
                               this->locations->targ[i]-puckthick,
                               morph::Vector<float,3>({1,0,0}), morph::Vector<float,3>({0,1,0}),
                               this->locations->clr[i], this->locations->clr[i],
                               this->radiusFixed, 16);
        }
        // Connections
        for (auto c : this->locations->c) {
            morph::Vector<Flt, 3> c1 = this->locations->targ[c[0]];
            morph::Vector<Flt, 3> c2 = this->locations->targ[c[1]];
            std::array<float, 3> clr1 = this->locations->clr[c[0]];
            std::array<float, 3> clr2 = this->locations->clr[c[1]];
            //if ((c1-c2).length() < Flt{0.2}) {
            this->computeLine (idx, c1, c2, this->uz, clr1, clr2, this->linewidth, this->linewidth/4);
            //}
        }
    }

    //! Set this->radiusFixed, then re-compute vertices.
    void setRadius (float fr)
    {
        this->radiusFixed = fr;
        this->reinit();
    }

    //! Pointer to a vector of locations to visualise
    net<Flt>* locations = (net<Flt>*)0;
    Flt radiusFixed = 0.01;
    Flt linewidth = 0.004;
    morph::Vector<float,3> puckthick = { 0, 0, 0.002 };
    // What to show
    netvisual_viewmode viewmode = netvisual_viewmode::actual;
    //! A normal vector, fixed as pointing up
    morph::Vector<float, 3> uz = {0,0,1};
};
