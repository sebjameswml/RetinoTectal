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
#include <morph/vec.h>
#include <morph/Scale.h>
#include <vector>
#include <array>

#include "net.h"

enum class netvisual_viewmode
{
    actual, // Show a visualisation of actual locations and relations
    actual_nolines, // actual locations, but not the relationship lines
    target, // Show a visualisation of target locations and relations
    targetplus // Show target locations (with relationship net), and actual locations and lines connecting
};

template <typename Flt>
class NetVisual : public morph::VisualModel
{
public:
    NetVisual(const morph::vec<float, 3> _offset, net<Flt>* _locations)
    {
        this->locations = _locations;
        this->mv_offset = _offset;
        this->viewmatrix.translate (this->mv_offset);
    }

    void initializeVertices()
    {
        if (this->viewmode == netvisual_viewmode::target) {
            this->initv_target();
        } else if (this->viewmode == netvisual_viewmode::targetplus) {
            this->initv_target2();
        } else if (this->viewmode == netvisual_viewmode::actual
                   || this->viewmode == netvisual_viewmode::actual_nolines) {
            this->initv_actual();
        } else {
            throw std::runtime_error ("Unknown netvisual_viewmode");
        }
    }

    //! Draw vertices for the net's actual locations p
    void initv_actual()
    {
        // Discs at the net vertices
        morph::vec<float,3> puckthick = { 0, 0, 0.002 };
        for (unsigned int i = 0; i < this->locations->p.size(); ++i) {
            this->computeTube (this->idx,
                               (this->locations->p[i]+puckthick)*zoom,
                               (this->locations->p[i]-puckthick)*zoom,
                               morph::vec<float,3>({1,0,0}), morph::vec<float,3>({0,1,0}),
                               this->locations->clr[i], this->locations->clr[i],
                               this->radiusFixed*zoom, 16);
        }
        // Connection lines
        if (this->viewmode == netvisual_viewmode::actual) {
            for (auto c : this->locations->c) {
                morph::vec<Flt, 3> c1 = this->locations->p[c[0]] * zoom;
                morph::vec<Flt, 3> c2 = this->locations->p[c[1]] * zoom;
                std::array<float, 3> clr1 = this->locations->clr[c[0]];
                std::array<float, 3> clr2 = this->locations->clr[c[1]];
                if ((c1-c2).length() < maxlen) {
                    // omit long lines for a bit more clarity
                    this->computeLine (this->idx, c1, c2, this->uz, clr1, clr2,
                                       this->linewidth*zoom, this->linewidth/4*zoom);
                }
            }
        }
        this->drawBoundary();
    }

    //! Draw where the net vertices are EXPECTED (according to their targ attribute)
    void initv_target()
    {
        for (unsigned int i = 0; i < this->locations->p.size(); ++i) {
            this->computeTube (this->idx,
                               (this->locations->targ[i]+puckthick)*zoom,
                               (this->locations->targ[i]-puckthick)*zoom,
                               morph::vec<float,3>({1,0,0}), morph::vec<float,3>({0,1,0}),
                               this->locations->clr[i], this->locations->clr[i],
                               this->radiusFixed*zoom, 16);
        }
        // Connections
        for (auto c : this->locations->c) {
            morph::vec<Flt, 3> c1 = this->locations->targ[c[0]] * zoom;
            morph::vec<Flt, 3> c2 = this->locations->targ[c[1]] * zoom;
            std::array<float, 3> clr1 = this->locations->clr[c[0]];
            std::array<float, 3> clr2 = this->locations->clr[c[1]];
            this->computeLine (this->idx, c1, c2, this->uz, clr1, clr2, this->linewidth*zoom, this->linewidth/4*zoom);
        }

        this->drawBoundary();
    }

    //! Draw where the net vertices are EXPECTED (according to their targ attribute) AND
    //! where the centroid is using an arrow for the difference
    void initv_target2()
    {
        for (unsigned int i = 0; i < this->locations->p.size(); ++i) {

            // The puck for the target position
            this->computeTube (this->idx,
                               (this->locations->targ[i]+puckthick)*zoom,
                               (this->locations->targ[i]-puckthick)*zoom,
                               morph::vec<float,3>({1,0,0}), morph::vec<float,3>({0,1,0}),
                               this->locations->clr[i], this->locations->clr[i],
                               this->radiusFixed*zoom, 16);

            if (this->draw_actual == true) {
                // A line (cyl. tube) from target to actual position
                this->computeTube (this->idx,
                                   this->locations->targ[i]*zoom,
                                   (this->locations->p[i]+actualpuckoffs)*zoom,
                                   this->locations->clr[i], this->locations->clr[i],
                                   this->puckthick[2]*zoom, 8);

                // A slightly smaller puck for actual position
                this->computeTube (this->idx,
                                   (this->locations->p[i]+actualpuckoffs+puckthick)*zoom,
                                   (this->locations->p[i]+actualpuckoffs-puckthick)*zoom,
                                   morph::vec<float,3>({1,0,0}), morph::vec<float,3>({0,1,0}),
                                   this->locations->clr[i], this->locations->clr[i],
                                   this->radiusFixed*0.667*zoom, 16);
            }
        }
        // Connections
        for (auto c : this->locations->c) {
            morph::vec<Flt, 3> c1 = this->locations->targ[c[0]] * zoom;
            morph::vec<Flt, 3> c2 = this->locations->targ[c[1]] * zoom;
            std::array<float, 3> clr1 = this->locations->clr[c[0]];
            std::array<float, 3> clr2 = this->locations->clr[c[1]];
            this->computeLine (this->idx, c1, c2, this->uz, clr1, clr2, this->linewidth*zoom, this->linewidth/4*zoom);
        }

        // Finally, draw a square around the domain
        this->drawBoundary();
    }

    // The tissue boundary, indicated by dashed lines
    void drawBoundary()
    {
        std::array<float, 3> gry = { 0.2, 0.2, 0.2 };
        float w = zoom, h = zoom;
        if (this->locations->p.size() > 1) {
            w *= (this->locations->domain_w-1) * this->locations->dx[0];
            h *= (this->locations->domain_h-1) * this->locations->dx[1];
        }
        float _z = puckthick[2]*float{0.5001}; // Ensure boundary is visible above rest of drawing

        // Spins in here... (FIXME)
        this->computeFlatDashedLine (this->idx,
                                     morph::vec<Flt, 3>({0, 0, _z}),
                                     morph::vec<Flt, 3>({w, 0, _z}),
                                     this->uz,
                                     gry,
                                     this->linewidth*zoom, 0.0f,
                                     this->linewidth*5.0f, 0.4f);

        this->computeFlatDashedLine (this->idx, morph::vec<Flt, 3>({w, 0, _z}), morph::vec<Flt, 3>({w, h, _z}),
                                     this->uz,
                                     gry,
                                     this->linewidth*zoom, 0.0f,
                                     this->linewidth*5.0f, 0.4f);

        this->computeFlatDashedLine (this->idx, morph::vec<Flt, 3>({w, h, _z}), morph::vec<Flt, 3>({0, h, _z}),
                                     this->uz,
                                     gry,
                                     this->linewidth*zoom, 0.0f,
                                     this->linewidth*5.0f, 0.4f);

        this->computeFlatDashedLine (this->idx, morph::vec<Flt, 3>({0, h, _z}), morph::vec<Flt, 3>({0, 0, _z}),
                                     this->uz,
                                     gry,
                                     this->linewidth*zoom, 0.0f,
                                     this->linewidth*5.0f, 0.4f);
    }

    //! Set this->radiusFixed, then re-compute vertices.
    void setRadius (float fr)
    {
        this->radiusFixed = fr;
        this->reinit();
    }

    //! Pointer to a vector of locations to visualise
    net<Flt>* locations = nullptr;
    Flt radiusFixed = 0.01;
    Flt linewidth = 0.004;
    //! Zoom the size of the netvisual
    float zoom = float{1};
    //! If true, draw small pucks and lines showing actual position imposed on the 'expt suggests' view
    bool draw_actual = false;
    //! The maximum length of a line betwen two vertices for it to be visualised
    Flt maxlen = 1e9;
    morph::vec<float,3> puckthick = { 0, 0, 0.002 };
    // An offset to the puck for the actual position
    morph::vec<float,3> actualpuckoffs = { 0, 0, 0.01 };
    // What to show
    netvisual_viewmode viewmode = netvisual_viewmode::actual;
    //! A normal vector, fixed as pointing up
    morph::vec<float, 3> uz = {0,0,1};
};
