/*
 * A net of locations, with information about what their ideal nearest neighbours are.
 */

#pragma once

#include <vector>
#include <morph/Vector.h>

template<typename T>
struct net // : public tissue?
{
    size_t w = 0;
    size_t h = 0;

    //! Initialize a rectangular net of width w and height h. Note that this resizes p
    //! and arranges c, but does not fill p with positions.
    void init (size_t _w, size_t _h)
    {
        this->w = _w;
        this->h = _h;
        this->p.resize(_w*_h);
        this->targ.resize(_w*_h);
        this->clr.resize(_w*_h);
        // Set up colours. Hack. Hardcoded.
        size_t i = 0;
        for (size_t y = 0; y < _h; ++y) {
            for (size_t x = 0; x < _w; ++x) {
                this->clr[i++] = {(float)x/(float)w, (float)y/(float)h, 0.0f};
            }
        }
        // Set up connections
        for (size_t y = 0; y < _h-1; ++y) {
            for (size_t x = 0; x < _w; ++x) {
                // Connect down
                this->c.insert (morph::Vector<size_t, 2>({x+y*_w, x+(y+1)*_w}));
            }
        }
        for (size_t x = 0; x < _w-1; ++x) {
            for (size_t y = 0; y < _h; ++y) {
                // Connect right
                this->c.insert (morph::Vector<size_t, 2>({x+y*_w, 1+x+y*_w}));
            }
        }
    }

    //! Return sum of squared distances between p and targ
    T sos() { return (p-targ).sos()[0]; }

    //! Positions of the vertices of the net
    morph::vVector<morph::Vector<T, 3>> p;
    //! Colours of the vertices of the net
    std::vector<std::array<float, 3>> clr;
    //! Connections of the net. The indices into p that are the ends of line segments
    std::set<morph::Vector<size_t, 2>> c;
    //! Target positions of the net, to allow the net to be used to provide a SOS metric for how close the pattern is to ideal
    morph::vVector<morph::Vector<T, 3>> targ;

    //! Apply a graft rotation to the targs
    void targ_graftrotate (morph::Vector<size_t,2> x1, size_t sz, size_t quarter_rotations)
    {
        if (quarter_rotations%4 == 0) { return; }

        // Potential for runtime crash here, so check x1 and sz.
        if (x1[0] + sz > this->w || x1[1] + sz > this->h) {
            throw std::runtime_error ("Patch for rotation is off the edge of the net");
        }

        morph::vVector<morph::Vector<T,3>> targ_graft (sz*sz);

        for (size_t rot = 1; rot <= quarter_rotations; ++rot) {
            // Copy the square
            for (size_t j = 0; j < sz; ++j) { // y
                for (size_t i = 0; i < sz; ++i) { // x
                    targ_graft[i+sz*j] = this->targ[(j+x1[1])*this->w+(i+x1[0])];
                }
            }
            // Write back, rotating 90 clockwise once
            for (size_t j = 0; j < sz; ++j) {
                for (size_t i = 0; i < sz; ++i) {
                    size_t ridx = (sz-1-j) + (i*sz);
                    this->targ[(j+x1[1])*this->w+(i+x1[0])] = targ_graft[ridx];
                }
            }
        }
    }

};
