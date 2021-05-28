/*
 * A net of locations, with information about what their ideal nearest neighbours are.
 */

#pragma once

#include <vector>
#include <morph/Vector.h>

template<typename T>
struct net
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
};

// A retinal ganglion cell net specialised to know about graft rotations/swaps etc
template<typename T>
struct rgcnet : public net<T>
{
    // Set true if the x/y coordinates in targ have been mirrrored
    bool mirrored = false;

    //! Expand the top half down across the full tissue area
    void targ_expand_topdown() { for (auto& t : this->targ) { t[1] = t[1] * T{2} - T{1}; } }

    //! Tectum has reduced in size so squish the targs up
    void targ_squish_bottomup() { for (auto& t : this->targ) { t[1] = t[1] * T{0.5} + T{0.5}; } }

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

    // Tell this function if the target positions have been transformed with respect to
    // their location in the array, because it needs to know.
    void targ_graftswap (morph::Vector<size_t,2> x1,
                         morph::Vector<size_t,2> sz,
                         morph::Vector<size_t,2> x2)
    {
        // I copied this from tissue, which does the swap by index and so I need to
        // cover two situations; one where I swap row by row, and one where, if
        // mirrored, I swap col by col.
        if (mirrored == true) {
            // move transformed section of targ row by row.
            morph::vVector<morph::Vector<T,3>> graft (sz[1]);

            size_t idx1 = x1[1] + this->w * x1[0];
            size_t idx2 = x2[1] + this->w * x2[0];
            for (size_t i = 0; i < sz[0]; ++i) {

                if (idx1 + sz[1] > this->targ.size() || idx2 + sz[1] > this->targ.size()) {
                    std::cout << "Can't make that patchswap\n";
                    throw std::runtime_error ("Can't make that patchswap");
                }

                auto oi = this->targ.begin();
                oi += idx1;
                auto ni = this->targ.begin();
                ni += idx2;
                std::copy_n (ni, sz[1], graft.begin());
                std::copy_n (oi, sz[1], ni);
                std::copy_n (graft.begin(), sz[1], oi);

                idx1 += this->w;
                idx2 += this->w;
            }

        } else {
            // move section of targ row by row
            morph::vVector<morph::Vector<T,3>> graft (sz[0]);

            size_t idx1 = x1[0] + this->w * x1[1];
            size_t idx2 = x2[0] + this->w * x2[1];
            for (size_t i = 0; i < sz[1]; ++i) {

                if (idx1 + sz[0] > this->targ.size() || idx2 + sz[0] > this->targ.size()) {
                    throw std::runtime_error ("Can't make that patchswap");
                }

                auto oi = this->targ.begin();
                oi += idx1;
                auto ni = this->targ.begin();
                ni += idx2;
                std::copy_n (ni, sz[0], graft.begin());
                std::copy_n (oi, sz[0], ni);
                std::copy_n (graft.begin(), sz[0], oi);

                idx1 += this->w;
                idx2 += this->w;
            }
        }
    }
};
