/*
 * A net of locations, with information about what their ideal nearest neighbours are.
 */

#pragma once

#include <vector>
#include <set>
#include <morph/vec.h>
#include <morph/vvec.h>
#include <morph/MathAlgo.h>
#include <morph/ColourMap.h>
#include "tissue.h"

// Comparison for std::set of vecs
struct ncmp
{
    bool operator() (morph::vec<size_t, 2> a, morph::vec<size_t, 2> b) const
    {
        return a.lexical_lessthan(b);
    }
};

template<typename T>
struct net
{
    size_t w = 0;
    size_t h = 0;

    // To aid with visualisation, store the distance for each unit step in the grid
    size_t domain_w = 0;
    size_t domain_h = 0;
    morph::vec<float, 2> dx = {1,1};
    morph::vec<float, 2> domain_offs = {0,0};

    //! Initialize a rectangular net of width w and height h. Note that this resizes p
    //! and arranges c, but does not fill p with positions.
    void init (size_t _w, size_t _h)
    {
        this->w = _w;
        this->h = _h;
        this->p.resize(_w*_h);
        this->targ.resize(_w*_h);
        this->clr.resize(_w*_h);
        // Encode colours in the net, using a ColourMap for help
        morph::ColourMap<float> cmap;
        if (use_hsv_cmaps == true) {
            cmap.setType (morph::ColourMapType::HSV);
            cmap.setHueRotation (-morph::mathconst<float>::pi);
        } else {
            cmap.setType (morph::ColourMapType::Duochrome);
            cmap.setHueRG();
        }
        size_t i = 0;
        for (size_t y = 0; y < _h; ++y) {
            for (size_t x = 0; x < _w; ++x) {
                this->clr[i++] = cmap.convert ((float)x/(float)w, (float)y/(float)h);
            }
        }
        // Set up connections
        for (size_t y = 0; y < _h-1; ++y) {
            for (size_t x = 0; x < _w; ++x) {
                // Connect down
                this->c.insert (morph::vec<size_t, 2>({x+y*_w, x+(y+1)*_w}));
            }
        }
        for (size_t x = 0; x < _w-1; ++x) {
            for (size_t y = 0; y < _h; ++y) {
                // Connect right
                this->c.insert (morph::vec<size_t, 2>({x+y*_w, 1+x+y*_w}));
            }
        }
    }

    static constexpr bool debug_crossing = false;
    // Count the number of crossings of the links in the net using naive algorithm (O(n^2) time)
    // Algorithm taken from
    // https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
    size_t crosscount()
    {
        size_t cc = 0;

        // Have to iterate the std::set
        typename std::set<morph::vec<size_t, 2>, ncmp>::iterator ci = this->c.begin();

        typename std::set< std::set<morph::vec<size_t, 2>, ncmp>> pairs_tested;

        // Can this be OpenMPed?
        while (ci != this->c.end()) {
            typename std::set<morph::vec<size_t, 2>, ncmp>::iterator cj = this->c.begin();
            while (cj != this->c.end()) {

                if (cj == ci) { // then skip
                    ++cj;
                    continue;
                }

                morph::vec<float, 2> p1 = this->p[(*ci)[0]].less_one_dim();
                morph::vec<float, 2> q1 = this->p[(*ci)[1]].less_one_dim();
                morph::vec<float, 2> p2 = this->p[(*cj)[0]].less_one_dim();
                morph::vec<float, 2> q2 = this->p[(*cj)[1]].less_one_dim();

                if (p1 == p2 || p1 == q2 || q1 == p2 || q1 == q2) {
                    // Don't count intersections between line segments that are connected in the net ANYWAY
                    ++cj;
                    continue;
                }

                // Have we already considered ci, cj as cj, ci?
                std::set<morph::vec<size_t, 2>, ncmp> linepair;
                linepair.insert (*ci);
                linepair.insert (*cj);
                if (pairs_tested.count (linepair)) { // skip this pair as we already considered it in another order
                    ++cj;
                    continue;
                }
                pairs_tested.insert (linepair);

                morph::rotation_sense p1q1p2 = morph::MathAlgo::orientation (p1,q1,p2);
                morph::rotation_sense p1q1q2 = morph::MathAlgo::orientation (p1,q1,q2);
                morph::rotation_sense p2q2p1 = morph::MathAlgo::orientation (p2,q2,p1);
                morph::rotation_sense p2q2q1 = morph::MathAlgo::orientation (p2,q2,q1);

                if (p1q1p2 != p1q1q2 && p2q2p1 != p2q2q1) {
                    // ci and cj intersect
                    ++cc;
                } else {
                    // Are they colinear?
                    if (p1q1p2 == morph::rotation_sense::colinear && morph::MathAlgo::onsegment (p1, p2, q1)) { ++cc; }
                    else if (p1q1q2 == morph::rotation_sense::colinear && morph::MathAlgo::onsegment (p1, q2, q1)) { ++cc; }
                    else if (p2q2p1 == morph::rotation_sense::colinear && morph::MathAlgo::onsegment (p2, p1, q2)) { ++cc; }
                    else if (p2q2q1 == morph::rotation_sense::colinear && morph::MathAlgo::onsegment (p2, q1, q2)) { ++cc; }
                }
                ++cj;
            }
            ++ci;
        }

        return cc;
    }

    //! Return sum of squared distances between p and targ
    T sos() const { return ((p-targ).sos()[0]); }

    //! Return sum of squared distances between p and targ for a patch inside the (2D) square defined by 2 corners
    T sos_inside (const morph::vec<T, 3>& corner1, const morph::vec<T, 3>& corner2) const
    {
        morph::vvec<T> mask (p.size());
        for (unsigned int i = 0; i < p.size(); ++i) {
            // Mask set to 1 for targets inside the corners
            mask[i] = (targ[i][0] >= corner1[0]
                       && targ[i][1] >= corner1[1]
                       && targ[i][0] <= corner2[0]
                       && targ[i][1] <= corner2[1]) ? T{1} : T{0};
        }
        morph::vvec<morph::vec<T, 3>> p_minus_targ = (p - targ) * mask;
        return (p_minus_targ.sos()[0]);
    }

    //! Return sum of squared distances between p and targ for a patch outside the (2D) square defined by 2 corners
    T sos_outside (const morph::vec<T, 3>& corner1, const morph::vec<T, 3>& corner2) const
    {
        morph::vvec<T> mask (p.size());
        for (unsigned int i = 0; i < p.size(); ++i) {
            // Mask set to 1 for targets outside the corners
            mask[i] = (targ[i][0] < corner1[0]
                       || targ[i][1] < corner1[1]
                       || targ[i][0] > corner2[0]
                       || targ[i][1] > corner2[1]) ? T{1} : T{0};
        }
        morph::vvec<morph::vec<T, 3>> p_minus_targ = (p - targ) * mask;
        return (p_minus_targ.sos()[0]);
    }

    //! Return RMS error of each agent from its target
    T rms() const { return std::sqrt(((p-targ).sos()[0]) / this->targ.size()); }

    //! Return mean distances between p and targ for a patch inside the (2D) square defined by 2 corners
    T rms_inside (const morph::vec<T, 3>& corner1, const morph::vec<T, 3>& corner2) const
    {
        morph::vvec<T> mask (p.size());
        for (unsigned int i = 0; i < p.size(); ++i) {
            // Mask set to 1 for targets inside the corners
            mask[i] = (targ[i][0] >= corner1[0]
                       && targ[i][1] >= corner1[1]
                       && targ[i][0] <= corner2[0]
                       && targ[i][1] <= corner2[1]) ? T{1} : T{0};
        }
        morph::vvec<morph::vec<T, 3>> p_minus_targ = (p - targ) * mask;
        return std::sqrt((p_minus_targ.sos()[0])/mask.sum());
    }

    //! Return mean distances between p and targ for a patch outside the (2D) square defined by 2 corners
    T rms_outside (const morph::vec<T, 3>& corner1, const morph::vec<T, 3>& corner2) const
    {
        morph::vvec<T> mask (p.size());
        for (unsigned int i = 0; i < p.size(); ++i) {
            // Mask set to 1 for targets outside the corners
            mask[i] = (targ[i][0] < corner1[0]
                       || targ[i][1] < corner1[1]
                       || targ[i][0] > corner2[0]
                       || targ[i][1] > corner2[1]) ? T{1} : T{0};
        }
        morph::vvec<morph::vec<T, 3>> p_minus_targ = (p - targ) * mask;
        return std::sqrt((p_minus_targ.sos()[0])/mask.sum());
    }

    //! Positions of the vertices of the net
    morph::vvec<morph::vec<T, 3>> p;
    //! Colours of the vertices of the net
    std::vector<std::array<float, 3>> clr;

    //! Connections of the net. The indices into p that are the ends of line segments
    std::set<morph::vec<size_t, 2>, ncmp> c;

    //! Target positions of the net, to allow the net to be used to provide a SOS metric for how close the pattern is to ideal
    morph::vvec<morph::vec<T, 3>> targ;
};

// A retinal ganglion cell net specialised to know about graft rotations/swaps etc
template<typename T>
struct rgcnet : public net<T>
{
    // Set true if the x/y coordinates in targ have been mirrrored
    bool mirrored = false;

    //! Expand the top half down across the full tissue area
    void targ_expand_topdown() { for (auto& t : this->targ) { t[1] = t[1] * T{2} - T{1}; } }
    //! Expand bottom half up across the full tissue area
    void targ_expand_bottomup() { for (auto& t : this->targ) { t[1] = t[1] * T{2}; } }

    //! Tectum has reduced in size so squish the targs down
    void targ_squish_topdown() { for (auto& t : this->targ) { t[1] = t[1] * T{0.5}; } }
    //! Tectum has reduced in size so squish the targs up
    void targ_squish_bottomup() { for (auto& t : this->targ) { t[1] = t[1] * T{0.5} + T{0.5}; } }

    //! Apply a graft rotation to the targs
    void targ_graftrotate (morph::vec<size_t,2> x1, size_t sz, size_t quarter_rotations)
    {
        if (quarter_rotations%4 == 0) { return; }

        // Potential for runtime crash here, so check x1 and sz.
        if (x1[0] + sz > this->w || x1[1] + sz > this->h) {
            throw std::runtime_error ("Patch for rotation is off the edge of the net");
        }

        morph::vvec<morph::vec<T,3>> targ_graft (sz*sz);

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
    void targ_graftswap (morph::vec<size_t,2> x1,
                         morph::vec<size_t,2> sz,
                         morph::vec<size_t,2> x2)
    {
        // I copied this from tissue, which does the swap by index and so I need to
        // cover two situations; one where I swap row by row, and one where, if
        // mirrored, I swap col by col.
        if (mirrored == true) {
            // move transformed section of targ row by row.
            morph::vvec<morph::vec<T,3>> graft (sz[1]);

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
            morph::vvec<morph::vec<T,3>> graft (sz[0]);

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

    // NOT the same as tissue::compound tissue. For the *targets* we expect two overlaid nets.
    void targ_compound_tissue (tissue_region retina_mirrored = tissue_region::left_half)
    {
        morph::vec<T,3> fac2y = {1, 2, 1};
        morph::vec<T,3> facm2y = {1, -2, 1};
        morph::vec<T,3> offstwoy = {0, 2, 0};
        switch (retina_mirrored) {
        case tissue_region::right_half:
        case tissue_region::top_half:
        case tissue_region::bottom_half:
            std::cout << "implement\n";
            break;
        case tissue_region::left_half:
        default:
            // Retina left was copied to right. So See overlay on the y axis of targ

            // Bottom half has locations stretched up
            for (size_t j = 0; j < this->h/2; ++j) {
                for (size_t i = 0; i < this->w; ++i) {
                    this->targ[j+this->w*i] *= fac2y;
                }
            }
            // Top part is folded over and stretched out!
            for (size_t j = this->h/2; j < this->h; ++j) {
                for (size_t i = 0; i < this->w; ++i) {
                    this->targ[j+this->w*i] = offstwoy + (this->targ[j+this->w*i] * facm2y);
                    //this->targ[j+this->w*i] = oney - this->targ[j+this->w*i];
                }
            }
            break;
        }
    }

    // Set up the expected experimental layout for the Reber or Brown manipulations. The
    // genetic manipulations are carried out by
    // guidingtissue::receptor_knockin/out/down. I'm not sure this is actually used.
    void targ_reber()
    {
        morph::vec<T,3> fac2y = {1, 2, 1};
        morph::vec<T,3> offshalfy = {0, 0.5, 0};

        // Reber expt: knockdown one EphR, selectively knockin another. Both are equivalent to our rcpt[0].
        // Brown expt: Just do the selective knockin on its own.

        // Bottom half has locations stretched down
        for (size_t j = 0; j < this->h/2; ++j) {
            for (size_t i = 0; i < this->w; ++i) {
                this->targ[j+this->w*i] *= fac2y;
            }
        }
        // Top part is shifted/stretched up
        for (size_t j = this->h/2; j < this->h; ++j) {
            for (size_t i = 0; i < this->w; ++i) {
                this->targ[j+this->w*i] -= offshalfy;
                this->targ[j+this->w*i] *= fac2y;
            }
        }
    }
};
