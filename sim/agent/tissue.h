#pragma once

#include <morph/vVector.h>
#include <morph/Vector.h>
#include <algorithm>

//! Some 2D brain tissue, arranged as a cartesian grid of neurons of width w and height h.
template<typename T>
struct tissue
{
    size_t w;
    size_t h;
    const size_t num() const { return w * h; }
    //! The distance from neuron to neuron; horizontally and vertically.
    morph::Vector<T,2> dx;
    //! Position of bottom left neuron
    morph::Vector<T,2> x0;
    //! Positions stored in raster fashion from bottom left
    morph::vVector<morph::Vector<T,2>> posn;

    tissue(size_t _w, size_t _h, morph::Vector<T,2> _dx, morph::Vector<T,2> _x0)
        : w(_w), h(_h), dx(_dx), x0(_x0)
    {
        this->posn.resize(w*h);
        for (size_t i = 0; i < this->h; ++i) {
            for (size_t j = 0; j < this->w; ++j) {
                posn[i*w+j] = {dx[0]*j, dx[1]*i};
            }
        }
    }
};

//! Tissue which has guidance parameters in the form of a 2D vector which stands for the
//! expression of 2 or possibly 4 receptor molecules. If positive/negative values are
//! allowed, then it's 4.
template<typename T, size_t N>
struct guidingtissue : public tissue<T>
{
    //! A left-right interaction parameter and an up-down interaction parameter for each
    //! piece of guiding tissue requires 4 receptors and 4 ligands. Holds receptor
    //! expressions for each cell.
    morph::vVector<morph::Vector<T,N>> rcpt;

    //! The tissue also has an expression of ligands to interact with receptors of other
    //! cells/agents
    morph::vVector<morph::Vector<T,N>> lgnd;

    //! Init the wildtype tissue's receptors and ligands.
    guidingtissue(size_t _w, size_t _h, morph::Vector<T,2> _dx, morph::Vector<T,2> _x0)
        : tissue<T> (_w, _h, _dx, _x0)
    {
        this->rcpt.resize (this->posn.size());
        this->lgnd.resize (this->posn.size());

        T max_x = this->w * this->dx[0] + this->x0[0];
        T max_y = this->h * this->dx[1] + this->x0[1];

        // Ligands and receptors are set up as a function of their cell's position in the tissue.
        for (size_t ri = 0; ri < this->posn.size(); ++ri) {
            if constexpr (N == 4) {
                // First orthogonal pair of receptors
                this->rcpt[ri][0] = T{1.05} + (T{0.26} * std::exp (T{2.3} * this->posn[ri][0])); // R(x) = 0.26e^(2.3x) + 1.05,
                this->rcpt[ri][1] = T{1.05} + (T{0.26} * std::exp (T{2.3} * this->posn[ri][1]));
                // First orthogonal pair of ligands
                this->lgnd[ri][0] = T{1.05} + (T{0.26} * std::exp (T{2.3} * this->posn[ri][0]));
                this->lgnd[ri][1] = T{1.05} + (T{0.26} * std::exp (T{2.3} * this->posn[ri][1]));

                // Second orthogonal pair in opposite sense
                this->rcpt[ri][2] = T{1.05} + (T{0.26} * std::exp (T{2.3} * max_x-this->posn[ri][0]));
                this->rcpt[ri][3] = T{1.05} + (T{0.26} * std::exp (T{2.3} * max_y-this->posn[ri][1]));

                this->lgnd[ri][2] = T{1.05} + (T{0.26} * std::exp (T{2.3} * max_x-this->posn[ri][0]));
                this->lgnd[ri][3] = T{1.05} + (T{0.26} * std::exp (T{2.3} * max_y-this->posn[ri][1]));

            } else if constexpr (N == 2) {
                // Just one orthogonal pair of receptors and ligands:
                this->rcpt[ri][0] = T{1.05} + (T{0.26} * std::exp (T{2.3} * this->posn[ri][0]));
                this->rcpt[ri][1] = T{1.05} + (T{0.26} * std::exp (T{2.3} * this->posn[ri][1]));
                this->lgnd[ri][0] = T{1.05} + (T{0.26} * std::exp (T{2.3} * this->posn[ri][0]));
                this->lgnd[ri][1] = T{1.05} + (T{0.26} * std::exp (T{2.3} * this->posn[ri][1]));
            } // else error
        }
    }

    //! Move a rectangular patch of tissue at indexed location locn1, of size
    //! sz, and swap it with another region of the tissue at locn2.
    void graftswap (morph::Vector<size_t,N> locn1,
                    morph::Vector<size_t,N> sz,
                    morph::Vector<size_t,N> locn2)
    {
        // This is all about moving a section of rcpt. Have to move row by row though.
        morph::vVector<morph::Vector<T,N>> tmp_rcpt; // Temporarily holds a section of receptors

        tmp_rcpt.resize (sz[0]);

        size_t idx1 = locn1[0] + this->w * locn1[1];
        size_t idx2 = locn2[0] + this->w * locn2[1];
        for (size_t i = 0; i < sz[1]; ++i) {

            if (idx1 + sz[0] > this->rcpt.size()
                || idx2 + sz[0] > this->rcpt.size()) {
                throw std::runtime_error ("Can't make that patchswap");
            }

            // Swap rcpt parameters
            auto oi = this->rcpt.begin();
            oi += idx1;
            auto ni = this->rcpt.begin();
            ni += idx2;
            std::copy_n (ni, sz[0], tmp_rcpt.begin());
            std::copy_n (oi, sz[0], ni);
            std::copy_n (tmp_rcpt.begin(), sz[0], oi);

            idx1 += this->w;
            idx2 += this->w;
        }
    }
};
