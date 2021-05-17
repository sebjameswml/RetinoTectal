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
    // Set false to get linear receptor and ligand expression; set true to get exponential curves
    bool exp_expression = true;

    //! A left-right interaction parameter and an up-down interaction parameter for each
    //! piece of guiding tissue requires 4 receptors and 4 ligands. Holds receptor
    //! expressions for each cell.
    morph::vVector<morph::Vector<T,N>> rcpt;
    morph::vVector<morph::Vector<T,N>> rcpt_grad;

    //! The tissue also has an expression of ligands to interact with receptors of other
    //! cells/agents
    morph::vVector<morph::Vector<T,N>> lgnd;
    morph::vVector<morph::Vector<T,N>> lgnd_grad;

    //! Init the wildtype tissue's receptors and ligands.
    guidingtissue(size_t _w, size_t _h, morph::Vector<T,2> _dx, morph::Vector<T,2> _x0, bool ee = true)
        : tissue<T> (_w, _h, _dx, _x0)
        , exp_expression (ee)
    {
        this->rcpt.resize (this->posn.size());
        this->lgnd.resize (this->posn.size());
        this->rcpt_grad.resize (this->posn.size());
        this->lgnd_grad.resize (this->posn.size());

        T max_x = this->w * this->dx[0] + this->x0[0];
        T max_y = this->h * this->dx[1] + this->x0[1];
        // A cout avoids some unused variable compiler warnings:
        std::cout << "Max x position: " << max_x << " and max y position: " << max_y << std::endl;

        // Ligands and receptors are set up as a function of their cell's position in the tissue.
        for (size_t ri = 0; ri < this->posn.size(); ++ri) {
            if constexpr (N == 4) {
                // First orthogonal pair of receptors
                this->rcpt[ri][0] = this->expression_function (this->posn[ri][0]);
                this->rcpt[ri][1] = this->expression_function (this->posn[ri][1]);
                this->rcpt_grad[ri][0] = this->gradient_function (this->posn[ri][0]);
                this->rcpt_grad[ri][1] = this->gradient_function (this->posn[ri][1]);
                // First orthogonal pair of ligands have same expression as receptors
                this->lgnd[ri][0] = this->rcpt[ri][0];
                this->lgnd[ri][1] = this->rcpt[ri][1];
                this->lgnd_grad[ri][0] = this->rcpt_grad[ri][0];
                this->lgnd_grad[ri][1] = this->rcpt_grad[ri][1];

                // Second orthogonal pair in opposite sense
                this->rcpt[ri][2] = this->expression_function (max_x-this->posn[ri][0]);
                this->rcpt[ri][3] = this->expression_function (max_y-this->posn[ri][1]);
                this->rcpt_grad[ri][2] = this->gradient_function (max_x-this->posn[ri][0]);
                this->rcpt_grad[ri][3] = this->gradient_function (max_y-this->posn[ri][1]);
                this->lgnd[ri][2] = this->rcpt[ri][2];
                this->lgnd[ri][3] = this->rcpt[ri][3];
                this->lgnd_grad[ri][2] = this->rcpt_grad[ri][2];
                this->lgnd_grad[ri][3] = this->rcpt_grad[ri][3];

            } else if constexpr (N == 2) {
                // Just one orthogonal pair of receptors and ligands:
                this->rcpt[ri][0] = this->expression_function (this->posn[ri][0]);
                this->rcpt[ri][1] = this->expression_function (this->posn[ri][1]);
                this->rcpt_grad[ri][0] = this->gradient_function (this->posn[ri][0]);
                this->rcpt_grad[ri][1] = this->gradient_function (this->posn[ri][1]);
                this->lgnd[ri][0] = this->rcpt[ri][0];
                this->lgnd[ri][1] = this->rcpt[ri][1];
                this->lgnd_grad[ri][0] = this->rcpt[ri][0];
                this->lgnd_grad[ri][1] = this->rcpt_grad[ri][1];
            } // else error
        }

        // Now build up the tissue's enclosing border - imagine a border of signalling
        // molecules that push axons back into the region in which they're supposed to
        // be moving. Or perhaps, I should have a branching scheme where the axons
        // branch more in the middle of the tissue and less at the periphery?  Problem
        // with this idea is that the tissue is only the square. So, either I make the
        // tissue larger, or go with an algorithm based thing in branch.h
    }

    // R(x) = 1.31 + 2.3333x gives a gradient of 2.3333
    T linear_expression (const T& x) const { return T{1.31} + T{2.3333} * x; }
    T linear_gradient (const T& x) const { return T{2.3333}; }
    // R(x) = 0.26e^(2.3x) + 1.05,
    T exponential_expression (const T& x) const { return T{1.05} + T{0.26} * std::exp (T{2.3} * x); }
    // 0.26*2.3 e^(2.3x)
    T exponential_gradient (const T& x) const { return T{0.598} * std::exp (T{2.3} * x); }

    T expression_function (const T& x) const
    {
        if (this->exp_expression) {
            return this->exponential_expression (x);
        } else {
            return this->linear_expression (x);
        }
    }
    T gradient_function (const T& x) const
    {
        if (this->exp_expression) {
            return this->exponential_gradient (x);
        } else {
            return this->linear_gradient (x);
        }
    }

    //! With the passed-in location, find the closest gradient in lgnd_grad and return
    //! this.  Note: Not a function (like linear_gradient()) but a lookup, because this
    //! makes it easy to perform tissue manipulations on target (tectal) tissue in the
    //! same way as for source (retinal) tissue. This can also potentially return border
    //! gradients for regions outside the nominal tissue area.
    morph::Vector<T,N> lgnd_grad_at (const morph::Vector<T,2>& x) const
    {
        morph::vVector<morph::Vector<T,2>> p = this->posn - x;
        size_t idx = p.argshortest();
        return this->lgnd_grad[idx];
    }

    //! Move a rectangular patch of tissue at indexed location x1, of size
    //! sz, and swap it with another region of the tissue at x2.
    void graftswap (morph::Vector<size_t,N> x1,
                    morph::Vector<size_t,N> sz,
                    morph::Vector<size_t,N> x2)
    {
        // This is all about moving a section of rcpt. Have to move row by row though.
        morph::vVector<morph::Vector<T,N>> graft (sz[0]);

        size_t idx1 = x1[0] + this->w * x1[1];
        size_t idx2 = x2[0] + this->w * x2[1];
        for (size_t i = 0; i < sz[1]; ++i) {

            if (idx1 + sz[0] > this->rcpt.size() || idx2 + sz[0] > this->rcpt.size()) {
                throw std::runtime_error ("Can't make that patchswap");
            }

            // Swap rcpt and lgnd parameters
            auto oi = this->rcpt.begin();
            oi += idx1;
            auto ni = this->rcpt.begin();
            ni += idx2;
            std::copy_n (ni, sz[0], graft.begin());
            std::copy_n (oi, sz[0], ni);
            std::copy_n (graft.begin(), sz[0], oi);

            oi = this->lgnd.begin();
            oi += idx1;
            ni = this->lgnd.begin();
            ni += idx2;
            std::copy_n (ni, sz[0], graft.begin());
            std::copy_n (oi, sz[0], ni);
            std::copy_n (graft.begin(), sz[0], oi);

            oi = this->rcpt_grad.begin();
            oi += idx1;
            ni = this->rcpt_grad.begin();
            ni += idx2;
            std::copy_n (ni, sz[0], graft.begin());
            std::copy_n (oi, sz[0], ni);
            std::copy_n (graft.begin(), sz[0], oi);

            oi = this->lgnd_grad.begin();
            oi += idx1;
            ni = this->lgnd_grad.begin();
            ni += idx2;
            std::copy_n (ni, sz[0], graft.begin());
            std::copy_n (oi, sz[0], ni);
            std::copy_n (graft.begin(), sz[0], oi);

            idx1 += this->w;
            idx2 += this->w;
        }
    }
};
