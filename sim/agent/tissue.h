#pragma once

#include <morph/vVector.h>
#include <morph/Vector.h>
#include <algorithm>

//! Some 2D brain tissue, arranged as a Cartesian grid of neurons of width w and height h.
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
    //! Each receptor has a 2D gradient field, hence 2*N values here
    morph::vVector<morph::Vector<T,2*N>> rcpt_grad;

    //! The tissue also has an expression of ligands to interact with receptors of other
    //! cells/agents
    morph::vVector<morph::Vector<T,N>> lgnd;
    morph::vVector<morph::Vector<T,2*N>> lgnd_grad;

    static constexpr bool debug_gradf = false;

    void spacegrad2D (morph::vVector<morph::Vector<T,N>>& f,
                      morph::vVector<morph::Vector<T, (2*N)>>& gradf)
    {
        if constexpr (debug_gradf) { std::cout << "Main block\n"; }
        //#pragma omp parallel for schedule(static)
        for (size_t x = 0; x < this->w-1; ++x) {
            for (size_t y = 0; y < this->h-1; ++y) {
                size_t i = x + y*this->w;
                for (size_t n = 0; n < N; n++) {
                    size_t nn = 2*n;
                    gradf[i][nn] = (f[i+1][n] - f[i][n]) / this->dx[0];
                    gradf[i][nn+1] = (f[i+this->w][n] - f[i][n]) / this->dx[1];
                }
                if constexpr (debug_gradf) {
                    std::cout << "("<<x<<","<<y<<") gradf=" << gradf[i] << std::endl;
                }
            }
        }
        // Last row is set to same gradient as penultimate row
        if constexpr (debug_gradf) { std::cout << "Last row\n"; }
        for (size_t x = 0; x < this->w-1; ++x) {
            size_t i = x + (this->h-1)*this->w;
            if constexpr (debug_gradf) {
                std::cout << "x=" << x << ", set gradf[" << i << "] to value in gradf[" << (i-this->w) << "]"<< std::endl;
            }
            gradf[i] = gradf[i-this->w];
        }
        // Last col
        if constexpr (debug_gradf) { std::cout << "Last col\n"; }
        for (size_t y = 0; y < this->h; ++y) {
            size_t i = (this->w-1) + (y*this->w);
            if constexpr (debug_gradf) {
                std::cout << "y=" << y << ", set gradf[" << i << "] to value in gradf[" << (i-1) << "] = " << gradf[i-1] << std::endl;
            }
            gradf[i] = gradf[i-1];
        }
    }

    //! Numerically differentiate rcpt and lgnd
    void compute_gradients()
    {
        this->rcpt_grad.resize (this->rcpt.size());
        this->lgnd_grad.resize (this->lgnd.size());
        this->spacegrad2D (this->rcpt, this->rcpt_grad);
        this->spacegrad2D (this->lgnd, this->lgnd_grad);
    }

    //! Init the wildtype tissue's receptors and ligands.
    guidingtissue(size_t _w, size_t _h, morph::Vector<T,2> _dx, morph::Vector<T,2> _x0, bool ee = true)
        : tissue<T> (_w, _h, _dx, _x0)
        , exp_expression (ee)
    {
        std::cout << "This guiding tissue has "
                  << (this->exp_expression ? "exponential" : "linear") << " expression" << std::endl;
        this->rcpt.resize (this->posn.size());
        this->lgnd.resize (this->posn.size());

        T max_x = this->w * this->dx[0] + this->x0[0];
        T max_y = this->h * this->dx[1] + this->x0[1];
        // A cout avoids some 'unused variable' compiler warnings:
        std::cout << "Max x position: " << max_x << " and max y position: " << max_y << std::endl;

        // Ligands and receptors are set up as a function of their cell's position in the tissue.
        for (size_t ri = 0; ri < this->posn.size(); ++ri) {
            if constexpr (N == 4) {
                // First orthogonal pair of receptors
                this->rcpt[ri][0] = this->expression_function (this->posn[ri][0]);
                this->rcpt[ri][1] = this->expression_function (this->posn[ri][1]);
                // First orthogonal pair of ligands have same expression as receptors
                this->lgnd[ri][0] = this->rcpt[ri][0];
                this->lgnd[ri][1] = this->rcpt[ri][1];

                // Second orthogonal pair in opposite sense
                this->rcpt[ri][2] = this->expression_function (max_x-this->posn[ri][0]);
                this->rcpt[ri][3] = this->expression_function (max_y-this->posn[ri][1]);
                this->lgnd[ri][2] = this->rcpt[ri][2];
                this->lgnd[ri][3] = this->rcpt[ri][3];

            } else if constexpr (N == 2) {
                // Just one orthogonal pair of receptors and ligands:
                this->rcpt[ri][0] = this->expression_function (this->posn[ri][0]);
                this->rcpt[ri][1] = this->expression_function (this->posn[ri][1]);
                this->lgnd[ri][0] = this->rcpt[ri][0];
                this->lgnd[ri][1] = this->rcpt[ri][1];
            } // else error
        }

        this->compute_gradients();

        // Now build up the tissue's enclosing border - imagine a border of signalling
        // molecules that push axons back into the region in which they're supposed to
        // be moving. Or perhaps, I should have a branching scheme where the axons
        // branch more in the middle of the tissue and less at the periphery?  Problem
        // with this idea is that the tissue is only the square. So, either I make the
        // tissue larger, or go with an algorithm based thing in branch.h
    }

    // R(x) = 1.31 + 2.3333x gives a gradient of 2.3333
    T linear_expression (const T& x) const { return T{1.31} + T{2.3333} * x; }
    T exponential_expression (const T& x) const { return T{1.05} + T{0.26} * std::exp (T{2.3} * x); }

    T expression_function (const T& x) const
    {
        if (this->exp_expression) {
            return this->exponential_expression (x);
        } else {
            return this->linear_expression (x);
        }
    }

    //! With the passed-in location, find the closest gradient in lgnd_grad and return
    //! this.  Note: Not a function (like linear_gradient()) but a lookup, because this
    //! makes it easy to perform tissue manipulations on target (tectal) tissue in the
    //! same way as for source (retinal) tissue. This can also potentially return border
    //! gradients for regions outside the nominal tissue area.
    morph::Vector<T,2*N> lgnd_grad_at (const morph::Vector<T,2>& x) const
    {
        morph::vVector<morph::Vector<T,2>> p = this->posn - x;
        size_t idx = p.argshortest();
        return this->lgnd_grad[idx];
    }

    //! Move a rectangular patch of tissue at indexed location x1, of size
    //! sz, and swap it with another region of the tissue at x2.
    void graftswap (morph::Vector<size_t,2> x1,
                    morph::Vector<size_t,2> sz,
                    morph::Vector<size_t,2> x2)
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

            idx1 += this->w;
            idx2 += this->w;
        }

        this->compute_gradients();
    }
};
