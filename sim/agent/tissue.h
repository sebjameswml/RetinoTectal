#pragma once

#include <morph/vVector.h>
#include <morph/Vector.h>
#include <morph/Random.h>
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

    const T x_min() const { return this->posn[0][0]; }
    const T x_max() const { return this->posn[this->posn.size()-1][0]; }
    const T y_min() const { return this->posn[0][1]; }
    const T y_max() const { return this->posn[this->posn.size()-1][1]; }

    // Find the coordinate closest to the indexed position xx, yy and return
    morph::Vector<T,2> coord (size_t xx, size_t yy) const
    {
        morph::Vector<T,2> c = {0,0};

        // The x/y position to find has to be constructed from the index args xx and yy
        morph::Vector<T,2> xy = { x0[0] + xx * dx[0], x0[1] + yy * dx[1] };

        // Now find the first posn in the tissue grid which is close to xy
        for (auto p : this->posn) {
            if ((xy-p).length() < dx.shortest()/T{2}) {
                c = p;
                break;
            }
        }

        return c;
    }

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

enum class expression_form
{
    lin,  // A linear function
    quad, // A quadratic function
    exp,  // An exponential function
    log   // A logarithmic function
};

// Which sense is an expression pattern becoming stronger?
enum class expression_direction
{
    increasing, // More expression as x or y increases
    decreasing  // Less expression as x or y increases
};

enum class tissue_region
{
    top_half,    // Refers to the top part of the tissue (from y_max/2 to y_max)
    bottom_half, // The bottom part of the tissue (0 to y_max/2)
    left_half,   // The left half of the tissue (from 0 to x_max/2)
    right_half   // The right half
};

// What interaction does a forward signal transmit when ligand binds to receptor?  This
// is either attraction - the branch climbs the ligand gradient - or repulsion - the
// branch descends the ligand gradient.
enum class interaction
{
    attraction,
    repulsion
};

/*!
 * Tissue which has guidance parameters in the form of a 2D vector which stands for the
 * expression of 2 or possibly 4 receptor molecules. If positive/negative values are
 * allowed, then it's 4.
 *
 * tparam T: Data type for processing. Usually float unless double seems necessary.
 *
 * tparam N: There are N receptors and N ligands in the guidingtissue.
 */
template<typename T, size_t N>
struct guidingtissue : public tissue<T>
{
    //! With what kind of function are receptors expressed?
    morph::Vector<expression_form, N> rcpt_form;
    //! With what kind of function are ligands expressed?
    morph::Vector<expression_form, N> lgnd_form;
    //! Directions of receptor expression gradients
    morph::Vector<expression_direction, N> rcpt_dirns;
    //! Directions of ligand expression gradients
    morph::Vector<expression_direction, N> lgnd_dirns;

    // Forward signalling interactions when receptors in this tissue are activated
    morph::Vector<interaction, N> forward_interactions;
    // Reverse signalling interactions when ligands in this tissue are activated
    morph::Vector<interaction, N> reverse_interactions;

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

    //! Init the wildtype tissue's receptors and ligands.
    guidingtissue(size_t _w, size_t _h,
                  morph::Vector<T,2> _dx, morph::Vector<T,2> _x0,
                  morph::Vector<expression_form, N> _rcpt_form,
                  morph::Vector<expression_form, N> _lgnd_form,
                  morph::Vector<expression_direction, N> _rcpt_dirn,
                  morph::Vector<expression_direction, N> _lgnd_dirn,
                  morph::Vector<interaction, N> _for_int,
                  morph::Vector<interaction, N> _rev_int)
        : tissue<T> (_w, _h, _dx, _x0)
        , rcpt_form(_rcpt_form)
        , lgnd_form(_lgnd_form)
        , rcpt_dirns(_rcpt_dirn)
        , lgnd_dirns(_lgnd_dirn)
        , forward_interactions(_for_int)
        , reverse_interactions(_rev_int)
    {
        this->rcpt.resize (this->posn.size());
        this->lgnd.resize (this->posn.size());

        T max_x = this->w * this->dx[0] + this->x0[0];
        T max_y = this->h * this->dx[1] + this->x0[1];
        // A cout avoids some 'unused variable' compiler warnings:
        std::cout << "Max x position: " << max_x << " and max y position: " << max_y << std::endl;

        // Ligands and receptors are set up as a function of their cell's position in the tissue.
        for (size_t ri = 0; ri < this->posn.size(); ++ri) {
            T xin = T{0};
            T yin = T{0};
            if constexpr (N == 4 || N == 2) {
                // First orthogonal pair of receptors.
                xin = (this->rcpt_dirns[0] == expression_direction::increasing) ? this->posn[ri][0] : (max_x-this->posn[ri][0]);
                this->rcpt[ri][0] = this->rcpt_expression_function (xin, 0);
                yin = (this->rcpt_dirns[1] == expression_direction::increasing) ? this->posn[ri][1] : (max_y-this->posn[ri][1]);
                this->rcpt[ri][1] = this->rcpt_expression_function (yin, 1);
                // First orthogonal pair of ligands
                xin = (this->lgnd_dirns[0] == expression_direction::increasing) ? this->posn[ri][0] : (max_x-this->posn[ri][0]);
                this->lgnd[ri][0] = this->lgnd_expression_function (xin, 0);
                yin = (this->lgnd_dirns[1] == expression_direction::increasing) ? this->posn[ri][1] : (max_y-this->posn[ri][1]);
                this->lgnd[ri][1] = this->lgnd_expression_function (yin, 1);
            } else {
                // C++-20 mechanism to trigger a compiler error for the else case. Not user friendly!
                []<bool flag = false>() { static_assert(flag, "no match"); }();
            }
            if constexpr (N == 4) {
                // Add a second orthogonal pair for N==4
                xin = (this->rcpt_dirns[2] == expression_direction::increasing) ? this->posn[ri][0] : (max_x-this->posn[ri][0]);
                this->rcpt[ri][2] = this->rcpt_expression_function (xin, 2);
                yin = (this->rcpt_dirns[3] == expression_direction::increasing) ? this->posn[ri][1] : (max_y-this->posn[ri][1]);
                this->rcpt[ri][3] = this->rcpt_expression_function (yin, 3);

                // Second orthogonal pair of ligands
                xin = (this->lgnd_dirns[2] == expression_direction::increasing) ? this->posn[ri][0] : (max_x-this->posn[ri][0]);
                this->lgnd[ri][2] = this->lgnd_expression_function (xin, 2);
                yin = (this->lgnd_dirns[3] == expression_direction::increasing) ? this->posn[ri][1] : (max_y-this->posn[ri][1]);
                this->lgnd[ri][3] = this->lgnd_expression_function (yin, 3);
            }
        }

        this->compute_gradients();

        // Now build up the tissue's enclosing border - imagine a border of signalling
        // molecules that push axons back into the region in which they're supposed to
        // be moving. Or perhaps, I should have a branching scheme where the axons
        // branch more in the middle of the tissue and less at the periphery?  Problem
        // with this idea is that the tissue is only the square. So, either I make the
        // tissue larger, or go with an algorithm based thing in branch.h
    }

    //! Do you want to debug spacegrad2D?
    static constexpr bool debug_gradf = false;
    /*!
     * Computes the gradient at an element from the expression in the neighbouring elements
     *
     * \param f Variable vector (vVector) containing Vectors which encode N scalar fields.
     *
     * \param gradf Output. Each Vector in gradf has twice as many fields as the Vectors
     * in f. The gradients are:
     *
     *  gradf[0]: f[0] x gradient          gradf[1]: f[0] y gradient
     *  gradf[2]: f[1] x gradient          gradf[3]: f[1] y gradient
     *  ... etc
     */
    void spacegrad2D (morph::vVector<morph::Vector<T,N>>& f,
                      morph::vVector<morph::Vector<T,(2*N)>>& gradf)
    {
        //#pragma omp parallel for schedule(static)
        for (size_t x = 0; x < this->w; ++x) {
            for (size_t y = 0; y < this->h; ++y) {

                size_t i = x + y*this->w;

                if (x == 0) {
                    if (y == 0) {
                        // Both on left col and on bottom row.
                        for (size_t n = 0; n < N; n++) {
                            gradf[i][2*n] = (f[i+1][n] - f[i][n]) / this->dx[0]; // x
                            gradf[i][2*n+1] = (f[i+this->w][n] - f[i][n]) / this->dx[1]; // y
                        }
                    } else if (y == this->h-1) {
                        // On left col and on top row
                        for (size_t n = 0; n < N; n++) {
                            gradf[i][2*n] = (f[i+1][n] - f[i][n]) / this->dx[0]; // x
                            gradf[i][2*n+1] = (f[i][n] - f[i-this->w][n]) / this->dx[1]; // y
                        }
                    } else {
                        // Somewhere in middle of left col
                        for (size_t n = 0; n < N; n++) {
                            gradf[i][2*n] = (f[i+1][n] - f[i][n]) / this->dx[0]; // x
                            gradf[i][2*n+1] = (f[i+this->w][n] - f[i-this->w][n]) / (2*this->dx[1]); // y
                        }
                    }

                } else if (x == this->w-1) {
                    if (y == 0) {
                        // On right col and on bottom row
                        for (size_t n = 0; n < N; n++) {
                            gradf[i][2*n] = (f[i][n] - f[i-1][n]) / this->dx[0]; // x
                            gradf[i][2*n+1] = (f[i+this->w][n] - f[i][n]) / this->dx[1]; // y
                        }

                    } else if (y == this->h-1) {
                        // On right col and on top row
                        for (size_t n = 0; n < N; n++) {
                            gradf[i][2*n] = (f[i][n] - f[i-1][n]) / this->dx[0]; // x
                            gradf[i][2*n+1] = (f[i][n] - f[i-this->w][n]) / this->dx[1]; // y
                        }
                    } else {
                        // Somewhere in middle of right col
                        for (size_t n = 0; n < N; n++) {
                            gradf[i][2*n] = (f[i][n] - f[i-1][n]) / this->dx[0]; // x
                            gradf[i][2*n+1] = (f[i+this->w][n] - f[i-this->w][n]) / (2*this->dx[1]); // y
                        }
                    }

                } else {
                    // Not on far left or far right cols
                    if (y == 0) {
                        // Somewhere on bottom row
                        for (size_t n = 0; n < N; n++) {
                            gradf[i][2*n] = (f[i+1][n] - f[i-1][n]) / (2*this->dx[0]); // x
                            gradf[i][2*n+1] = (f[i+this->w][n] - f[i][n]) / this->dx[1]; // y
                        }
                    } else if (y == this->h-1) {
                        // Somewhere on top row
                        for (size_t n = 0; n < N; n++) {
                            gradf[i][2*n] = (f[i+1][n] - f[i-1][n]) / (2*this->dx[0]); // x
                            gradf[i][2*n+1] = (f[i][n] - f[i-this->w][n]) / this->dx[1]; // y
                        }
                    } else {
                        // Somewhere in middle: we can assume that there are elements to l, r, u and d.
                        for (size_t n = 0; n < N; n++) {
                            gradf[i][2*n] = (f[i+1][n] - f[i-1][n]) / (2*this->dx[0]); // x
                            gradf[i][2*n+1] = (f[i+this->w][n] - f[i-this->w][n]) / (2*this->dx[1]); // y
                        }
                    }
                }

                if constexpr (debug_gradf) {
                    std::cout << "("<<x<<","<<y<<") gradf=" << gradf[i] << std::endl;
                }
            }
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

    T linear_expression (const T& x) const { return T{1.31} + T{2.3333} * x; }
    T quadratic_expression (const T& x) const { return T{1.31} + T{2.3333} * x * x; }
    T exponential_expression (const T& x) const { return T{1.05} + T{0.26} * std::exp (T{2.3} * x); }
    T logarithmic_expression (const T& x) const { return T{2.32} + T{1.29} * std::log (T{2.3} * (x+T{0.2})); }

    T expression_function (const T& x, const expression_form& ef) const
    {
        T rtn = T{0};
        switch (ef) {
        case expression_form::lin:
            rtn = this->linear_expression (x);
            break;
        case expression_form::quad:
            rtn = this->quadratic_expression (x);
            break;
        case expression_form::exp:
            rtn = this->exponential_expression (x);
            break;
        case expression_form::log:
            rtn = this->logarithmic_expression (x);
            break;
        default:
            break;
        }
        return rtn;
    }

    T rcpt_expression_function (const T& x, const size_t i) const { return this->expression_function (x, this->rcpt_form[i]); }
    T lgnd_expression_function (const T& x, const size_t i) const { return this->expression_function (x, this->lgnd_form[i]); }

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

    //! Remove half of the tissue
    void ablate (tissue_region region)
    {
        switch (region) {
        case tissue_region::top_half:
            this->ablate_top_half();
            break;
        case tissue_region::bottom_half:
            this->ablate_bottom_half();
            break;
        case tissue_region::left_half:
            this->ablate_left_half();
            break;
        case tissue_region::right_half:
        default:
            this->ablate_right_half();
            break;
        }
    }

    void ablate_right_half()
    {
        size_t newsize = (this->w/2)*this->h;

        morph::vVector<morph::Vector<T,2>> posn_new(newsize);
        morph::vVector<morph::Vector<T,N>> rcpt_new(newsize);
        morph::vVector<morph::Vector<T,2*N>> rcpt_grad_new(newsize);
        morph::vVector<morph::Vector<T,N>> lgnd_new(newsize);
        morph::vVector<morph::Vector<T,2*N>> lgnd_grad_new(newsize);

        // Copy the left hand half to keep
        size_t k = 0;
        for (size_t i = 0; i < this->h; ++i) {
            for (size_t j = 0; j < this->w/2; ++j) {
                posn_new[k] = this->posn[i*this->w+j];
                rcpt_new[k] = this->rcpt[i*this->w+j];
                rcpt_grad_new[k] = this->rcpt_grad[i*this->w+j];
                lgnd_new[k] = this->lgnd[i*this->w+j];
                lgnd_grad_new[k] = this->lgnd_grad[i*this->w+j];
                k++;
            }
        }
        // Now switch the 'news' to olds and update this->w.
        this->posn.swap(posn_new);
        this->rcpt.swap(rcpt_new);
        this->rcpt_grad.swap(rcpt_grad_new);
        this->lgnd.swap(lgnd_new);
        this->lgnd_grad.swap(lgnd_grad_new);
        this->w = this->w/2;
    }

    void ablate_left_half()
    {
        size_t newsize = (this->w/2)*this->h;

        morph::vVector<morph::Vector<T,2>> posn_new(newsize);
        morph::vVector<morph::Vector<T,N>> rcpt_new(newsize);
        morph::vVector<morph::Vector<T,2*N>> rcpt_grad_new(newsize);
        morph::vVector<morph::Vector<T,N>> lgnd_new(newsize);
        morph::vVector<morph::Vector<T,2*N>> lgnd_grad_new(newsize);

        // Copy the right hand half to keep
        size_t k = 0;
        for (size_t i = 0; i < this->h; ++i) {
            for (size_t j = this->w/2; j < this->w; ++j) {
                posn_new[k] = this->posn[i*this->w+j];
                rcpt_new[k] = this->rcpt[i*this->w+j];
                rcpt_grad_new[k] = this->rcpt_grad[i*this->w+j];
                lgnd_new[k] = this->lgnd[i*this->w+j];
                lgnd_grad_new[k] = this->lgnd_grad[i*this->w+j];
                k++;
            }
        }
        // Now switch the 'news' to olds and update this->w.
        this->posn.swap(posn_new);
        this->rcpt.swap(rcpt_new);
        this->rcpt_grad.swap(rcpt_grad_new);
        this->lgnd.swap(lgnd_new);
        this->lgnd_grad.swap(lgnd_grad_new);
        this->w = this->w/2;
    }

    void ablate_top_half()
    {
        size_t newsize = this->w*(this->h/2);

        morph::vVector<morph::Vector<T,2>> posn_new(newsize);
        morph::vVector<morph::Vector<T,N>> rcpt_new(newsize);
        morph::vVector<morph::Vector<T,2*N>> rcpt_grad_new(newsize);
        morph::vVector<morph::Vector<T,N>> lgnd_new(newsize);
        morph::vVector<morph::Vector<T,2*N>> lgnd_grad_new(newsize);

        // Copy the bottom half to keep
        size_t k = 0;
        for (size_t i = 0; i < this->h/2; ++i) {
            for (size_t j = 0; j < this->w; ++j) {
                posn_new[k] = this->posn[i*this->w+j];
                rcpt_new[k] = this->rcpt[i*this->w+j];
                rcpt_grad_new[k] = this->rcpt_grad[i*this->w+j];
                lgnd_new[k] = this->lgnd[i*this->w+j];
                lgnd_grad_new[k] = this->lgnd_grad[i*this->w+j];
                k++;
            }
        }

        // Now switch the 'news' to olds and update this->w.
        this->posn.swap(posn_new);
        this->rcpt.swap(rcpt_new);
        this->rcpt_grad.swap(rcpt_grad_new);
        this->lgnd.swap(lgnd_new);
        this->lgnd_grad.swap(lgnd_grad_new);
        this->h = this->h/2;
    }

    void ablate_bottom_half()
    {
        size_t newsize = this->w*(this->h/2);

        morph::vVector<morph::Vector<T,2>> posn_new(newsize);
        morph::vVector<morph::Vector<T,N>> rcpt_new(newsize);
        morph::vVector<morph::Vector<T,2*N>> rcpt_grad_new(newsize);
        morph::vVector<morph::Vector<T,N>> lgnd_new(newsize);
        morph::vVector<morph::Vector<T,2*N>> lgnd_grad_new(newsize);

        // Copy the top half to keep
        size_t k = 0;
        for (size_t i = this->h/2; i < this->h; ++i) {
            for (size_t j = 0; j < this->w; ++j) {
                posn_new[k] = this->posn[i*this->w+j];
                rcpt_new[k] = this->rcpt[i*this->w+j];
                rcpt_grad_new[k] = this->rcpt_grad[i*this->w+j];
                lgnd_new[k] = this->lgnd[i*this->w+j];
                lgnd_grad_new[k] = this->lgnd_grad[i*this->w+j];
                k++;
            }
        }

        // Now switch the 'news' to olds and update this->w.
        this->posn.swap(posn_new);
        this->rcpt.swap(rcpt_new);
        this->rcpt_grad.swap(rcpt_grad_new);
        this->lgnd.swap(lgnd_new);
        this->lgnd_grad.swap(lgnd_grad_new);
        this->h = this->h/2;
    }

    // Receptor knockout should be easy, right?
    void receptor_knockout (size_t idx)
    {
        if (idx >= N) { throw std::runtime_error ("receptor index out of range"); }
        for (auto& r : this->rcpt) { r[idx] = T{0}; }
        this->compute_gradients();
    }

    void ligand_knockout (size_t idx)
    {
        if (idx >= N) { throw std::runtime_error ("ligand index out of range"); }
        for (auto& l : this->lgnd) { l[idx] = T{0}; }
        this->compute_gradients();
    }

    void receptor_knockdown (size_t idx, T amount)
    {
        if (idx >= N) { throw std::runtime_error ("receptor index out of range"); }
        for (auto& r : this->rcpt) { r[idx] -= amount; }
        this->compute_gradients();
    }

    void ligand_knockdown (size_t idx, T amount)
    {
        if (idx >= N) { throw std::runtime_error ("ligand index out of range"); }
        for (auto& l : this->lgnd) { l[idx] -= amount; }
        this->compute_gradients();
    }

    // For receptor index idx, knock-in affected proportion (range [0,1]) of cells,
    // adding amount as a constant increase in receptor expression
    static constexpr bool random_knockin = false;
    void receptor_knockin (size_t idx, T affected, T amount)
    {
        if (idx >= N) { throw std::runtime_error ("receptor index out of range"); }
        if (affected < T{0} || affected > T{1}) { throw std::runtime_error ("'affected' proportion out of range"); }
        if constexpr (random_knockin == true) {
            morph::RandUniform<T> rng(0,1);
            std::vector<T> rns = rng.get(this->rcpt.size());
            size_t ri = 0;
            for (auto& r : this->rcpt) {
                if (rns[ri++] < affected) { r[idx] += amount; }
            }
        } else {
            size_t ri = 0;
            for (auto& r : this->rcpt) {
                if (ri++%2 == 0) { r[idx] += amount; }
            }
        }
        this->compute_gradients();
    }

    void ligand_knockin (size_t idx, T affected, T amount)
    {
        if (idx >= N) { throw std::runtime_error ("ligand index out of range"); }
        if (affected < T{0} || affected > T{1}) { throw std::runtime_error ("'affected' proportion out of range"); }
        if constexpr (random_knockin == true) {
            morph::RandUniform<T> rng(0,1);
            std::vector<T> rns = rng.get(this->lgnd.size());
            size_t ri = 0;
            for (auto& l : this->lgnd) {
                if (rns[ri++] < affected) { l[idx] += amount; }
            }
        } else {
            size_t ri = 0;
            for (auto& l : this->lgnd) {
                if (ri++%2 == 0) { l[idx] += amount; }
            }
        }
        this->compute_gradients();
    }

    //! Cut a square patch (side legnth sz) of the tissue out and rotate it by a number of quarter_rotations
    void graftrotate (morph::Vector<size_t,2> x1, size_t sz, size_t quarter_rotations)
    {
        if (quarter_rotations%4 == 0) { return; }

        // Potential for runtime crash here, so check x1 and sz.
        if (x1[0] + sz > this->w || x1[1] + sz > this->h) {
            throw std::runtime_error ("Patch for rotation is off the edge of the tissue");
        }

        morph::vVector<morph::Vector<T,N>> rcpt_graft (sz*sz);
        morph::vVector<morph::Vector<T,N>> lgnd_graft (sz*sz);

        for (size_t rot = 1; rot <= quarter_rotations; ++rot) {
            // Copy the square
            for (size_t j = 0; j < sz; ++j) { // y
                for (size_t i = 0; i < sz; ++i) { // x
                    rcpt_graft[i+sz*j] = this->rcpt[(j+x1[1])*this->w+(i+x1[0])];
                    lgnd_graft[i+sz*j] = this->lgnd[(j+x1[1])*this->w+(i+x1[0])];
                }
            }
            // Write back, rotating 90 clockwise once
            for (size_t j = 0; j < sz; ++j) {
                for (size_t i = 0; i < sz; ++i) {
                    size_t ridx = (sz-1-j) + (i*sz);
                    this->rcpt[(j+x1[1])*this->w+(i+x1[0])] = rcpt_graft[ridx];
                    this->lgnd[(j+x1[1])*this->w+(i+x1[0])] = lgnd_graft[ridx];
                }
            }
        }

        this->compute_gradients();
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

    //! Take tissue and mirror the receptor/ligand expressions from left to right (or
    //! from right to left, top to bottom, etc)
    void compound_tissue (tissue_region to_mirror = tissue_region::left_half)
    {
        switch (to_mirror) {
        case tissue_region::right_half:
            // Copy right to left, mirroring.
            for (size_t j = 0; j < this->h; ++j) {
                for (size_t i = 0; i < this->w/2; ++i) {
                    this->rcpt[j*this->w + i] = this->rcpt[j*this->w + (this->w-i-1)];
                    this->lgnd[j*this->w + i] = this->lgnd[j*this->w + (this->w-i-1)];
                }
            }
            break;
        case tissue_region::top_half:
            // Copy top to bottom, mirroring.
            for (size_t j = 0; j < this->h/2; ++j) {
                for (size_t i = 0; i < this->w; ++i) {
                    this->rcpt[this->w*j + i] = this->rcpt[this->w*this->h - (this->w*(j+1)) + i];
                    this->lgnd[this->w*j + i] = this->lgnd[this->w*this->h - (this->w*(j+1)) + i];
                }
            }
            break;
        case tissue_region::bottom_half:
            // Copy bottom to top, mirroring.
            for (size_t j = 0; j < this->h/2; ++j) {
                for (size_t i = 0; i < this->w; ++i) {
                    this->rcpt[this->w*this->h - (this->w*(j+1)) + i] = this->rcpt[this->w*j + i];
                    this->lgnd[this->w*this->h - (this->w*(j+1)) + i] = this->lgnd[this->w*j + i];
                }
            }
            break;
        case tissue_region::left_half:
        default:
            // Copy left to right, mirroring.
            std::cout << "Copy left to right, from x=0 to x=" << this->w/2 << std::endl;
            for (size_t j = 0; j < this->h; ++j) {
                for (size_t i = 0; i < this->w/2; ++i) {
                    std::cout << "i=" << i << ", j=" << j
                              << ". Copying rcpt["<<(j*this->w + i)<<"] to rcpt["
                              <<(j*this->w + (this->w-i-1))<<"]." << std::endl;
                    this->rcpt[j*this->w + (this->w-i-1)] = this->rcpt[j*this->w + i];
                    this->lgnd[j*this->w + (this->w-i-1)] = this->lgnd[j*this->w + i];
                }
            }
            break;
        }

        this->compute_gradients();
    }
};
