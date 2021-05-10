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

//! Tissue which has guidance parameters in the form of an interaction parameter vector
template<typename T>
struct guidingtissue : public tissue<T>
{
    //! A left-right interaction parameter and an up-down interaction parameter for each
    //! retinal neuron (hence Vectors with 2 elements, each of which can be +/-)
    morph::vVector<morph::Vector<T,2>> interaction;

    //! Init the wildtype retina's interaction parameters, which stand in for the Ephrins.
    guidingtissue(size_t _w, size_t _h, morph::Vector<T,2> _dx, morph::Vector<T,2> _x0)
        : tissue<T> (_w, _h, _dx, _x0)
    {
        this->interaction = this->posn - this->x0;
    }

    //! Move a rectangular patch of tissue at indexed location original, of size
    //! patchsize, and swap it with another region of the retina at newlocn.
    void patchswap (morph::Vector<size_t,2> original,
                    morph::Vector<size_t,2> patchsize,
                    morph::Vector<size_t,2> newlocn)
    {
        // This is all about moving a section of interaction. Have to move row by row though.
        morph::vVector<morph::Vector<T,2>> tmp; // Temporarily holds a section of interaction

        size_t patchrowlen = patchsize[0];
        size_t patchnumrows = patchsize[1];

        tmp.resize (patchrowlen);
        std::cout << "tmp size: " << tmp.size() << ", from " << patchsize << std::endl;

        size_t orig_idx = original[0] + this->w * original[1];
        size_t new_idx = newlocn[0] + this->w * newlocn[1];
        for (size_t i = 0; i < patchnumrows; ++i) {

            if (orig_idx + patchrowlen > this->interaction.size()
                || new_idx + patchrowlen > this->interaction.size()) {
                throw std::runtime_error ("Can't make that patchswap");
            }
            auto oi = this->interaction.begin();
            oi += orig_idx;
            auto ni = this->interaction.begin();
            ni += new_idx;
            std::copy_n (ni, patchrowlen, tmp.begin());
            std::copy_n (oi, patchrowlen, ni);
            std::copy_n (tmp.begin(), patchrowlen, oi);

            orig_idx += this->w;
            new_idx += this->w;
        }
    }
};

#if 0 // In case I need to differentiate retina and tectum. Eph receptors?
//! A retina, for visualisation and for simulated experimental manipulations. This class holds data
template<typename T>
struct retina : public guidingtissue<T>
{
};

//! A tectum, for we will manipulate the tectum, too
template<typename T>
struct tectum : public guidingtissue<T>
{
};
#endif
