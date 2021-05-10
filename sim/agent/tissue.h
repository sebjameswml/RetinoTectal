#pragma once

#include <morph/Vector.h>

//! Some 2D brain tissue, arranged as a cartesian grid of neurons of width w and height h.
template<typename T>
struct tissue
{
    size_t w;
    size_t h;
    size_t num() { return w * h; }
    // The distance from neuron to neuron; horizontally and vertically.
    morph::Vector<T,2> dx;
    // Position of bottom left neuron
    morph::Vector<T,2> x0;
    // Positions stored in raster fashion from bottom left
    morph::vVector<morph::Vector<T,2>> posn;

    tissue(size_t _w, size_t _h, morph::Vector<T,2> _dx, morph::Vector<T,2> _x0)
        : w(_w), h(_h), dx(_dx), x0(_x0)
    {
        this->posn.resize(w*h);
        for (size_t i = 0; i < this->h; ++i) {
            for (size_t j = 0; j < this->w; ++j) {
                posn[i*w+j] = {dx*j, dy*i};
            }
        }
    }
};

//! A retina, for visualisation and for simulated experimental manipulations. This class holds data
template<typename T>
struct retina : public tissue<T>
{
    // A left-right interaction parameter and an up-down interaction parameter for each
    // retinal neuron (hence Vectors with 2 elements, each of which can be +/-)
    morph::vVector<morph::Vector<T,2>> interaction;

    // Init the wildtype retina's interaction parameters, which stand in for the Ephrins.
    retina(size_t _w, size_t _h, morph::Vector<T,2> _dx, morph::Vector<T,2> _x0)
        : tissue (_w, _h, _dx, _x0)
    {
        this->interaction.resize(w*h);
        for (size_t i = 0; i < this->h; ++i) {
            for (size_t j = 0; j < this->w; ++j) {
                //interaction[i*w+j] = {dx[0]*j-x0[0], dx[1]*i-x0[1]};
                this->interaction[i*w+j] = (dx*j)-x0;

            }
        }
    }
};

//! A tectum, for we will manipulate the tectum, too
template<typename T>
struct tectum : public tissue<T>
{
};
