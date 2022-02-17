/*
 * Base class for branches
 */

#pragma once

#include <morph/Vector.h>

// A retinotectal axon branch class. Holds current and historical positions, a preferred
// termination zone, and the algorithm for computing the next position. Could derive
// from an 'agent' base class. The branches express N Ephrin receptor types (EphA, EphB etc)
template<typename T, size_t N>
struct branch_base
{
    virtual void init() = 0;

    // Place the next computed location for path in 'next' so that while computing, we
    // don't modify the numbers we're working from. After looping through all branches,
    // place current into next and move on to next time step.
    morph::Vector<T, 2> current;
    morph::Vector<T, 2> next;

    // The 'target' location for the axon/branch. This is the origin location in the source tissue (retina)
    morph::Vector<T, 2> target;

    // Interaction parameters for this branch, taken from the soma in the source
    // tissue. This is the N receptor expressions at the growth cone.
    morph::Vector<T, N> rcpt;
    // Ligands on the RGCs may interact with receptors in the tectal tissue
    morph::Vector<T, N> lgnd;

    // By default, axons reckoned to be outside tissue. Can also be used to mean "active".
    bool entered = false;

    // A sequence id
    int id = 0;
    // An id for the parent axon (there are many branches per axon)
    int aid = 0; // this is id/bpa (computed with integer division)
    // Noise gain if required
    T noise_gain = T{1};
};
