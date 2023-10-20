/*
 * Base class for branches
 */

#pragma once

#include <morph/vec.h>

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
    morph::vec<T, 2> current;
    morph::vec<T, 2> next;

    // The 'target' location for the axon/branch. This is the origin location in the source tissue (retina)
    morph::vec<T, 2> target;

    // Interaction parameters for this branch, taken from the soma in the source
    // tissue. This is the N receptor expressions at the growth cone.
    morph::vec<T, N> rcpt;
    // Ligands on the RGCs may interact with receptors in the tectal tissue
    morph::vec<T, N> lgnd;

    // The base level of EphA4 in the retina, copied from the source tissue. Not dynamic.
    T rcpt0_EphA4_base = T{0}; // Base level, before any cis interactions
    T rcpt0_EphA4_free = T{0}; // Free, unbound EphA4 receptors
    T rcpt0_EphA4_cis = T{0};  // phosphorylised (cis-bound)
    T rcpt0_EphA4_cis_min = T{0}; // min phosphorylised amount

    // By default, axons reckoned to be outside tissue. Can also be used to mean "active".
    bool entered = false;

    // A sequence id
    int id = 0;
    // An id for the parent axon (there are many branches per axon)
    int aid = 0; // this is id/bpa (computed with integer division)
    // Noise gain if required
    T noise_gain = T{1};

    // Parameters for special EphA4 model
    T Ax_thresh = T{2.5};
    T A4_thresh = T{1.5};

    T getr() { return this->r; }
    void setr (T _r) { this->r = _r; this->two_r = _r*T{2}; }

protected:
    T r = T{0.04};   // A radius for a growth cone (may be used for visualisation)
    T two_r = T{2}*r;
};
