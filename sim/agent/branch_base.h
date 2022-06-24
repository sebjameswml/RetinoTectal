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

    // The base level of EphA4 in the retina, copied from the source tissue. Not dynamic.
    T rcpt0_EphA4;

    // By default, axons reckoned to be outside tissue. Can also be used to mean "active".
    bool entered = false;

    // A sequence id
    int id = 0;
    // An id for the parent axon (there are many branches per axon)
    int aid = 0; // this is id/bpa (computed with integer division)
    // Noise gain if required
    T noise_gain = T{1};

    // Probability of a cluster having EphA4 'side-attachments', hypothesised to be a 'normal' EphA cluster
    T side_attach_prob = T{1};
    T normal_cluster_gain = T{1};
    T enhanced_cluster_gain = T{1}; // An enhanced cluster is one with no EphA4 sides, and could be a larger cluster
    T epha4_attachment_proportion = T{0};
};
