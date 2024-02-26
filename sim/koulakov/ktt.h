// Implementation of Koulakov/Tsigankov model as used in Trippett 2011
#pragma once

#include <cmath>
#include <memory>
#include <morph/vvec.h>
#include <morph/Random.h>

// Define possible experiments
enum class experiment {
    wildtype,
    knockin_hetero,
    knockin_homo,
    knockin_steeper_ret,
    knockin_reduced_lig,
    knockin_steeper_lig,
};

// Koulakov/Tsigankov competition model as found in Trippett et al 2011
template <experiment E = experiment::wildtype, size_t N = 100>
struct ktt1d
{
    // Experience will suggest a limit here, it's for memory reservation
    static constexpr int synapse_upper_limit = 100;

    ktt1d() // Constructor does init
    {
        // Resize arrays
        this->la.resize (N);
        this->ra.resize (N);
        this->synapses.resize (N);

        this->rgc_for_sc_idx.resize (N);

        // Make a randomly shuffled initial index array
        this->rgc_for_sc_idx.arange (0, N, 1);
        this->rgc_for_sc_idx.shuffle();

        // Set values for LA and RA and initial projections of rgcs
        for (int k = 0; k < N; ++k) {
            this->la[k] = this->la_k (k);
            this->ra[k] = this->ra_i (k);
            this-synapses[k].reserve (synapse_upper_limit);
        }
        // Compute gradients
        this->grad_la = la.diff() / this->d;
        this->grad_ra = ra.diff() / this->d;

        // Set up random number generators
        this->rng_idx = std::make_unique<morph::RandUniform<int>>(0, N-2);
        this->rng_prob = std::make_unique<morph::RandUniform<float>>(0, 1.0f);
    }

    // Ligand expression function
    float la_k (int k)
    {
        if constexpr (E == experiment::knockin_steeper_lig) {
            return std::exp (-(static_cast<float>(2*k)/N));
        } else if constexpr (E == experiment::knockin_reduced_lig) {
            return 0.75f * std::exp (-(static_cast<float>(k)/N));
        } else {
            return std::exp (-(static_cast<float>(k)/N));
        }
    }

    // RGC receptor expression function
    float ra_i (int i)
    {
        if constexpr (E == experiment::knockin_hetero
                      || E == experiment::knockin_steeper_lig
                      || E == experiment::knockin_reduced_lig) {
            return (i%2==0 ? 0.0f : 0.25f) + std::exp (-(static_cast<float>(N-i)/N));
        } else if constexpr (E == experiment::knockin_steeper_ret) {
            return (i%2==0 ? 0.0f : 0.25f) + std::exp (-(static_cast<float>(2*(N-i))/N));
        } else if constexpr (E == experiment::knockin_homo) {
            return (i%2==0 ? 0.0f : 0.5f) + std::exp (-(static_cast<float>(N-i)/N));
        } else if constexpr (E == experiment::wildtype) {
            return std::exp (-(static_cast<float>(N-i)/N));
        } else {
            []<bool flag = false>() { static_assert(flag, "Unexpected experiment E"); }();
        }
    }

    void step()
    {
        // choose a random i from 0 to N-2
        int i = this->rng_idx->get();
        int j = i+1;

        // i and j here are the index for two selected sc locations (neighbouring)
        int ra_i = this->rgc_for_sc_idx[i];
        int ra_j = this->rgc_for_sc_idx[j];

        // Compute probability of exchange here
        //float deltaE = std::exp (-0.5f * alpha * grad_ra[ra_i] * grad_la[i] * this->dsq);

        float deltaA_chem = 0.0f;
        float deltaA_act = 0.0f;
        float deltaA_comp = 0.0f;

        float deltaA = deltaA_chem + deltaA_act + deltaA_comp;

        // probability of accepting a change (addition or removal of synapse)
        float p_accept = 1.0f / (1.0f + std::exp (4.0f * deltaA));

        float p = this->rng_prob->get();

        if (p < p_accept) { // then accept change
            this->rgc_for_sc_idx[j] = ra_i;
            this->rgc_for_sc_idx[i] = ra_j;
        }
    }

    // Collicular ligand expression, by index on SC. ephrin A only as 1D. May become
    // available ligand expression after subtraction of la_ret.
    morph::vvec<float> la;
    morph::vvec<float> grad_la;
    // Collicular receptor expression
    morph::vvec<float> ra_col;
    // Retinal ganglion cell receptor expression by index on retina
    morph::vvec<float> ra;
    morph::vvec<float> grad_ra;
    // Retinal ligand expression
    morph::vvec<float> la_ret;

    // synapse[i] contains a vector of the SC dendrite indices to which axon i makes synapses.
    morph::vvec<morph::vvec<int>> synapses_form1;
    // In this form, synapse[i] contains the source axon and the destination
    // dendrite. Better for randomly selecting synapses.
    morph::vvec<morph::vec<int, 2>> synapses_form2;

    // The alpha parameter default value
    float alpha = 120.0f;
    // The distance d between adjacent cells
    float d = 0.01f;
    // d squared
    float dsq = 0.0001f;
    // Random number gen to select indices
    std::unique_ptr<morph::RandUniform<int>> rng_idx;
    // Random number gen to test probabilities
    std::unique_ptr<morph::RandUniform<float>> rng_prob;
};
