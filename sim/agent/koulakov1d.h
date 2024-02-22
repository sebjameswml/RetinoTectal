// Implementation of Koulakov 1D model, in simplest form.
#pragma once

#include <cmath>
#include <memory>
#include <morph/vvec.h>
#include <morph/Random.h>

// Define possible experiments
enum class experiment {
    wildtype,            // Fig 4, Fig 5A
    knockin_hetero,      // Fig 5B
    knockin_homo,        // Fig 5C
    knockin_steeper_ret, // Fig 6A
    knockin_reduced_lig, // Fig 6B
    knockin_steeper_lig, // Fig 6C
};

// Number of elements in arrays is compile-time fixed and global
static constexpr int N = 100;

// Koulakov 1D model class
template <experiment E = experiment::wildtype>
struct k1d
{
    k1d() // Constructor does init
    {
        // Resize arrays
        this->la.resize (N);
        this->ra.resize (N);
        this->rgc_for_sc_idx.resize (N);

        // Make a randomly shuffled initial index array
        this->rgc_for_sc_idx.arange (0, N, 1);
        this->rgc_for_sc_idx.shuffle();

        // Set values for LA and RA and initial projections of rgcs
        for (int k = 0; k < N; ++k) {
            this->la[k] = this->la_k (k);
            this->ra[k] = this->ra_i (k);
        }

        // Set up random number generators
        this->rng_idx = std::make_unique<morph::RandUniform<int>>(0, N-2);
        this->rng_prob = std::make_unique<morph::RandUniform<float>>(0, 1.0f);
    }

    // Ligand expression function
    float la_k (int k)
    {
        if constexpr (E == experiment::knockin_steeper_lig) { // Fig 6C
            return std::exp (-(static_cast<float>(2*k)/N));
        } else if constexpr (E == experiment::knockin_reduced_lig) { // Fig 6B
            //return -0.25f + std::exp (-(static_cast<float>(k)/N));
            return 0.75f * std::exp (-(static_cast<float>(k)/N)); // I think this one?
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
            // 'Value is increased by 25% of the maximum'. Max is 1, so +0.25
            return (i%2==0 ? 0.0f : 0.25f) + std::exp (-(static_cast<float>(N-i)/N));
        } else if constexpr (E == experiment::knockin_steeper_ret) {
            // Like knockin_hetero but with a steeper retinal expression
            return (i%2==0 ? 0.0f : 0.25f) + std::exp (-(static_cast<float>(2*(N-i))/N));
        } else if constexpr (E == experiment::knockin_homo) {
            // 'Value is increased by 50% of the maximum'. Max is 1, so +0.5
            return (i%2==0 ? 0.0f : 0.5f) + std::exp (-(static_cast<float>(N-i)/N));
        } else if constexpr (E == experiment::wildtype) {
            return std::exp (-(static_cast<float>(N-i)/N));
        } else {
            []<bool flag = false>() { static_assert(flag, "Unexpected experiment E"); }();
        }
    }

    // In the initial description of the model Eq1 describes probabilty of exchange. In
    // methods, Eq 4 describes an alternative probability of exchange, bounded between 0
    // and 1.
    static constexpr bool use_eq1_pex_expression = false;

    void step()
    {
        // choose a random i from 0 to N-2
        int i = this->rng_idx->get();
        int j = i+1;

        // i and j here are the index for two selected sc locations (neighbouring)
        int ra_i = this->rgc_for_sc_idx[i];
        int ra_j = this->rgc_for_sc_idx[j];

        // Determine probability of exchange
        float pex = 0.5f;
        if constexpr (use_eq1_pex_expression) {  // Use paper Eq 1
            pex += alpha * (this->ra[ra_i] - this->ra[ra_j]) * (this->la[i] - this->la[j]);
        } else { // Use Eq 4 from paper methods
            pex += 0.5f * std::tanh (2.0f * alpha * (this->ra[ra_i] - this->ra[ra_j]) * (this->la[i] - this->la[j]) );
        }

        float p = this->rng_prob->get();

        if (p < pex) { // then exchange
            this->rgc_for_sc_idx[j] = ra_i;
            this->rgc_for_sc_idx[i] = ra_j;
        }
    }

    // Collicular ligand expressions, by index on SC
    morph::vvec<float> la;
    // Retinal ganglion cell receptor expression by index on retina
    morph::vvec<float> ra;
    // Holds the origin index i of the RGC that terminates at SC location k
    morph::vvec<int> rgc_for_sc_idx;
    // The alpha parameter default value
    float alpha = 30.0f;
    // Random number gen to select indices
    std::unique_ptr<morph::RandUniform<int>> rng_idx;
    // Random number gen to test probabilities
    std::unique_ptr<morph::RandUniform<float>> rng_prob;
};
