// Implementation of Koulakov/Tsigankov model as used in Trippett 2011
#pragma once

#include <cmath>
#include <memory>
#include <list>
#include <morph/vec.h>
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
template <experiment E = experiment::wildtype, int N = 100, typename F = float>
struct ktt1d
{
    ktt1d()
    {
        // Resize arrays
        this->col_la.resize (N);
        this->col_ra.resize (N);
        this->ret_la.resize (N);
        this->ret_ra.resize (N);
        this->ret_synapses.resize (N);
        this->ret_synapse_density.resize (N*N, 0.0f);
        this->den_synapse_counts.resize (N, 0);

        // Set values for LA and RA and initial projections of rgcs
        for (int k = 0; k < N; ++k) {
            this->col_la[k] = this->collicular_ligand (k);
            this->col_ra[k] = this->collicular_receptor (k);
            this->ret_la[k] = this->retinal_ligand (k);
            this->ret_ra[k] = this->retinal_receptor (k);
        }
        // Compute gradients
        this->grad_col_la = col_la.diff() / this->d;
        this->grad_ret_ra = ret_ra.diff() / this->d;

        // Set up random number generators
        this->rng_idx = std::make_unique<morph::RandUniform<int>>(0, N-2);
        this->rng_synidx = std::make_unique<morph::RandUniform<int>>(0, synapse_upper_limit);
        this->rng_prob = std::make_unique<morph::RandUniform<F>>(0, F{1});
    }

    // Ligand expression function
    F collicular_ligand (int i)
    {
        if constexpr (E == experiment::knockin_steeper_lig) {
            return std::exp (-(static_cast<F>(2*i)/N));
        } else if constexpr (E == experiment::knockin_reduced_lig) {
            return 0.75f * std::exp (-(static_cast<F>(i)/N));
        } else {
            return std::exp (-(static_cast<F>(i)/N));
        }
    }

    // Collicular receptor expression
    F collicular_receptor (int i)
    {
        return 0.0f;
    }

    // Retinal ligand expression
    F retinal_ligand (int i)
    {
        return 0.0f;
    }

    // RGC receptor expression function
    F retinal_receptor (int i)
    {
        if constexpr (E == experiment::knockin_hetero
                      || E == experiment::knockin_steeper_lig
                      || E == experiment::knockin_reduced_lig) {
            return (i%2==0 ? 0.0f : 0.25f) + std::exp (-(static_cast<F>(N-i)/N));
        } else if constexpr (E == experiment::knockin_steeper_ret) {
            return (i%2==0 ? 0.0f : 0.25f) + std::exp (-(static_cast<F>(2*(N-i))/N));
        } else if constexpr (E == experiment::knockin_homo) {
            return (i%2==0 ? 0.0f : 0.5f) + std::exp (-(static_cast<F>(N-i)/N));
        } else if constexpr (E == experiment::wildtype) {
            return std::exp (-(static_cast<F>(N-i)/N));
        } else {
            []<bool flag = false>() { static_assert(flag, "Unexpected experiment E"); }();
        }
    }

    // On each step, make an attempt to create a synapse between two randomly chosen end
    // points and make one attempt to destroy an existing synapse.
    void step()
    {
        this->attempt_creation();
        this->attempt_destruction();
    }

    // Experience will suggest a limit here, it's for memory reservation
    static constexpr int synapse_upper_limit = 1000;

    // Set true for debug output to stdout
    static constexpr bool debug_synapse_changes = false;

    // Competition component parameters
    static constexpr F comp_param_A = F{500};      // The competition parameter 'A'. 500 in paper.
    static constexpr F comp_param_B = F{1};
    static constexpr F comp_param_D = F{1};

    // Chemotaxis component parameters
    static constexpr F chem_param_alpha = F{120};  // 120 in paper

    // Activity component parameters
    static constexpr F act_param_a = F{3};         // 3 in paper
    static constexpr F act_param_b = F{11};        // 11 in paper
    static constexpr F act_param_gamma = F{0.05};  // 0.05 in paper

    // Derived constants:
    static constexpr F act_param_minus1_over_2a_squared = F{-1} / (F{2} * act_param_a * act_param_a);
    static constexpr F act_param_minus1_over_b = F{-1} / act_param_b;
    static constexpr F act_param_minusgamma_over_2 = -act_param_gamma / F{2};

    // Compute the three components of delta A for a connection from retinal index ret_i to SC index
    // sc_j where there are currently axonsyns_for_i retinal synapses from i and densyns_for_j
    // synapses connecting through to SC target location j (dendrites of j). synchange controls how
    // many synapses we are adding/removing for connection i->j.
    morph::vec<F, 3> compute_deltaA (const int ret_i, const int sc_j,
                                     const int axonsyns_for_i, const int densyns_for_j,
                                     const int synchange = 1)
    {
        /*
         * 1. change in energy due to chemotactic changes
         */

        //F deltaA_chem = synchange * this->alpha * this->grad_ret_ra[ret_i] * this->grad_col_la[sc_j];
        F deltaA_chem = -synchange * this->alpha * this->ret_ra[ret_i] * this->col_la[sc_j];

        /*
         * 2. change in energy due to activity changes
         */

        // To compute activity change, we need to compute the new synapse (call it s2) against all
        // the existing synapses (s1).
        F sum = F{0};
        for (int ri = 0; ri < N; ++ri) {
            auto si = ret_synapses[ri].begin();
            F ret_i_to_ret_ri_spacing = std::abs (ret_i-ri) * this->d;
            while (si != ret_synapses[ri].end()) {
                // *si is an index on the SC
                F syn_spacing = std::abs (*si - sc_j) * this->d_sc;
                F Ca1a2 = std::exp (act_param_minus1_over_b * ret_i_to_ret_ri_spacing);
                F U = std::exp (act_param_minus1_over_2a_squared * syn_spacing);
                sum += (Ca1a2 * U);
                ++si;
            }
        }
        F deltaA_act = synchange * act_param_minusgamma_over_2 * sum;

        /*
         * 3. change in energy due to competition changes
         */

        int new_densyns_for_j = densyns_for_j + synchange;
        int new_axonsyns_for_i = axonsyns_for_i + synchange;

        F A_comp_before = (-comp_param_A * std::sqrt(static_cast<F>(axonsyns_for_i))
                           + comp_param_B * static_cast<F>(axonsyns_for_i * axonsyns_for_i)
                           + comp_param_D * static_cast<F>(densyns_for_j * densyns_for_j));

        F A_comp_after = (-comp_param_A * std::sqrt(static_cast<F>(new_axonsyns_for_i))
                           + comp_param_B * static_cast<F>(new_axonsyns_for_i * new_axonsyns_for_i)
                           + comp_param_D * static_cast<F>(new_densyns_for_j * new_densyns_for_j));

        F deltaA_comp = A_comp_after - A_comp_before;

        // Pack and return
        return morph::vec<F, 3>({deltaA_chem, deltaA_act, deltaA_comp});
    }

    // Attempt to create a synapse
    void attempt_creation()
    {
        // choose a random i from 0 to N-2. i indexes retina, j indexes SC.
        int i = this->rng_idx->get();
        int j = this->rng_idx->get();

        // How many synapses for axon i at present?
        int axonsyns_for_i = static_cast<int>(this->ret_synapses[i].size());

        // And how many synapses for dendrite j?
        // Go through all ret_synapses[i] and count how many instances of j there are!
        int densyns_for_j = den_synapse_counts[j]; // ?? how to get. Need to record and store.

        morph::vec<F, 3> deltaA_cmpts = compute_deltaA (i, j, axonsyns_for_i, densyns_for_j, +1);

        // probability of accepting a change (addition of synapse)
        F p_accept = F{1} / (F{1} + std::exp (F{4} * deltaA_cmpts.sum()));

        F p = this->rng_prob->get();

        if (this->ret_synapses[i].size() > synapse_upper_limit) {
            throw std::runtime_error ("Hit synapse_upper_limit");
        }

        if (p < p_accept) {
            // then add synapse
            if constexpr (debug_synapse_changes == true) {
                std::cout << "Add synapse for SC dendrite " << j << " to ret axon " << i
                          << " where chem + act + comp = "
                          << deltaA_cmpts[0] << " + "
                          << deltaA_cmpts[1] << " + "
                          << deltaA_cmpts[2] << " = "  << deltaA_cmpts.sum() << " => p_acc = " << p_accept <<   std::endl;
            }
            this->ret_synapses[i].push_back (j);
            this->den_synapse_counts[j]++;
            ++this->n_syn;

        } else {
            if constexpr (debug_synapse_changes == true) {
                std::cout << "   NO ADD   for SC dendrite " << j << " to ret axon " << i
                          << " where chem + act + comp = "
                          << deltaA_cmpts[0] << " + "
                          << deltaA_cmpts[1] << " + "
                          << deltaA_cmpts[2] << " = "  << deltaA_cmpts.sum() << " => p_acc = " << p_accept <<   std::endl;
            }
        }
    }

    // Attempt to remove a synapse
    void attempt_destruction()
    {
        // choose a random i from 0 to N-2. i indexes retina, j indexes SC.
        int i = this->rng_idx->get();
        // j is different from attempt_creation. Here, we choose a random one of the
        // elements of the vvec ret_synapses[i]. The modulus ensures we choose an index
        // in range.
        int axonsyns_for_i = static_cast<int>(this->ret_synapses[i].size());
        // Can't destruct if no synapses
        if (axonsyns_for_i == 0) { return; }
        int j_idx = this->rng_synidx->get() % axonsyns_for_i;

        auto j_iter = this->ret_synapses[i].begin();
        for (int jj = 0; jj < j_idx && j_iter != this->ret_synapses[i].end(); ++jj) { ++j_iter; }

        int j = *j_iter;

        // How many synapses for dendrite j_iter?
        int densyns_for_j = den_synapse_counts[j];

        morph::vec<F, 3> deltaA_cmpts = compute_deltaA (i, j, axonsyns_for_i, densyns_for_j, -1);

        // probability of accepting a change (removal of synapse)
        F p_accept = F{1} / (F{1} + std::exp (F{4} * deltaA_cmpts.sum()));

        F p = this->rng_prob->get();
        if (p < p_accept) { // then accept change
            // Remove synapse
            //std::cout << "Erase synapse for SC dendrite " << j << " from ret axon " << i << std::endl;
            this->ret_synapses[i].erase (j_iter);
            this->den_synapse_counts[j]--;
            --this->n_syn;
        //} else {
            //std::cout << "No erase.\n";
        }
    }

    void compute_ret_synapse_density()
    {
        this->ret_synapse_density.zero();
        for (int i = 0; i < N; ++i) {
            for (auto syn : ret_synapses[i]) {
                // With this index, when viewed on a Grid with GridOrder::bottomleft_to_topright
                // (row major), the retinal address will be given by the vertical (y) axis and the
                // Synapse location by the horizontal (x) axis.
                int idx = i * N + syn;
                this->ret_synapse_density[idx] += 1.0f;
            }
        }
        this->ret_synapse_density /= this->ret_synapse_density.max();
    }

    // SC
    // Collicular ligand expression, by index on SC. ephrin A only as 1D. May become
    // available ligand expression after subtraction of ret_la.
    morph::vvec<F> col_la;
    // Gradient of ligand expression
    morph::vvec<F> grad_col_la;
    // Collicular receptor expression
    morph::vvec<F> col_ra;

    // RETINA
    // Retinal ganglion cell receptor expression by index on retina
    morph::vvec<F> ret_ra;
    // Gradient of retinal receptor expression
    morph::vvec<F> grad_ret_ra;
    // Retinal ligand expression
    morph::vvec<F> ret_la;
    // synapse[i] contains a vector of the SC dendrite indices to which axon i makes
    // synapses. The actual branches of the dendrites (and axons) are ignored. Has a
    // ret_ prefix, because the list of synapses is indexed by the retinal origin.
    //morph::vvec<morph::vvec<int>> ret_synapses;
    morph::vvec<std::list<int>> ret_synapses;
    morph::vvec<float> ret_synapse_density; // float, not F because? Because of Grids?

    // Have to keep a record of how many synapses there are for each dendrite
    morph::vvec<int> den_synapse_counts;

    // Track total number of synapses
    unsigned long long int n_syn = 0;

    // The current activation of the system, if it's required to keep it around
    F A_total = F{0};
    F A_chem = F{0};
    F A_comp = F{0};
    F A_act = F{0};
    // The alpha parameter default value
    F alpha = chem_param_alpha;
    // The distance d between adjacent cells on the retina
    F d = F{1};
    // The SC distances are the same
    F d_sc = F{1};
    // Random number gen to select indices
    std::unique_ptr<morph::RandUniform<int>> rng_idx;
    // To select indices from synapse lists
    std::unique_ptr<morph::RandUniform<int>> rng_synidx;
    // Random number gen to test probabilities
    std::unique_ptr<morph::RandUniform<F>> rng_prob;
};
