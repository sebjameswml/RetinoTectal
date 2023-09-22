#include "branch.h"

// Stochastic branches which estimate ligand gradient based on Mortimer paper.

// \tparam T floating point type for the numbers
// \tparam N number of receptor-ligand pairs in the system
// \tparam R The number of locations at which to sample the ligand gradient.
//           For realism, select R=1001 or so
template<typename T, size_t N, size_t R>
struct branch_stochastic : public branch<T,N>
{
    // Params for the receptor binding model
    static constexpr T gc_width = T{1}; // growth cone width. 10 um is 1e-5 and is a realistic value

    // Integrate estimates of the ligand gradient in which this branch is moving
    static constexpr size_t gh_size = 30;
    morph::vec<morph::vec<T, 2*N>, gh_size> gradient_history;
    size_t gh_i = 0;

    // Estimate the ligand gradient, given the true ligand gradient and receptor noise
    // \param mu true ligand gradient
    // \param gamma true ligand expression
    // \return Estimate of ligand gradient
    virtual morph::vec<T, 2*N> estimate_ligand_gradient (morph::vec<T,2*N>& mu,
                                                         morph::vec<T,N>& gamma)
    {
        // Compile-time assertion that R is odd
        static_assert (R%2 == 1);

        // Choose a set of R locations across the growth cone.
        morph::vec<T, R> r; // The positions, r across the growth cone. Must be -ve to +ve
        morph::vec<T, R> c; // The concentration of ligand, c
        morph::vec<T, R> p; // The probability of binding, p
        morph::vec<bool, R> b; // The state of binding
        morph::vec<T, 2*N> lg; // Current ligand gradient estimate container

        //morph::RandUniform<T, std::mt19937> rng;
        morph::vec<T, R> rns;

        static constexpr bool debug_stochastic = false;
        static constexpr bool debug_stochastic2 = false;

        // Used to compute prob. of binding
        T kp = T{1};
        // Artificially amplify the concentration gradient
        T amplify = T{100};

        for (size_t i = 0; i < N; ++i) { // For each receptor-ligand pairing
            for (size_t d = 0; d < 2; ++d) { // For each dimension

                morph::vvec<T> bound; // posns of bound receptors
                brng::i()->get(rns);

                for (size_t ir = 0; ir < R; ++ir) { // For each of R evenly spaced samples

                    // Long winded for now, until it's sorted
                    r[ir] = gc_width / T{R} * (T{0.5}+ir) - gc_width / T{R} * (T{R}/T{2}) ; // position
                    c[ir] = gamma[i] * (1 + mu[2*i+d] * r[ir] * amplify); // concentration
                    p[ir] = c[ir] / (kp + c[ir]);               // probability of binding
                    b[ir] = (rns[ir] < p[ir]) ? true : false;   // sampled state of binding

                    if (b[ir] == true) { bound.push_back (r[ir]); }
                }

                if constexpr (debug_stochastic) {
                    std::cout << "r=" << r << std::endl;
                    std::cout << "c=" << c << std::endl;
                    std::cout << "p=" << p << std::endl;
                    std::cout << "b=" << b << std::endl;
                    std::cout << "rns=" << rns << std::endl;
                    std::cout << "bound=" << bound << std::endl;
                }
                // Linear fit to the bound receptors, weighted by position
                morph::vvec<T> bound_w = bound.abs();
                if (bound.size() < 2) {
                    lg[2*i+d] = T{0};
                } else {
                    morph::vec<T,2> mc = morph::MathAlgo::linregr (bound, bound_w);
                    lg[2*i+d] = mc[0];
                }
                if constexpr (debug_stochastic2) {
                    std::cout << "Gradient " << mu[2*i+d] << " estimated as: " << lg[2*i+d] << "; bound size: " << bound.size() << std::endl;
                    std::cout << "c[0] - c[end] = " << (c.front() - c.back()) << std::endl;
                }
            }
        }

        //std::cout << "Ligand gradient is " << mu << "\nEstimate gradient is " << lg << std::endl;
        this->gradient_history[gh_i] = lg;
        //std::cout << "\n\n\ngh_i = " << gh_i << std::endl << std::endl << std::endl;
        gh_i = (gh_i+1)%gh_size;

        return gradient_history.mean();
    }

    // This function is duplicated (cf struct branch) so I can pass in a vector<branch_stochastic<T, N, R>&
    void compute_next (const std::vector<branch_stochastic<T, N, R>>& branches,
                       const guidingtissue<T, N>* source_tissue,
                       const guidingtissue<T, N>* tissue,
                       const morph::vec<T, 5>& m,
                       const morph::vec<T, 2*N>& rns)
    {
        // Current location is named b
        morph::vec<T, 2> b = this->current;
        // Chemoaffinity, graded by origin position (i.e termination zone) of each retinal axon
        morph::vec<T, 2> G = this->compute_chemo (source_tissue, tissue);

        // Competition, C, and Axon-axon interactions, I&J, computed during the same loop
        // over the other branches
        morph::vec<T, 2> C = {0, 0};
        morph::vec<T, 2> I = {0, 0};
        morph::vec<T, 2> J = {0, 0};

        // Other branches are called k, making a set (3 sets) B_b, with a number of members that I call n_k etc
        T n_k = T{0};
        T n_ki = T{0};
        T n_kj = T{0};
        for (auto k : branches) {
            if (k.id == this->id) { continue; } // Don't interact with self
            std::bitset<3> cij_added = this->compute_for_branch (source_tissue, (&k), C, I, J);
            n_k += cij_added[0] ? T{1} : T{0};
            n_ki += cij_added[1] ? T{1} : T{0};
            n_kj += cij_added[2] ? T{1} : T{0};
        }

        // Do the 1/|B_b| multiplication to normalize C and I(!!)
        if (n_k > T{0}) { C = C/n_k; } // else C will be {0,0} still
        I = n_ki > T{0} ? I/n_ki : I;
        J = n_kj > T{0} ? J/n_kj : J;

        if constexpr (this->store_interaction_history == true) {
            // Add I to history & rotate, reducing effect by 90% due to one time step passing
            for (size_t i = 0; i>this->ihs-2; ++i) {
                this->Ihist[i] = this->Ihist[i+1] * T{0.98};
            }
            this->Ihist[this->ihs-1] = I;
            // Reset I and then sum Ihist to get effective I
            I = {0,0};
            for (size_t i = 0; i<this->ihs; ++i) { I += this->Ihist[i]; }
        }
        // Collected non-border movement components
        morph::vec<T, 2> nonB = G * m[0] + J * m[1] + I * m[2]  + C * m[3];

        // Border effect. A 'force' to move agents back inside the tissue boundary
        morph::vec<T, 2> B = this->apply_border_effect (tissue, nonB);

        // The change in b from the model and border effects:
        morph::vec<T, 2> db = (nonB + B * m[4]);

        // Finally add b and db to get next
        this->next = b + db;
    }
};
