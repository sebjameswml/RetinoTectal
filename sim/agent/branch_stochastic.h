#include "branch.h"

// \tparam T floating point type for the numbers
// \tparam N number of receptor-ligand pairs in the system
// \tparam R The number of locations at which to sample the ligand gradient.
//           For realism, select R=1000 or so
template<typename T, size_t N, size_t R>
struct branch_stochastic : public branch<T,N>
{
    // Params for the receptor binding model
    static constexpr T gc_width = T{1}; // growth cone width. 10 um is 1e-5 and is a realistic value

    // Estimate the ligand gradient, given the true ligand gradient and receptor noise
    virtual morph::Vector<T, 2*N> estimate_ligand_gradient (morph::Vector<T,2*N>& mu,
                                                            morph::Vector<T,N>& gamma)
    {
        // Choose a set of R locations across the growth cone.
        morph::Vector<T, R> r; // The positions, r across the growth cone. Must be -ve to +ve
        morph::Vector<T, R> c; // The concentration of ligand, c
        morph::Vector<T, R> p; // The probability of binding, p
        morph::Vector<bool, R> b; // The state of binding
        morph::Vector<T, 2*N> lg; // rtn container

        morph::RandUniform<T, std::mt19937> rng;
        morph::Vector<T, R> rns;

        static constexpr bool debug_stochastic = false;

        for (size_t i = 0; i < N; ++i) { // For each receptor-ligand pairing
            for (size_t d = 0; d < 2; ++d) { // For each dimension

                morph::vVector<T> bound; // posns of bound receptors
                rng.get(rns);

                for (size_t ir = 0; ir < R; ++ir) { // For each of R evenly spaced samples

                    // Long winded for now, until it's sorted
                    r[ir] = gc_width / T{R} * ir - gc_width / T{R} * (R/2) ; // position
                    c[ir] = gamma[i] * (1 + mu[2*i+d] * r[ir]); // concentration
                    p[ir] = c[ir] / (T{1} + c[ir]); // probability of binding
                    b[ir] = (rns[ir] < p[ir]) ? true : false; // sampled state of binding

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
                morph::vVector<T> bound_w = bound.abs();
                if (bound.size() < 2) {
                    lg[2*i+d] = T{0};
                } else {
                    std::pair<T,T> mc = morph::MathAlgo::linregr (bound, bound_w);
                    lg[2*i+d] = mc.first;
                }
                if constexpr (debug_stochastic) {
                    std::cout << "Gradient " << mu[2*i+d] << " estimated as: " << lg[2*i+d] << "; bound size: " << bound.size() << std::endl;
                }
            }
        }

        //std::cout << "Ligand gradient is " << mu << "\nEstimate gradient is " << lg << std::endl;

        return lg;
    }

    // This function is duplicated (cf struct branch) so I can pass in a vector<branch_stochastic<T, N, R>&
    void compute_next (const std::vector<branch_stochastic<T, N, R>>& branches,
                       const guidingtissue<T, N>* source_tissue,
                       const guidingtissue<T, N>* tissue,
                       const morph::Vector<T, 5>& m,
                       const morph::Vector<T, 2*N>& rns)
    {
        // Current location is named b
        morph::Vector<T, 2> b = this->current;
        // Chemoaffinity, graded by origin position (i.e termination zone) of each retinal axon
        morph::Vector<T, 2> G = this->compute_chemo (source_tissue, tissue);

        // Competition, C, and Axon-axon interactions, I, computed during the same loop
        // over the other branches
        morph::Vector<T, 2> C = {0, 0};
        morph::Vector<T, 2> I = {0, 0};
        morph::Vector<T, 2> J = {0, 0};

        // Other branches are called k, making a set B_b, with a number of members that I call n_k
        T n_k = T{0};
        for (auto k : branches) {
            if (k.id == this->id) { continue; } // Don't interact with self
            n_k += this->compute_for_branch (source_tissue, &k, C, I, J);
        }

        // Do the 1/|B_b| multiplication
        if (n_k > T{0}) {
            C = C/n_k;
            I = I/n_k;
        } // else C and I will be {0,0} still

        // Collected non-border movement components
        morph::Vector<T, 2> nonB = G * m[0] + C * m[1] + I * m[2];

        // Option for how movement is dealt with near the border; how to get axons back inside domain
        border_effect be = border_effect::gradients;

        // Border effect. A 'force' to move agents back inside the tissue boundary
        morph::Vector<T, 2> B = this->apply_border_effect (tissue, nonB, be);

        // The change in b from the model and border effects:
        morph::Vector<T, 2> db = (nonB + B * m[3]);

        // Finally add b and db to get next (uncomment to apply a speedlimit)
        this->next = b + db /* * this->speedlimit(db) */;

        // If we're penning the agent in, then check this->next and change as necessary
        if (be == border_effect::penned && this->entered == true) { this->pen_in (tissue); }
    }

};
