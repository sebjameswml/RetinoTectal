#pragma once

#include "branch_base.h"

#include <vector>
#include <morph/vVector.h>
#include <morph/Vector.h>
#include <morph/MathAlgo.h>
#include "tissue.h"

// A choice of schemes to deal with axon branches that try to move outside the tissue domain
enum class border_effect
{
    original, // The original scheme proposed by Simpson & Goodhill
    gradients, // A scheme which acts like there are gradients pushing axon branchs back inside tissue
    penned,     // Axon branches are 'penned in' - once inside the tissue, they can't move outside the tissue region
    penned_pushed // Like 'penned in' but also pushed away from the border with 1/r component as in original
};

// The branch which does chemoaffinity + comp + axon-axon interaction
//
// For receptor/ligand interactions we define:
// [see compute_next()]
//
template<typename T, size_t N>
struct branch : public branch_base<T,N>
{
    // Distance parameter r is used as 2r
    void setr (T _r) { this->r = _r; this->two_r = _r*T{2}; }
protected:
    T r = T{0.04};  // Note - arrangement is quite strongly dependent on this interaction radius
    T two_r = T{2}*r;
public:
    // Signalling ratio parameter for S&G-type interaction (but on 4 receptors, not 1)
    static constexpr bool s_g_interaction = false;
    static constexpr T s = T{1.1};

    // Interaction effect may have a history - after a cone's receptors are activated,
    // assume that the movement induced has some lifetime (hardcoded here as size of
    // Ihist container)
    static constexpr bool store_interaction_history = false;
    static constexpr size_t ihs = 1; // History size
    morph::Vector<morph::Vector<T, 2>, ihs> Ihist;
    size_t ihp = 0;

    // zero out the interaction history in init
    virtual void init() { this->Ihist.zero(); }

    // Estimate the ligand gradient, given the true ligand gradient
    virtual morph::Vector<T, 2*N> estimate_ligand_gradient(morph::Vector<T,2*N>& true_lgnd_grad,
                                                           morph::Vector<T,N>& true_lgnd_exp)
    {
        // This is the non-stochastic implementation
        morph::Vector<T, 2*N> lg = true_lgnd_grad;
        return lg;
    }

    // A subroutine of compute_next
    morph::Vector<T, 2> compute_chemo (const guidingtissue<T, N>* source_tissue,
                                       const guidingtissue<T, N>* tissue)
    {
        morph::Vector<T, 2> b = this->current;
        morph::Vector<T, 2> G;

        if constexpr (N==4) {
            // First, find the ligand gradients close to the current location b.  Adding
            // noise could be as simple as multiplying by a Vector of numbers drawn from
            // a Gaussian distribution with mean 1. But need to think carefully and be
            // prepared to justify it with reference to the Mortimer paper. Updated: I
            // have a possibly better way to do this. See branch_stochastic.h.
            morph::Vector<T, 8> lg0 = tissue->lgnd_grad_at (b);
            morph::Vector<T, 4> l0 = tissue->lgnd_at (b);
            morph::Vector<T, 8> lg = this->estimate_ligand_gradient (lg0, l0);

            // 4 receptors and 4 ligand gradients.
            // let receptor 0 interact primarily with ligand 0 [gradients (0,1)]
            // let receptor 1 interact primarily with ligand 1 [gradients (2,3)]
            // let receptor 2 interact primarily with ligand 2 [gradients (4,5)]
            // let receptor 3 interact primarily with ligand 3 [gradients (6,7)]
            G[0] = this->rcpt[0] * (source_tissue->forward_interactions[0] == interaction::repulsion ? -lg[0] : lg[0])
            + this->rcpt[1] * (source_tissue->forward_interactions[1] == interaction::repulsion ? -lg[2] : lg[2])
            + this->rcpt[2] * (source_tissue->forward_interactions[2] == interaction::repulsion ? -lg[4] : lg[4])
            + this->rcpt[3] * (source_tissue->forward_interactions[3] == interaction::repulsion ? -lg[6] : lg[6]);

            G[1] = this->rcpt[0] * (source_tissue->forward_interactions[0] == interaction::repulsion ? -lg[1] : lg[1])
            + this->rcpt[1] * (source_tissue->forward_interactions[1] == interaction::repulsion ? -lg[3] : lg[3])
            + this->rcpt[2] * (source_tissue->forward_interactions[2] == interaction::repulsion ? -lg[5] : lg[5])
            + this->rcpt[3] * (source_tissue->forward_interactions[3] == interaction::repulsion ? -lg[7] : lg[7]);
        } else if constexpr (N==2) {
            morph::Vector<T, 4> lg = tissue->lgnd_grad_at (b);
            G[0] = this->rcpt[0] * (source_tissue->forward_interactions[0] == interaction::repulsion ? -lg[0] : lg[0])
            + this->rcpt[1] * (source_tissue->forward_interactions[1] == interaction::repulsion ? -lg[2] : lg[2]);
            G[1] = this->rcpt[0] * (source_tissue->forward_interactions[0] == interaction::repulsion ? -lg[1] : lg[1])
            + this->rcpt[1] * (source_tissue->forward_interactions[1] == interaction::repulsion ? -lg[3] : lg[3]);
        }

        return G;
    }

    // A subroutine of compute_next
    //
    // Within a set distance of two_r, consider that receptors and ligands on
    // the two growth cones interact and produce an attractive (and/or repulsive)
    // interaction. Use the strength of this interaction to weight a vector between the
    // two cones.
    T compute_for_branch (const guidingtissue<T, N>* source_tissue, branch_base<T, N>* kp,
                          morph::Vector<T, 2>& C, morph::Vector<T, 2>& I)
    {
        // Paper deals with U_C(b,k) - the vector from branch b to branch k - and
        // sums these. However, that gives a competition term with a sign error. So
        // here, sum up the unit vectors kb.
        morph::Vector<T, 2> kb = this->current - kp->current;
        T d = kb.length();
        kb.renormalize(); // vector kb is a unit vector

        // Space-based competition, C
        //
        // W is a distance-dependent weight, which is 0 outside a distance of two_r and
        // linearly increases to 1 when d=0.
        T W = d <= this->two_r ? (T{1} - d/this->two_r) : T{0};
        C += kb * W;

        // Receptor-ligand axon-axon interaction, I
        //
        // Q collects the overall signal transmitted by ligands binding to
        // receptors. Repulsive interactions add to Q; attractive interactions make Q
        // more negative.
        T Q = T{0};
        if constexpr (s_g_interaction == true) {
            // The S & G interaction is based on the receptor expression only
            if constexpr (N == 4) {
                Q = kp->rcpt[0] / this->rcpt[0]
                + kp->rcpt[1] / this->rcpt[1]
                + kp->rcpt[2] / this->rcpt[2]
                + kp->rcpt[3] / this->rcpt[3];
            } else if constexpr (N == 2) {
                Q = kp->rcpt[0] / this->rcpt[0] + kp->rcpt[1] / this->rcpt[1];
            }
            Q /= N;
            morph::Vector<T, 2> nullvec = {0, 0};
            I += Q > this->s ? kb * W : nullvec;
        } else {
            if constexpr (N == 4) {
                // Forward signalling is activation of the receptor by the ligand. Treat as a multiplicative signal.
                Q = kp->lgnd[0] * this->rcpt[0] * (source_tissue->forward_interactions[0] == interaction::repulsion ? 1 : -1)
                + kp->lgnd[1] * this->rcpt[1] * (source_tissue->forward_interactions[1] == interaction::repulsion ? 1 : -1)
                + kp->lgnd[2] * this->rcpt[2] * (source_tissue->forward_interactions[2] == interaction::repulsion ? 1 : -1)
                + kp->lgnd[3] * this->rcpt[3] * (source_tissue->forward_interactions[3] == interaction::repulsion ? 1 : -1);
            } else if constexpr (N == 2) {
                Q = kp->lgnd[0] * this->rcpt[0] * (source_tissue->forward_interactions[0] == interaction::repulsion ? 1 : -1)
                + kp->lgnd[1] * this->rcpt[1] * (source_tissue->forward_interactions[1] == interaction::repulsion ? 1 : -1);
            }
            Q /= N;
            I += kb * Q * (d <= this->two_r ? T{1} : T{0});
        }

        // In client code, should we add to n_k or not?
        return (d <= this->two_r ? T{1} : T{0});
    }

    // A subroutine of compute_next
    morph::Vector<T, 2> apply_border_effect (const guidingtissue<T, N>* tissue, morph::Vector<T, 2>& nonB, const border_effect be)
    {
        morph::Vector<T, 2> b = this->current;

        morph::Vector<T, 2> B = {0,0};

        if (be == border_effect::original) {
            // A piecewise linear contribution to movement of agents lying outside the
            // boundary towards the inside of the tissue boundary. Includes a movement
            // away from the tissue edge up to an inner distance r from the inside
            // boundary. Also zeros the other components I (G, I, C) that contribute to
            // the agent's next position if agent is outside. This is the scheme used in
            // the Simpson & Goodhill paper and leads to quite large oscillating
            // dynamics around the tissue edge.
            if (b[0] < tissue->x_min()) {
                nonB = {0,0};
                B[0] = T{1};
            } else if (b[0] < tissue->x_min()+r) {
                B[0] = T{1} * (T{1} - b[0]/r); // B[0] prop 1 - b/r
            } else if (b[0] > tissue->x_max()) {
                nonB = {0,0};
                B[0] = T{-1};
            } else if (b[0] > (tissue->x_max()-r)) {
                B[0] = -(b[0] + r - T{1})/r; // B[0] prop (b+r-1)/r
            }

            if (b[1] < tissue->y_min()) {
                nonB = {0,0};
                B[1] = T{1};
            } else if (b[1] < tissue->y_min()+r) {
                B[1] = T{1} - b[1]/r;
            } else if (b[1] > tissue->y_max()) {
                nonB = {0,0};
                B[1] = T{-1};
            } else if (b[1] > (tissue->y_max()-r)) {
                B[1] = -(b[1] + r - T{1})/r; // B[1] prop (b+r-1)/r
            }

        } else if (be == border_effect::gradients) {
            //std::cout << "gradients\n";
            // Updated border effect. Don't change competition, interaction or chemoaffinity. B is effectively a gradient effect.
            if (b[0] < (tissue->x_min()+r)) {
                B[0] = (tissue->x_min()+r) - b[0];
            } else if (b[0] > (tissue->x_max()-r)) {
                B[0] = -(b[0] - (tissue->x_max()-r));
            }

            if (b[1] < (tissue->y_min()+r)) {
                B[1] = (tissue->y_min()+r) - b[1];
            } else if (b[1] > (tissue->y_max()-r)) {
                B[1] = -(b[1] - (tissue->y_max()-r));
            }

        } else if (be == border_effect::penned) {
            //std::cout << "penned\n";
            // Try looking at the next location and preventing the agent from moving
            // outside the domain. Have to deal with incoming agents. Give each branch
            // an 'inside' marker, which if true, says the axon has entered the domain
            // and may never leave (mu-ha-ha)
            if (this->entered == false) {
                // Use the S-G scheme before branches have entered the domain
                if (b[0] < tissue->x_min()) {
                    nonB = {0,0};
                    B[0] = T{1};
                } else if (b[0] < tissue->x_min()+r) {
                    B[0] = T{1} * (T{1} - b[0]/r); // B[0] prop 1 - b/r
                    this->entered = true;
                } else if (b[0] > tissue->x_max()) {
                    nonB = {0,0};
                    B[0] = T{-1};
                } else if (b[0] > (tissue->x_max()-r)) {
                    B[0] = -(b[0] + r - T{1})/r; // B[0] prop (b+r-1)/r
                    this->entered = true;
                } else {
                    this->entered = true;
                }

                if (b[1] < tissue->y_min()) {
                    nonB = {0,0};
                    B[1] = T{1};
                } else if (b[1] < tissue->y_min()+r) {
                    B[1] = T{1} - b[1]/r;
                    this->entered = true;
                } else if (b[1] > tissue->y_max()) {
                    nonB = {0,0};
                    B[1] = T{-1};
                } else if (b[1] > (tissue->y_max()-r)) {
                    B[1] = -(b[1] + r - T{1})/r; // B[1] prop (b+r-1)/r
                    this->entered = true;
                } else {
                    this->entered = true;
                }

                //b += nonB + B * m[3];
                //this->next = b;
            } else {
                // Border effect when near edge
                if (b[0] < tissue->x_min()+r) {
                    B[0] = T{1} * (T{1} - b[0]/r);
                } else if (b[0] > (tissue->x_max()-r)) {
                    B[0] = -(b[0] + r - T{1})/r;
                }

                if (b[1] < tissue->y_min()+r) {
                    B[1] = T{1} - b[1]/r;
                } else if (b[1] > (tissue->y_max()-r)) {
                    B[1] = -(b[1] + r - T{1})/r;
                }
            }
        }

        return B;
    }

    // Another subroutine of compute_next
    T speedlimit (const morph::Vector<T, 2>& db)
    {
        T slimit = T{1};

        morph::Vector<T, 2> b = this->current;

        //std::cout << "speedlimit\n";
        // Compute a 'speed limit' to slow movements down near to the edge of the tissue
        // and avoid the rather annoying oscillations around the edge.
        // Distance to the edge? Easy in a square. Find closest edge, find
        // perp. distance to that edge. Do this using b. Will assume that tissue is in
        // range 0 to 1. {t, b, r, l}
        morph::Vector<T, 4> distances = {1-b[1], -b[1], 1-b[0], -b[0]};
        size_t nearedge = distances.argshortest();

        bool goingin = false;
        switch (nearedge) {
        case 0: // top
        {
            if (db[1] < 0) { goingin = true; }
            break;
        }
        case 1: // bot
        {
            if (db[1] > 0) { goingin = true; }
            break;
        }
        case 2: // r
        {
            if (db[0] < 0) { goingin = true; }
            break;
        }
        case 3: // l
        {
            if (db[0] > 0) { goingin = true; }
            break;
        }
        default:
        {
            break;
        }
        }

        if (!goingin) {
            // Apply a slowdown based on distance toedge
            T toedge = distances.shortest();
            slimit = T{0.1} + (T{2} / (T{1} + std::exp(-T{50} * std::abs(toedge)))) - T{1};
            //std::cout << "To edge: " << toedge << " and speedlimit: " << slimit << std::endl;
        }

        return slimit;
    }

    // One of the border effect options is to pen the cones in.
    static constexpr bool apply_pen_in_effect = false;
    void pen_in (const guidingtissue<T, N>* tissue)
    {
        // Prevent agents from moving outside; limit this->next
        if (this->next[0] < tissue->x_min()) {
            this->next[0] = tissue->x_min();
        } else if (this->next[0] > tissue->x_max()) {
            this->next[0] = tissue->x_max();
        }
        if (this->next[1] < tissue->y_min()) {
            this->next[1] = tissue->y_min();
        } else if (this->next[1] > tissue->y_max()) {
            this->next[1] = tissue->y_max();
        }
    }

    // Compute the next position for this branch, using information from all other
    // branches and the parameters vector, m. Will also need location on target tissue,
    // to get its ephrin values.
    // rns are a set of random numbers to multiply the ligand gradient with.
    void compute_next (const std::vector<branch<T, N>>& branches,
                       const guidingtissue<T, N>* source_tissue,
                       const guidingtissue<T, N>* tissue,
                       const morph::Vector<T, 4>& m,
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

        // Other branches are called k, making a set B_b, with a number of members that I call n_k
        T n_k = T{0};
        for (auto k : branches) {
            if (k.id == this->id) { continue; } // Don't interact with self
            n_k += this->compute_for_branch (source_tissue, (&k), C, I);
        }

        // Do the 1/|B_b| multiplication to normalize C and I
        if (n_k > T{0}) { C = C/n_k; } // else C will be {0,0} still

        if constexpr (store_interaction_history == true) {
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
        //std::cout << "Movement for branch " << this->id
        //          << ": m[0]*G:" << (G*m[0])
        //          << "+ m[1]*C:" << (C*m[1])
        //          << "+ m[2]*I:" << (I*m[2]) << std::endl;
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
        if constexpr (apply_pen_in_effect == true) {
            if (be == border_effect::penned && this->entered == true) { this->pen_in (tissue); }
        }
    }
};

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
                       const morph::Vector<T, 4>& m,
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

        // Other branches are called k, making a set B_b, with a number of members that I call n_k
        T n_k = T{0};
        for (auto k : branches) {
            if (k.id == this->id) { continue; } // Don't interact with self
            n_k += this->compute_for_branch (source_tissue, &k, C, I);
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
