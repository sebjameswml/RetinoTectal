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
    T getr() { return this->r; }
    T getrc() { return this->rc; }
    T getrrl() { return this->rrl; }
    T getrrr() { return this->rrr; }
    void setrc (T _r) { this->rc = _r; this->two_rc = _r*T{2}; }
    void setrrl (T _r) { this->rrl = _r; this->two_rrl = _r*T{2}; }
    void setrrr (T _r) { this->rrr = _r; this->two_rrr = _r*T{2}; }
protected:
    T r = T{0.04};   // Actual radius of a growth cone
    T two_r = T{2}*r;
    T rc = T{0.04};  // competition interaction distance (arrangement is quite strongly dependent on this interaction radius)
    T two_rc = T{2}*r;
    T rrl = T{0.04}; // receptor-ligand interaction distance
    T two_rrl = T{2}*r;
    T rrr = T{0.04}; // receptor-receptor interaction distance
    T two_rrr = T{2}*r;
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
    // Within a set distance of two_r, consider that receptors and ligands on the two
    // growth cones interact and produce an attractive (and/or repulsive)
    // interaction. Use the strength of this interaction to weight a vector between the
    // two cones. I is the receptor-receptor interaction of S&G, J is my receptor-ligand
    // axon-axon interaction.
    T compute_for_branch (const guidingtissue<T, N>* source_tissue, branch_base<T, N>* kp,
                          morph::Vector<T, 2>& C, morph::Vector<T, 2>& I, morph::Vector<T, 2>& J)
    {
        // Paper deals with U_C(b,k) - the vector from branch b to branch k - and
        // sums these. However, that gives a competition term with a sign error. So
        // here, sum up the unit vectors kb.
        morph::Vector<T, 2> kb = this->current - kp->current;
        T d = kb.length();
        kb.renormalize(); // vector kb is a unit vector

        if constexpr (debug_compute_branch == true) {
            std::cout << "d=" << d << ", 2rrl=" << this->two_rrl << ", 2rrr=" << this->two_rrr <<  std::endl;
        }
        // Space-based competition, C
        //
        // W is a distance-dependent weight, which is 0 outside a distance of two_r and
        // linearly increases to 1 when d=0.
        T W = d <= this->two_rc ? (T{1} - d/this->two_rc) : T{0};
        C += kb * W;

        // QI/J collects the overall signal transmitted by ligands binding to
        // receptors. Repulsive interactions add to Q; attractive interactions make Q
        // more negative.

        // The S & G axon-axon interaction is based on the receptor expression only and
        // (guided by Reber) looks at the relative levels
        T QI = T{0};
        W = d <= this->two_rrr ? (T{1} - d/this->two_rrr) : T{0};
        if constexpr (N == 4) {
            QI = kp->rcpt[0] / this->rcpt[0]
            + kp->rcpt[1] / this->rcpt[1]
            + kp->rcpt[2] / this->rcpt[2]
            + kp->rcpt[3] / this->rcpt[3];
        } else if constexpr (N == 2) {
            QI = kp->rcpt[0] / this->rcpt[0] + kp->rcpt[1] / this->rcpt[1];
        }
        QI /= N;
        morph::Vector<T, 2> nullvec = {0, 0};
        if constexpr (debug_compute_branch == true) {
            std::cout << "QI=" << QI << ", this->s = " << this->s << std::endl;
        }
        I += QI > this->s ? kb * W : nullvec;

        // Receptor-ligand axon-axon interaction, J
        T QJ = T{0};
        if constexpr (N == 4) {
            // Forward signalling is activation of the receptor by the ligand. Treat as a multiplicative signal.
            QJ = kp->lgnd[0] * this->rcpt[0] * (source_tissue->forward_interactions[0] == interaction::repulsion ? 1 : -1)
            + kp->lgnd[1] * this->rcpt[1] * (source_tissue->forward_interactions[1] == interaction::repulsion ? 1 : -1)
            + kp->lgnd[2] * this->rcpt[2] * (source_tissue->forward_interactions[2] == interaction::repulsion ? 1 : -1)
            + kp->lgnd[3] * this->rcpt[3] * (source_tissue->forward_interactions[3] == interaction::repulsion ? 1 : -1);
            // Alternative? Could have a kind of competition for the movement, rather than just adding 'em all up.
        } else if constexpr (N == 2) {
            QJ = kp->lgnd[0] * this->rcpt[0] * (source_tissue->forward_interactions[0] == interaction::repulsion ? 1 : -1)
            + kp->lgnd[1] * this->rcpt[1] * (source_tissue->forward_interactions[1] == interaction::repulsion ? 1 : -1);
        }
        QJ /= N;
        J += kb * QJ * (d <= this->two_rrl ? T{1} : T{0});

        if constexpr (debug_compute_branch == true) {
            std::cout << "For this branch, I=" << I << ", and J=" << J << std::endl;
        }

        // In client code, should we add to n_k or not? (only used for competition, hence d <= two_rc)
        return (d <= this->two_rc ? T{1} : T{0});
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
                       const morph::Vector<T, 5>& m,
                       const morph::Vector<T, 2*N>& rns)
    {
        // Current location is named b
        morph::Vector<T, 2> b = this->current;
        // Chemoaffinity, graded by origin position (i.e termination zone) of each retinal axon
        morph::Vector<T, 2> G = this->compute_chemo (source_tissue, tissue);

        // Competition, C, and Axon-axon interactions, I&J, computed during the same loop
        // over the other branches
        morph::Vector<T, 2> C = {0, 0};
        morph::Vector<T, 2> I = {0, 0};
        morph::Vector<T, 2> J = {0, 0};

        // Other branches are called k, making a set B_b, with a number of members that I call n_k
        T n_k = T{0};
        for (auto k : branches) {
            if (k.id == this->id) { continue; } // Don't interact with self
            n_k += this->compute_for_branch (source_tissue, (&k), C, I, J);
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
        morph::Vector<T, 2> nonB = G * m[0] + C * m[1] + I * m[2] + J * m[3];

        // Option for how movement is dealt with near the border; how to get axons back inside domain
        border_effect be = border_effect::gradients;

        // Border effect. A 'force' to move agents back inside the tissue boundary
        morph::Vector<T, 2> B = this->apply_border_effect (tissue, nonB, be);

        // The change in b from the model and border effects:
        morph::Vector<T, 2> db = (nonB + B * m[4]);

        // Finally add b and db to get next (uncomment to apply a speedlimit)
        this->next = b + db /* * this->speedlimit(db) */;

        // If we're penning the agent in, then check this->next and change as necessary
        if constexpr (apply_pen_in_effect == true) {
            if (be == border_effect::penned && this->entered == true) { this->pen_in (tissue); }
        }
    }
};
