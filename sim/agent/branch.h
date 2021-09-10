#pragma once

#include "branch_base.h"

#include <vector>
#include <bitset>
#include <morph/vVector.h>
#include <morph/Vector.h>
#include <morph/MathAlgo.h>
#include "tissue.h"

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
    T getr_c() { return this->r_c; }
    T getr_j() { return this->r_j; }
    T getr_i() { return this->r_i; }
    void setr_c (T _r) { this->r_c = _r; this->two_r_c = _r*T{2}; }
    void setr_j (T _r) { this->r_j = _r; this->two_r_j = _r*T{2}; }
    void setr_i (T _r) { this->r_i = _r; this->two_r_i = _r*T{2}; }
protected:
    T r = T{0.04};   // Actual radius of a growth cone (may be used for visualisation)
    T two_r = T{2}*r;
    T r_c = T{0.04};  // competition interaction distance (arrangement is quite strongly dependent on this interaction radius)
    T two_r_c = T{2}*r;
    T r_j = T{0.04}; // receptor-ligand interaction distance for axon-axon interactions
    T two_r_j = T{2}*r;
    T r_i = T{0.04}; // receptor-receptor interaction distance for axon-axon interactions
    T two_r_i = T{2}*r;
public:
    // Signalling ratio parameter for S&G-type interaction (but on 4 receptors, not 1)
    static constexpr bool s_g_interaction = false;
    T s = T{0.3};

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
    std::bitset<3>
    compute_for_branch (const guidingtissue<T, N>* source_tissue, branch_base<T, N>* kp,
                        morph::Vector<T, 2>& C, morph::Vector<T, 2>& I, morph::Vector<T, 2>& J)
    {
        // Holds flags to say if we should add to C, I or J.
        std::bitset<3> rtn;

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
        T W = d <= this->two_r_c ? (T{1} - d/this->two_r_c) : T{0};
        C += kb * W;
        rtn[0] = d <= this->two_r_c ? true : false;

        // QI/J collects the overall signal transmitted by ligands binding to
        // receptors. Repulsive interactions add to Q; attractive interactions make Q
        // more negative.
#define SIMPSON_GOODHILL_LIKE_RCPTRCPT_INTERACTIONS 1
#ifdef SIMPSON_GOODHILL_LIKE_RCPTRCPT_INTERACTIONS
        // The S & G axon-axon interaction is based on the receptor expression only and
        // (guided by Reber) looks at the relative levels
        T QI = T{0};
        W = d <= this->two_r_i ? (T{1} - d/this->two_r_i) : T{0};
        if constexpr (N == 4) {
            QI = kp->rcpt[0] / this->rcpt[0];
            //+ kp->rcpt[1] / this->rcpt[1]
            //+ kp->rcpt[2] / this->rcpt[2]
            //+ kp->rcpt[3] / this->rcpt[3];
        } else if constexpr (N == 2) {
            QI = kp->rcpt[0] / this->rcpt[0] + kp->rcpt[1] / this->rcpt[1];
        }
        //QI /= N;
        morph::Vector<T, 2> nullvec = {0, 0};
        I += QI > this->s ? kb * W : nullvec; // FIXME FIXME: Need to count up number of interactions and divide by this number (like for comp)
        rtn[1] = (QI > this->s && W > T{0}) ? true : false;

#else
        T QI = T{0};
        if constexpr (N == 4) {
            // Receptor-receptor activation according to rcptrcpt_interactions
            QI = kp->rcpt[0] * this->rcpt[0] * (source_tissue->rcptrcpt_interactions[0] == interaction::repulsion ? 1 :
                                                (source_tissue->rcptrcpt_interactions[0] == interaction::attraction ? -1 : 0))
            + kp->rcpt[1] * this->rcpt[1] * (source_tissue->rcptrcpt_interactions[1] == interaction::repulsion ? 1 :
                                             (source_tissue->rcptrcpt_interactions[1] == interaction::attraction ? -1 : 0))
            + kp->rcpt[2] * this->rcpt[2] * (source_tissue->rcptrcpt_interactions[2] == interaction::repulsion ? 1 :
                                             (source_tissue->rcptrcpt_interactions[2] == interaction::attraction ? -1 : 0))
            + kp->rcpt[3] * this->rcpt[3] * (source_tissue->rcptrcpt_interactions[3] == interaction::repulsion ? 1 :
                                             (source_tissue->rcptrcpt_interactions[3] == interaction::attraction ? -1 : 0));
        } else if constexpr (N == 2) {
            QI = kp->rcpt[0] * this->rcpt[0] * (source_tissue->rcptrcpt_interactions[0] == interaction::repulsion ? 1 :
                                                (source_tissue->rcptrcpt_interactions[0] == interaction::attraction ? -1 : 0))
            + kp->rcpt[1] * this->rcpt[1] * (source_tissue->rcptrcpt_interactions[1] == interaction::repulsion ? 1 :
                                             (source_tissue->rcptrcpt_interactions[1] == interaction::attraction ? -1 : 0));
        }
        QI /= N;
        I += kb * QI * (d <= this->two_r_i ? T{1} : T{0});
        rtn[1] = (d <= this->two_r_i && QI > T{0}) ? true : false;
#endif
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
        J += kb * QJ * (d <= this->two_r_j ? T{1} : T{0});
        rtn[2] = (d <= this->two_r_j && QJ > T{0}) ? true : false;

        return rtn;
    }

    // A subroutine of compute_next
    morph::Vector<T, 2> apply_border_effect (const guidingtissue<T, N>* tissue, morph::Vector<T, 2>& nonB)
    {
        morph::Vector<T, 2> B = {0,0};

        // Gradient based border effect. Don't change competition, interaction or chemoaffinity. B is effectively a gradient effect.
        if (this->current[0] < (tissue->x_min()+r)) {
            B[0] = (tissue->x_min()+r) - this->current[0];
        } else if (this->current[0] > (tissue->x_max()-r)) {
            B[0] = -(this->current[0] - (tissue->x_max()-r));
        }

        if (this->current[1] < (tissue->y_min()+r)) {
            B[1] = (tissue->y_min()+r) - this->current[1];
        } else if (this->current[1] > (tissue->y_max()-r)) {
            B[1] = -(this->current[1] - (tissue->y_max()-r));
        }

        return B;
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
        //std::cout << "C had " << n_k << " interactions, I had " << n_ki << " and J had " << n_kj << std::endl;
        //std::cout << "J=" << J << " and J/" << n_kj << "=" << (J/n_kj) << std::endl;

        // Do the 1/|B_b| multiplication to normalize C and I(!!)
        if (n_k > T{0}) { C = C/n_k; } // else C will be {0,0} still
        I = n_ki > T{0} ? I/n_ki : I;
        J = n_kj > T{0} ? J/n_kj : J;
        //std::cout << "C=" << C << std::endl;

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

#if 0
        morph::Vector<T, 2> R = {0, 0}; // A little random movement, too
        brng::i()->get(R);
#endif
        morph::Vector<T, 2> nonB = G * m[0] + J * m[1] + I * m[2]  + C * m[3]; // + R * m[5];

        // Border effect. A 'force' to move agents back inside the tissue boundary
        morph::Vector<T, 2> B = this->apply_border_effect (tissue, nonB);

        // The change in b from the model and border effects:
        morph::Vector<T, 2> db = (nonB + B * m[4]);

        // Finally add b and db to get next
        this->next = b + db;
    }
};
