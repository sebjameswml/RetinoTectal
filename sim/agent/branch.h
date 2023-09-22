#pragma once

#include "branch_base.h"

#include <vector>
#include <bitset>
#include <morph/vec.h>
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
    static constexpr size_t n_flags = 3;
    // Distance parameter r is used as 2r
    T getr() { return this->r; }
    T getr_c() { return this->r_c; }
    T getr_j() { return this->r_j; }
    T getr_i() { return this->r_i; }
    void setr (T _r) { this->r = _r; this->two_r = _r*T{2}; }
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
    // Signalling ratio parameter for S&G-type (relative receptor levels) interaction (but on 4 receptors, not 1)
    T s = T{0.3};
    // Simp and Goodhill signalling ratio param
    T si = T{1.1};
    // Track rcpt-rcpt interaction sizes
    morph::vec<T, N> minses;
    morph::vec<T, N> maxes;

    // Interaction effect may have a history - after a cone's receptors are activated,
    // assume that the movement induced has some lifetime (hardcoded here as size of
    // Ihist container)
    static constexpr bool store_interaction_history = false;
    static constexpr size_t ihs = 1; // History size
    morph::vec<morph::vec<T, 2>, ihs> Ihist;
    size_t ihp = 0;

    // zero out the interaction history in init
    virtual void init() { this->Ihist.zero(); }

    // Compile time decision as to whether to add small performance hit of random noise
    static constexpr bool add_simple_gradient_noise = true;
    // Storage for randomly generated numbers
    morph::vec<T, 2*N> rn;
    // Estimate the ligand gradient, given the true ligand gradient
    virtual morph::vec<T, 2*N> estimate_ligand_gradient (morph::vec<T,2*N>& true_lgnd_grad,
                                                         morph::vec<T,N>& true_lgnd_exp)
    {
        // This is the non-stochastic implementation...
        morph::vec<T, 2*N> lg = true_lgnd_grad;
        // ...unless this if is entered into
        if (add_simple_gradient_noise == true && this->noise_gain > T{0}) {
            this->rn.randomizeN (T{0}, this->noise_gain); // Randomize normally
            lg += this->rn;
        }
        return lg;
    }

    // The simplest expression for cluster size: 1/r_A4
    T clustersz() const { return T{1} / this->rcpt0_EphA4_free; }

    // Used in compute_chemo(). Compute signal size given rect0 expression and cluster size.
    T signal (const T rcpt0, const T csize) const { return rcpt0 * csize; }

    // A subroutine of compute_next. Return G as N individual component vectors
    morph::vec<morph::vec<T, 2>, N> compute_chemo (const guidingtissue<T, N>* source_tissue,
                                                   const guidingtissue<T, N>* tissue)
    {
        morph::vec<T, 2> b = this->current;
        morph::vec<morph::vec<T, 2>, N> G;

        if constexpr (N==4) {
            // First, find the ligand gradients close to the current location b.  Adding
            // noise could be as simple as multiplying by a Vector of numbers drawn from
            // a Gaussian distribution with mean 1. But need to think carefully and be
            // prepared to justify it with reference to the Mortimer paper. Updated: I
            // have a possibly better way to do this. See branch_stochastic.h.
            morph::vec<T, 8> lg0 = tissue->lgnd_grad_at (b);
            morph::vec<T, 4> l0 = tissue->lgnd_at (b);
            morph::vec<T, 8> lg = this->estimate_ligand_gradient (lg0, l0);

            // 4 receptors and 4 ligand gradients.
            // let receptor 0 interact primarily with ligand 0 [gradients (0,1)]
            // let receptor 1 interact primarily with ligand 1 [gradients (2,3)]
            // let receptor 2 interact primarily with ligand 2 [gradients (4,5)]
            // let receptor 3 interact primarily with ligand 3 [gradients (6,7)]

            // For rcpt[0], do a special thing if interaction is 'special_EphA', which means treat EphA4 as different from EphAx
            T r0 = this->rcpt[0]; // Default case
            T r2 = this->rcpt[2];
            if (source_tissue->forward_interactions[0] == interaction::special_EphA) {
                // Compute a signal for r0 that is based on a clustersize:
                r0 = this->signal (this->rcpt[0], this->clustersz());

                // Add a 'collapse condition'. Exp. observation is that when there is
                // knocked in EphA3 AND knocked down EphA4, then the push towards caudal
                // seems to vanish - that's a collapse of the r2 receptor strength under
                // these conditions.
                bool collapse_condition = this->rcpt[0] > this->Ax_thresh && this->rcpt0_EphA4_free < this->A4_thresh;
                // If condition is met, then the effectiveness of r2 becomes very
                // small. Otherwise, r2 has a slightly lower 'strength' in this
                // "Special_EphA" condition. r2 * 0.5 verified by running model with e_eph_wt-wt.json
                r2 = collapse_condition ? this->rcpt[2]/T{10} : this->rcpt[2] * T{0.5};

            } else if (source_tissue->forward_interactions[0] == interaction::special_EphA_simple) {
                r0 = this->signal (this->rcpt[0], this->clustersz());
            }

            bool repulse0 = source_tissue->forward_interactions[0] == interaction::repulsion
                            || source_tissue->forward_interactions[0] == interaction::special_EphA
                            || source_tissue->forward_interactions[0] == interaction::special_EphA_simple;

            G[0][0] = r0            * (repulse0                                                         ? -lg[0] : lg[0]);
            G[1][0] = this->rcpt[1] * (source_tissue->forward_interactions[1] == interaction::repulsion ? -lg[2] : lg[2]);
            G[2][0] = r2            * (source_tissue->forward_interactions[2] == interaction::repulsion ? -lg[4] : lg[4]);
            G[3][0] = this->rcpt[3] * (source_tissue->forward_interactions[3] == interaction::repulsion ? -lg[6] : lg[6]);

            G[0][1] = r0            * (repulse0                                                         ? -lg[1] : lg[1]);
            G[1][1] = this->rcpt[1] * (source_tissue->forward_interactions[1] == interaction::repulsion ? -lg[3] : lg[3]);
            G[2][1] = r2            * (source_tissue->forward_interactions[2] == interaction::repulsion ? -lg[5] : lg[5]);
            G[3][1] = this->rcpt[3] * (source_tissue->forward_interactions[3] == interaction::repulsion ? -lg[7] : lg[7]);

        } else if constexpr (N==2) {
            morph::vec<T, 4> lg = tissue->lgnd_grad_at (b);
            G[0][0] = this->rcpt[0] * (source_tissue->forward_interactions[0] == interaction::repulsion ? -lg[0] : lg[0]);
            G[1][0] = this->rcpt[1] * (source_tissue->forward_interactions[1] == interaction::repulsion ? -lg[2] : lg[2]);
            G[0][1] = this->rcpt[0] * (source_tissue->forward_interactions[0] == interaction::repulsion ? -lg[1] : lg[1]);
            G[1][1] = this->rcpt[1] * (source_tissue->forward_interactions[1] == interaction::repulsion ? -lg[3] : lg[3]);
        }

        return G;
    }

    // A subroutine of compute_next
    //
    // Within a set distance of two_r, consider that receptors and ligands on the two
    // growth cones interact and produce an attractive (and/or repulsive)
    // interaction. Use the strength of this interaction to weight a vector between the
    // two cones. I is the receptor-receptor interaction of S&G, J is my receptor-ligand
    // axon-axon interaction. Argument kp is "pointer to branch k"
    //
    // New 2023: 'Gradient jamming'. Receptor-receptor interactions interfere with the
    // RGCreceptor-TectalLigand 'G'. This introduces an effect similar to that proposed
    // in 'A stochastic model for retinocollicular map development', by Koulakov and
    // Tsigankov (BMC Neuroscience 2004, http://www.biomedcentral.com/1471-2202/5/30)
    //
    // \return bitset<n_flags> With element 0: Add to C; 1: add to I; 2: add to J; 3: jam G.
    std::bitset<n_flags>
    compute_for_branch (const guidingtissue<T, N>* source_tissue, branch_base<T, N>* kp,
                        morph::vec<T, 2>& C, morph::vec<T, 2>& I, morph::vec<T, 2>& J, morph::vec<T, 2*N>& lg0,
                        morph::vec<T, N>& p_jam)
    {
        // Holds flags to say if we should add to C, I or J.
        std::bitset<n_flags> rtn;

        // Paper deals with U_C(b,k) - the vector from branch b to branch k - and
        // sums these. However, that gives a competition term with a sign error. So
        // here, sum up the unit vectors kb.
        morph::vec<T, 2> kb = this->current - kp->current;
        T d = kb.length();
        kb.renormalize(); // vector kb is a unit vector

        //////
        // G jamming
        // Possibility of jamming if d < two_r_c
        //morph::vec<T, N> _rn;
        //brng::i()->get<N>(_rn); // brng is uniform rng
        if (d <= this->two_r_c) {
            // Compare:
            // kp->rcpt
            // this->rcpt; // Each [i] of these jams a different contribution to G?
            for (size_t i = 0; i < N; ++i) {
                // 1. probability of jamming on receptor i
                // rcpt[0]: interacts with lg[0] and lg[1] (x/y) For the jamming, just use the magnitude of
                // the grad. std::sqrt(lg[2i]*lg[2i]+lg[2i+1]*lg[2i+1])
                morph::vec<T, 2> lgi = { lg0[2*i] , lg0[2*i+1] };
                // Direction from one of the elements of lg.
                // rcpt[0]: Increases in -x direction on retina
                // rcpt[1]: Increases in -y direction on retina
                // rcpt[2]: Increases in +x direction on retina
                // rcpt[3]: Increases in +y direction on retina
                // rcpt[i]  r dirn    comparative lg direction from [element]
                //    0        -x         lgi[1] or (i/2>0?1:-1)*lgisig[(i+1)%2]
                //    1        -y         lgi[0]
                //    2        +x        -lgi[1]
                //    3        +y        -lgi[0]

                // Normalised receptor expressions
                T krcpt = (kp->rcpt[i] - source_tissue->rcpt_range[i].min) / (source_tissue->rcpt_range[i].max - source_tissue->rcpt_range[i].min);
                T thisrcpt = (this->rcpt[i] - source_tissue->rcpt_range[i].min) / (source_tissue->rcpt_range[i].max - source_tissue->rcpt_range[i].min);
                // Probability of jamming
                T _p_jam = T{0.5} - T{0.5} * (krcpt - thisrcpt) * (i/2>0?1:-1) * lgi[(i+1)%2];

#if 0
                std::cout << "\n**\n(krcpt-thisrcpt)=" << (krcpt - thisrcpt) << " * (i/2>0?1:-1)=" << (i/2>0?1:-1)
                          << " * lgi[i%2="<<((i+1)%2)<<"]: " << lgi[(i+1)%2]
                          << std::endl;
                std::cout << "-----> _p_jam["<<i<<"] = " << ".5 - .5*"
                          << ((krcpt - thisrcpt) * (i/2>0?1:-1) * lgi[(i+1)%2]) << " = " << _p_jam << std::endl;
#endif
                // Force range
                _p_jam = _p_jam < T{0} ? T{0} : _p_jam;
                _p_jam = _p_jam > T{1} ? T{1} : _p_jam;

                if (_p_jam > T{1} || _p_jam < T{0}) {
                    std::stringstream ee;
                    ee << "p_jam=" << p_jam << " oor. (kp->rcpt[i]=" << kp->rcpt[i] << " - this->rcpt[i]=" << this->rcpt[i] << ") * ("
                       << "(i/2>0?1:-1) = " << (i/2>0?1:-1) << " * lgi[(i+1)%2]=" << lgi[(i+1)%2] << ")";
                    throw std::runtime_error (ee.str());
                }
                // Add/multiplicatively modify a probability of jamming (fixme - combinatorial probability!)
                p_jam[i] = p_jam[i] == T{0} ? _p_jam : p_jam[i] * _p_jam;
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////////
        // Space-based competition, C
        //
        // W is a distance-dependent weight, which is 0 outside a distance of two_r and
        // linearly increases to 1 when d=0.
        T W = d <= this->two_r_c ? (T{1} - d/this->two_r_c) : T{0};
        C += kb * W;
        rtn[0] = d <= this->two_r_c ? true : false;

        // Note: QI(QJ) collect the overall signal transmitted by ligands binding to
        // receptors. Repulsive interactions add to QI(QJ); attractive interactions make
        // QI(QJ) more negative.

        /////////////////////////////////////////////////////////////////////////////////////////////
        // Axon-axon interaction, I
        //
        // The S & G axon-axon interaction is based on the receptor expression only and
        // (guided by Reber) looks at the relative levels
        morph::vec<T,N> QI;
        QI.zero();
        W = d <= this->two_r_i ? (T{1} - d/this->two_r_i) : T{0};
        QI = this->rcpt/kp->rcpt;
        rtn[1] = false;
        for (size_t i = 0; i < N; ++i) {
            if (source_tissue->rcptrcpt_interactions[i] != interaction::null && QI[i] > this->si) {
                I += kb * W;
                rtn[1] = true;
                break; // Because for ONE of the receptor types (i in N), QI[i] > s and one is all it takes.
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////////
        // Mass action receptor-ligand axon-axon interaction, J. This is really *the
        // same as I*. Or rather, this is the *mechanism behind I* (except that the I
        // effect implies that ligand expression is exactly 1/receptor expression)
        //
        // Here, we allow a branch to determine its direction of travel based on the
        // binding of receptors and ligands.
        //
        // One important difference is the lack of a distance based weighting (cf C and
        // I effects). In early 2022, I tried adding this, and as far as I could tell,
        // it makes little difference to the result.
        morph::vec<T,N> QJar;
        for (size_t i = 0; i < N; ++i) {
            bool repulsei = source_tissue->forward_interactions[i] == interaction::repulsion
            || source_tissue->forward_interactions[i] == interaction::special_EphA
            || source_tissue->forward_interactions[i] == interaction::special_EphA_simple;
            QJar[i] = (repulsei ? T{1} : (source_tissue->forward_interactions[i] == interaction::attraction ? T{-1} : T{0}));
        }
        QJar *= kp->lgnd * this->rcpt;

        bool subthreshold = false;
        if (d <= this->two_r_j) {

            // Add to competition if any of the mass-action interactions is true.
            // *Note1: use <= operator which is true if EVERY member of the vec is <=
            // threshold. *Note2: There's NO QJar.abs() here, so attractive elements of
            // QJar (which are -ve) will NEVER trigger the threshold, even if they are
            // strong.
            if ((subthreshold = QJar <= this->s) == false) {
                //std::cout << "subthreshold false so OVER threshold for interaction...\n";

                // Ideas that failed for a weighting:
                // W = QJar.sum()/N; // Sum the interactions to modulate the competition strength?
                // W = QJar.sum()/N * (T{1} - d/this->two_r_j); // Sum and also distance based?
                // W = QJar.longest(); // Use the maximum suprathreshold interaction element?

                // Just a simple, distance based weight, assuming a repulsive interaction (i.e. this is +ve)
                W = T{1} - d/this->two_r_j;
                J += kb * W;
                rtn[2] = true;

            } else {
                //std::cout << "subthreshold true (no interaction)...\n";
                rtn[2] = false;
            }
        } else {
            //std::cout << "d > 2rj (no interaction): " << d << " cf. two_r_j " << this->two_r_j << std::endl;
            rtn[2] = false;
        }
        return rtn;
    }

    // A subroutine of compute_next
    morph::vec<T, 2> apply_border_effect (const guidingtissue<T, N>* tissue, morph::vec<T, 2>& nonB)
    {
        morph::vec<T, 2> B = {0,0};

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
                       const morph::vec<T, 5>& m,
                       const morph::vec<T, 2*N>& rns)
    {
        // Current location is named b
        morph::vec<T, 2> b = this->current;
        // Chemotaxis, graded by origin position (i.e termination zone) of each retinal axon
        morph::vec<morph::vec<T, 2>, N> G = this->compute_chemo (source_tissue, tissue);

        // Competition, C, and Axon-axon interactions, I&J, computed during the same loop
        // over the other branches
        morph::vec<T, 2> C = {0, 0};
        morph::vec<T, 2> I = {0, 0};
        morph::vec<T, 2> J = {0, 0};

        // Other branches are called k, making a set (3 sets) B_b, with a number of members that I call n_k etc
        T n_k = T{0};
        T n_ki = T{0};
        T n_kj = T{0};

        // Use ligand gradient in compute_for_branch
        morph::vec<T, 2*N> lg0 = tissue->lgnd_grad_at (b);
        lg0 /= tissue->lgnd_grad_range[0].max; // Can assume max lgnd grad max for all is same. lg0 normalised

        // Build up a probability of jamming occurring from all the b-k interactions
        morph::vec<T, N> p_jam;
        p_jam.zero();
        for (auto k : branches) {
            if (k.id == this->id) { continue; } // Don't interact with self
            std::bitset<n_flags> cij_added = this->compute_for_branch (source_tissue, (&k), C, I, J, lg0, p_jam);
            n_k += cij_added[0] ? T{1} : T{0};
            n_ki += cij_added[1] ? T{1} : T{0};
            n_kj += cij_added[2] ? T{1} : T{0};
        }
        // Apply the probability of jamming
        morph::vec<T, N> r_jam;
        brng::i()->get<N>(r_jam);
        morph::vec<T, N> jam;
        for (size_t i = 0; i < N; ++i) {
            jam[i] = r_jam[i] < p_jam[i] ? T{0} : T{1};
            G[i] *= jam[i];
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
        morph::vec<T, 2> nonB = G.sum() * m[0] + J * m[1] + I * m[2]  + C * m[3];

        // Border effect. A 'force' to move agents back inside the tissue boundary
        morph::vec<T, 2> B = this->apply_border_effect (tissue, nonB);

        // The change in b from the model and border effects:
        morph::vec<T, 2> db = (nonB + B * m[4]);

        // Finally add b and db to get next
        this->next = b + db;
    }
};
