#pragma once

#include <vector>
#include <stdexcept>
#include <morph/vVector.h>
#include <morph/Vector.h>
#include <morph/MathConst.h>
#include "tissue.h"

// A choice of schemes to deal with axon branches that try to move outside the tissue domain
enum class border_effect
{
    original, // The original scheme proposed by Simpson & Goodhill
    gradients, // A scheme which acts like there are gradients pushing axon branchs back inside tissue
    penned,     // Axon branches are 'penned in' - once inside the tissue, they can't move outside the tissue region
    penned_pushed // Like 'penned in' but also pushed away from the border with 1/r component as in original
};

// A retinotectal axon branch class. Holds current and historical positions, a preferred
// termination zone, and the algorithm for computing the next position. Could derive
// from an 'agent' base class. The branches express N Ephrin receptor types (EphA, EphB etc)
template<typename T, size_t N>
struct branch_base
{
#if 0
    // Is branch k outside the region of the given tissue?
    bool is_outside (const branch_base* k, const guidingtissue<T, N>* tissue)
    {
        if (k->current[0] < tissue->x_min()
            || k->current[0] > tissue->x_max()
            || k->current[1] > tissue->y_max()
            || k->current[1] < tissue->y_min()) {
            return true;
        }
        return false;
    }
#endif

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

    // By default, axons reckoned to be outside tissue. Can also be used to mean "active".
    bool entered = false;

    // A sequence id
    int id = 0;
    // An id for the parent axon (there are many branches per axon)
    int aid = 0; // this is id/bpa (computed with integer division)
};

// A branch model inspired by Gebhardt et al 2012 (Balancing of ephrin/Eph forward and
// reverse signalling as the driving force of adaptive topographic mapping)
template<typename T, size_t N>
struct branch_geb : public branch_base<T,N>
{
    // The cone radius, hard coded
    static constexpr T r = T{0.05};

    // For random angle choice
    morph::RandUniform<T> rng;

    // To achieve rotation of tectum wrt retina, as in biology then, if N==4:
    // receptor 0 interacts primarily with ligand 1
    // receptor 1 interacts primarily with ligand 0
    // receptor 2 interacts primarily with ligand 3 // opposite to receptor 0
    // receptor 3 interacts primarily with ligand 2 // opposite to receptor 1
    void compute_next (const std::vector<branch_geb<T, N>>& branches,
                       const guidingtissue<T, N>* tissue,
                       const morph::Vector<T, 4>& m)
    {
        // Sums
        morph::Vector<T, N> LFRT; // reverse fibre-target signal (for N receptors, N ligands)
        morph::Vector<T, N> LFRF; // reverse cis (fibre-self) signal
        morph::Vector<T, N> LFRf; // reverse fibre - other fibres signal

        morph::Vector<T, N> RFLT; // forward fibre-target signal
        morph::Vector<T, N> RFLF; // forward cis (fibre-self) signal
        morph::Vector<T, N> RFLf; // forward fibre - other fibres signal

        // Will need to run through other branches to determine the R_f/L_f
        morph::vVector<morph::Vector<T, N>> Rf(tissue->posn.size());
        morph::vVector<morph::Vector<T, N>> Lf(tissue->posn.size());
        for (auto k : branches) {
            if (k.id == this->id) { continue; }
            // From location of k (k.current) figure out an idx on the target tissue and
            // add to Rf[idx] and Lf[idx].
            size_t idx = std::numeric_limits<size_t>::max();
            T dist = T{1e9};
            for (size_t i = 0; i < tissue->posn.size(); ++i) {
                if ((tissue->posn[i] - k.current).length() < dist) {
                    dist = (tissue->posn[i] - k.current).length();
                    idx = i;
                }
            }
            //std::cout << "Current location " << k.current << " adds to field index " << idx << std::endl;
            Rf[idx] += k.rcpt;
            Lf[idx] += k.lgnd;
        }

        // Branch has this->rcpt and this->lgnd corresponding to Gebhardt's R_F and L_F
        // Tissue has rcpt and lgnd corresponding to Gebhardt's R_T and L_T
        // Will need to compute R_f L_f accordint to Gebhardt's recipe.

        // First, consider own location, this->current. Consider own radius. Run through
        // tissue locations summing up R/L_F and R/L_T interactions at each that lies
        // under the cone.

        // Use this->current and a range of adjacent locations on a circle to determine next location.
        morph::vVector<morph::Vector<T,2>> candposns;
        candposns.push_back (this->current);
        // Choose a random direction to move a distance r (2r?) and evaluate that (G' in Gebhardt et al)
        T angle = this->rng.get() * T{morph::TWO_PI_D};
        morph::Vector<T,2> randvec = { this->r * std::sin(angle), this->r * std::cos(angle) };
        candposns.push_back (this->current+randvec);

        //std::cout << "Current and candidate positions: " << candposns[0] << ", " << candposns[1] << std::endl;

        // Candidate probabilities
        morph::vVector<T> pdG (2, T{0});;

        for (size_t j = 0; j<candposns.size(); ++j) {

            // Zero the sums
            LFRT.zero();
            LFRF.zero();
            LFRf.zero();

            RFLT.zero();
            RFLF.zero();
            RFLf.zero();

            const T a = T{100};
            for (unsigned int i = 0; i < tissue->posn.size(); ++i) {
                morph::Vector<T,2> pvec = tissue->posn[i] - candposns[j];
                T Gaussian = std::exp (-(a*(pvec[0]*pvec[0]) + a*(pvec[1]*pvec[1])));
                //std::cout << "Gaussian = " << Gaussian << std::endl;
                //if (pvec.length() <= this->r) { // Actually, in the Gebhardt work, r is defined by the Gaussian getting to 0.01
                if (Gaussian >= T{0.01}) {
                    // then posn[i] is under our cone!
                    morph::Vector<T, N> LgndGauss = this->lgnd * Gaussian;
                    morph::Vector<T, N> RcptGauss = this->rcpt * Gaussian;

                    // My scheme is that receptor 0 interacts with ligand 1 etc. to achieve
                    // rotation of tectum wrt retina.
                    morph::Vector<T, N> LgndGaussR = LgndGauss;
                    morph::Vector<T, N> RcptGaussR = RcptGauss;
                    morph::Vector<T, N> tissueRcptR = tissue->rcpt[i];
                    morph::Vector<T, N> tissueLgndR = tissue->lgnd[i];
                    morph::Vector<T, N> RfR = Rf[i];
                    morph::Vector<T, N> LfR = Lf[i];
                    LgndGaussR.rotate_pairs();
                    RcptGaussR.rotate_pairs();
                    tissueRcptR.rotate_pairs();
                    tissueLgndR.rotate_pairs();
                    RfR.rotate_pairs();
                    LfR.rotate_pairs();

                    // Reverse signalling
                    //std::cout << "Fibre lgnd * Gaussian: " << LgndGauss << std::endl;
                    LFRT += (LgndGauss * tissueRcptR); // operator* will give element-wise multiplication
                    LFRF += (LgndGauss * RcptGaussR);
                    LFRf += (LgndGauss * RfR);
                    // Forward signalling
                    //std::cout << "Fibre receptor * Gaussian: " << RcptGauss << std::endl;
                    RFLT += (RcptGauss * tissueLgndR);
                    RFLF += (RcptGauss * LgndGaussR);
                    RFLf += (RcptGauss * LfR);
                }
            }

            // Now calculate the 'potential for moving', G resulting from each receptor-ligand pair.
#if 0
            std::cout << "LFRT: " << LFRT
                      << ", LFRF: " << LFRF
                      << ", LFRf: " << LFRf
                      << ", RFLT: " << RFLT
                      << ", RFLF: " << RFLF
                      << ", RFLf: " << RFLf << std::endl;
#endif
            morph::Vector<T, N> G_num = LFRT + LFRF + LFRf;
            morph::Vector<T, N> G_denom = RFLT + RFLF + RFLf;
            morph::Vector<T, N> G;
            if (G_denom.sum()) {
                G = (G_num / G_denom);
            } else {
                G.zero();
            }

            // Final number is abs(log(G.sum())). This is G_{F,i} from paper
            T GFi = std::abs(std::log(G.sum()));
            std::cout << "G: " << G << ", abs(log(G.sum()): " << GFi << std::endl;
            T sigma = 0.12;
            // Store candidate G as probability density:
            pdG[j] = (T{1}/(sigma*T{morph::TWO_PI_D})) * std::exp ( - GFi * GFi / (2*sigma*sigma));
            std::cout << "pdG[candidate] = " << pdG[j] <<std::endl;
        }

        // Now, probabilistically determine what this->next will be. Can I do this
        // continuously? No, not without significant computation.  In Gebhardt, they
        // look at the 9 adjacent squares and compute G for each, then choose one
        // probabilistically, weighted by G. I'll do the same.

        // pdG[0] is the current location's potential probability density (pd(G))
        // pdG[1] is the candidate location's potential
        T p_change = pdG[1] / pdG.sum();
        T rn = this->rng.get();
        this->next = rn < p_change ? candposns[1] : this->current;
    }
};

// The branch which does chemoaffinity+comp+(forward)interaction
template<typename T, size_t N>
struct branch : public branch_base<T,N>
{
    // Distance parameter r is used as 2r
    static constexpr T two_r = T{0.1};
    static constexpr T r = T{0.05};
    // Signalling ratio parameter
    static constexpr T s = T{1.1};

    static constexpr bool biological_tectum = true;
    // Compute the next position for this branch, using information from all other
    // branches and the parameters vector, m. Will also need location on target tissue,
    // to get its ephrin values.
    void compute_next (const std::vector<branch<T, N>>& branches,
                       const guidingtissue<T, N>* tissue,
                       const morph::Vector<T, 4>& m)
    {
        // Current location is named b
        morph::Vector<T, 2> b = this->current;

        // Chemoaffinity, graded by origin position (i.e termination zone) of each retinal axon
        morph::Vector<T, 2> G;
        if constexpr (N==4) {
            // First, find the ligand gradients close to the current location b.
            morph::Vector<T, 8> lg = tissue->lgnd_grad_at (b);
            // 4 receptors and 4 ligand gradients. The 4 receptors are in order: r, u, l, d
            if constexpr (biological_tectum == false) {
                // If no rotation of tectum wrt retina then:
                // receptor 0 interacts primarily with ligand gradients (0,1)
                // receptor 1 interacts primarily with ligand gradients (2,3)
                // receptor 2 interacts primarily with ligand gradients (4,5) // opposite to receptor 0
                // receptor 3 interacts primarily with ligand gradients (6,7) // opposite to receptor 1
                G[0] = this->rcpt[0] * lg[0] + this->rcpt[1] * lg[2] + this->rcpt[2] * lg[4] + this->rcpt[3] * lg[6];
                G[1] = this->rcpt[0] * lg[1] + this->rcpt[1] * lg[3] + this->rcpt[2] * lg[5] + this->rcpt[3] * lg[7];
            } else {
                // To achieve rotation of tectum wrt retina, as in biology then:
                // receptor 0 interacts primarily with ligand gradients (2,3)
                // receptor 1 interacts primarily with ligand gradients (0,1)
                // receptor 2 interacts primarily with ligand gradients (6,7) // opposite to receptor 0
                // receptor 3 interacts primarily with ligand gradients (4,5) // opposite to receptor 1
                G[0] = this->rcpt[0] * lg[2] + this->rcpt[1] * lg[0] + this->rcpt[2] * lg[6] + this->rcpt[3] * lg[4];
                G[1] = this->rcpt[0] * lg[3] + this->rcpt[1] * lg[1] + this->rcpt[2] * lg[7] + this->rcpt[3] * lg[5];
            }
        } else if constexpr (N==2) {
            morph::Vector<T, 4> lg = tissue->lgnd_grad_at (b);
            G[0] = this->rcpt[0] * lg[0] + this->rcpt[1] * lg[2];
            G[1] = this->rcpt[1] * lg[1] + this->rcpt[1] * lg[3];
        }

        // Competition, C, and Axon-axon interactions, I, computed during the same loop
        // over the other branches
        morph::Vector<T, 2> C = {0, 0};
        morph::Vector<T, 2> I = {0, 0};
        morph::Vector<T, 2> nullvec = {0, 0}; // null vector
        // Other branches are called k, making a set B_b, with a number of members that I call n_k
        T n_k = T{0};
        for (auto k : branches) {
            if (k.id == this->id) { continue; } // Don't interact with self
            //if (this->is_outside (k, tissue)) { continue; } // Don't interact with the branches outside the tissue
            // Paper deals with U_C(b,k) - the vector from branch b to branch k - and
            // sums these. However, that gives a competition term with a sign error. So
            // here, sum up the unit vectors kb.
            morph::Vector<T, 2> kb = b - k.current;
            T d = kb.length();
            T W = d <= this->two_r ? (T{1} - d/this->two_r) : T{0};
            T Q = T{0};
            if constexpr (N == 4) {
                Q = k.rcpt[0] / this->rcpt[0]
                + k.rcpt[1] / this->rcpt[1]
                + k.rcpt[2] / this->rcpt[2]
                + k.rcpt[3] / this->rcpt[3]; // forward signalling (used predominantly in paper)
            } else if constexpr (N == 2) {
                Q = k.rcpt[0] / this->rcpt[0] + k.rcpt[1] / this->rcpt[1];
            }
            Q /= N; // to average the effect of each receptor type

            //T Q = this->rcpt[0] / k.rcpt[0] +... ; // reverse signalling
            //T Q = std::max(k.rcpt[0] / this->rcpt[0], this->rcpt[0] / k.rcpt[0]); // bi-dir signalling

            kb.renormalize(); // as in paper, vector bk is a unit vector
            I += Q > this->s ? kb * W : nullvec;
            C += kb * W;
            //if (W > T{0}) { n_k += T{1}; }
            n_k += W > T{0} ? T{1} : T{0};
        }

        // Do the 1/|B_b| multiplication
        if (n_k > T{0}) {
            C = C/n_k;
            I = I/n_k;
        } // else C and I will be {0,0} still

        // Collected non-border movement components
        morph::Vector<T, 2> nonB = G * m[0] + C * m[1] + I * m[2];

        // Border effect. A 'force' to move agents back inside the tissue boundary
        morph::Vector<T, 2> B = {0,0};

        // Two options for how movement is dealt with near the border. First is how to get axons back inside domain
        border_effect be = border_effect::gradients;
        // Second is an option to slow movement down near border to reduce the amplitude of the border dynamics
        constexpr bool limit_movement_near_border = false;

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
            if (b[0] < (tissue->x_min()+(2*r))) {
                B[0] = (tissue->x_min()+(2*r)) - b[0];
            } else if (b[0] > (tissue->x_max()-(2*r))) {
                B[0] = -(b[0] - (tissue->x_max()-(2*r)));
            }

            if (b[1] < (tissue->y_min()+(2*r))) {
                B[1] = (tissue->y_min()+(2*r)) - b[1];
            } else if (b[1] > (tissue->y_max()-(2*r))) {
                B[1] = -(b[1] - (tissue->y_max()-(2*r)));
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

        // The change in b from the model and border effects:
        morph::Vector<T, 2> db = (nonB + B * m[3]);

        T speedlimit = T{1};
        if constexpr (limit_movement_near_border) {
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
                speedlimit = T{0.1} + (T{2} / (T{1} + std::exp(-T{50} * std::abs(toedge)))) - T{1};
                //std::cout << "To edge: " << toedge << " and speedlimit: " << speedlimit << std::endl;
            }
        }

        // Finally add b and db to get next (with speedlimit)
        this->next = b + db*speedlimit;

        // If we're penning the agent in, then check this->next and change as necessary
        if (be == border_effect::penned && this->entered == true) {
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
    }
};
