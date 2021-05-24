#pragma once

#include <vector>
#include <stdexcept>
#include <morph/vVector.h>
#include <morph/Vector.h>
#include "tissue.h"

// A retinotectal axon branch class. Holds current and historical positions, a preferred
// termination zone, and the algorithm for computing the next position. Could derive
// from an 'agent' base class. The branches express N Ephrin receptor types (EphA, EphB etc)
template<typename T, size_t N>
struct branch
{
    // Is branch k outside the region of the given tissue?
    bool is_outside (const branch& k, const guidingtissue<T, N>* tissue)
    {
        if (k.current[0] < tissue->x_min()
            || k.current[0] > tissue->x_max()
            || k.current[1] > tissue->y_max()
            || k.current[1] < tissue->y_min()) {
            return true;
        }
        return false;
    }

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
                // receptor 0 interacts primarily with ligand (0,1)
                // receptor 1 interacts primarily with ligand (2,3)
                // receptor 2 interacts primarily with ligand (4,5) // opposite to receptor 0
                // receptor 3 interacts primarily with ligand (6,7) // opposite to receptor 1
                G[0] = this->rcpt[0] * lg[0] + this->rcpt[1] * lg[2] + this->rcpt[2] * lg[4] + this->rcpt[3] * lg[6];
                G[1] = this->rcpt[0] * lg[1] + this->rcpt[1] * lg[3] + this->rcpt[2] * lg[5] + this->rcpt[3] * lg[7];
            } else {
                // To achieve rotation of tectum wrt retina, as in biology then:
                // receptor 0 interacts primarily with ligand (2,3)
                // receptor 1 interacts primarily with ligand (0,1)
                // receptor 2 interacts primarily with ligand (6,7) // opposite to receptor 0
                // receptor 3 interacts primarily with ligand (4,5) // opposite to receptor 1
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
            // Note, for now, just as in Simpson & Goodhill, just choose a single receptor for axon-axon interactions.

            T Q = T{0};
            if constexpr (N == 4) {
                Q = k.rcpt[0] / this->rcpt[0]
                + k.rcpt[1] / this->rcpt[1]
                + k.rcpt[2] / this->rcpt[2]
                + k.rcpt[3] / this->rcpt[3]; // forward signalling (used predominantly in paper)
            } else if constexpr (N == 2) {
                Q = k.rcpt[0] / this->rcpt[0] + k.rcpt[1] / this->rcpt[1];
            }
            //T Q = this->rcpt[0] / k.rcpt[0] +... ; // reverse signalling
            //T Q = std::max(k.rcpt[0] / this->rcpt[0], this->rcpt[0] / k.rcpt[0]); // bi-dir signalling

            kb.renormalize(); // as in paper, vector bk is a unit vector
            I += Q > this->s ? kb * W : nullvec;
            C += kb * W;
            if (W > T{0}) { n_k += T{1}; }
        }

        // Do the 1/|B_b| multiplication
        if (n_k > T{0}) {
            C = C/n_k;
            I = I/n_k;
        } // else C and I will be {0,0} still

        // Border effect. A force perpendicular to the boundary, falling off over the
        // distance r.
        morph::Vector<T, 2> B = {0, 0};

        static constexpr bool original_border_effect = false;
        if constexpr (original_border_effect == true) {
            // Test b, to see if it's near the border. Use winding number to see if it's
            // inside? Then, if outside, find out which edge it's nearest and apply that
            // force. Too complex. Instead, look at b's location. If x<0, then add component
            // to B[0]; if y<0 then add component to B[1], etc.
            //std::cout << "tissue boundary: x: " << tissue->x_min() << " to " << tissue->x_max() << std::endl;
            if (b[0] < tissue->x_min()) {
                G = {0,0};
                I = {0,0};
                C = {0,0};
                B[0] = T{1};
            } else if (b[0] < tissue->x_min()+r) {
                B[0] = T{1} * (T{1} - b[0]/r); // B[0] prop 1 - b/r
            } else if (b[0] > tissue->x_max()) {
                G = {0,0};
                I = {0,0};
                C = {0,0};
                //std::cout << "B neg max!\n";
                B[0] = T{-1};
            } else if (b[0] > (tissue->x_max()-r)) {
                //std::cout << "B neg prop!\n";
                B[0] = -(b[0] + r - T{1})/r; // B[0] prop (b+r-1)/r
            }

            if (b[1] < tissue->y_min()) {
                G = {0,0};
                I = {0,0};
                C = {0,0};
                B[1] = T{1};
            } else if (b[1] < tissue->y_min()+r) {
                B[1] = T{1} - b[1]/r;
            } else if (b[1] > tissue->y_max()) {
                G = {0,0};
                I = {0,0};
                C = {0,0};
                B[1] = T{-1};
            } else if (b[1] > (tissue->y_max()-r)) {
                B[1] = -(b[1] + r - T{1})/r; // B[1] prop (b+r-1)/r
            }
        } else {
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
        }

        constexpr bool try_limit = false;
        if constexpr (try_limit) {
            // Compute a 'speed limit' to slow movements down near to the edge of the tissue
            // and avoid the rather annoying oscillations around the edge.
            // Distance to the edge? Easy in a square. Find closest edge, find
            // perp. distance to that edge. Do this using b. Will assume that tissue is in
            // range 0 to 1. {t, b, r, l}
            morph::Vector<T, 4> distances = {1-b[1], -b[1], 1-b[0], -b[0]};
            size_t nearedge = distances.argshortest();

            morph::Vector<T, 2> db = (G * m[0] + C * m[1] + I * m[2] + B * m[3]);
#if 0
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

            // If going inwards, then no slowdown, else slowdown
            if (goingin) {
                this->next = b + db;
            } else {
#endif
                // Apply a slowdown based on distance toedge
                T toedge = distances.shortest();
                T speedlimit = T{0.1} + (T{2} / (T{1} + std::exp(-T{50} * std::abs(toedge)))) - T{1};
                std::cout << "To edge: " << toedge << " and speedlimit: " << speedlimit << std::endl;
                b += db*speedlimit;
                this->next = b;
#if 0
            }
#endif

        } else {
            b += G * m[0] + C * m[1] + I * m[2] + B * m[3];
            this->next = b;
        }
    }
    // The location and all previous locations of this branch (selected axons only).
    //morph::vVector<morph::Vector<T, 2>> path;

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

    // A sequence id
    int id = 0;
    // An id for the parent axon (there are many branches per axon)
    int aid = 0; // this is id/bpa (computed with integer division)
    // Distance parameter r is used as 2r
    static constexpr T two_r = T{0.1};
    static constexpr T r = T{0.05};
    // Signalling ratio parameter
    static constexpr T s = T{1.1};
};
