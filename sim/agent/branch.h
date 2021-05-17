#pragma once

#include <vector>
#include <morph/vVector.h>
#include <morph/Vector.h>
#include "tissue.h"

// A retinotectal axon branch class. Holds current and historical positions, a preferred
// termination zone, and the algorithm for computing the next position. Could derive
// from an 'agent' base class. The branches express N Ephrin receptor types (EphA, EphB etc)
template<typename T, size_t N>
struct branch
{
    // Compute the next position for this branch, using information from all other
    // branches and the parameters vector, m. Will also need location on target tissue,
    // to get its ephrin values.
    void compute_next (const std::vector<branch<T, N>>& branches,
                       const guidingtissue<T, N>* tissue,
                       const morph::Vector<T, 4>& m)
    {
        // Current location is named b
        morph::Vector<T, 2> b = path.back();

        // Chemoaffinity, graded by origin position (i.e termination zone) of each retinal axon
        morph::Vector<T, 2> G;
        if constexpr (N==4) {
            // First, find the ligand gradients close to the current location b.
            morph::Vector<T, 4> lg = tissue->lgnd_grad_at (b);
            // With 4 receptor/ligand pairs, the x component is zeroth minus the 2nd
            G[0] = this->rcpt[0] * lg[0] - this->rcpt[2] * lg[2];
            // and y component is the first minus the 3rd
            G[1] = this->rcpt[1] * lg[1] - this->rcpt[3] * lg[3];
        } else if constexpr (N==2) {
            morph::Vector<T, 2> lg = tissue->lgnd_grad_at (b);
            G[0] = this->rcpt[0] * lg[0];
            G[1] = this->rcpt[1] * lg[1];
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
            // Paper deals with U_C(b,k) - the vector from branch b to branch k - and
            // sums these. However, that gives a competition term with a sign error. So
            // here, sum up the unit vectors kb.
            morph::Vector<T, 2> kb = b - k.path.back();
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
            //T Q = this->rcpt[0] / k.rcpt[0]; // reverse signalling
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
        // Test b, to see if it's near the border. Use winding number to see if it's
        // inside? Then, if outside, find out which edge it's nearest and apply that
        // force. Too complex. Instead, look at b's location. If x<0, then add component
        // to B[0]; if y<0 then add component to B[1], etc.
        if (b[0] < T{0}) {
            G = {0,0};
            I = {0,0};
            C = {0,0};
            B[0] = T{1};
        } else if (b[0] < r) {
            B[0] = T{1} * (T{1} - b[0]/r); // B[0] prop 1 - b/r
        } else if (b[0] > 1) {
            G = {0,0};
            I = {0,0};
            C = {0,0};
            B[0] = T{-1};
        } else if (b[0] > (1-r)) {
            B[0] = -(b[0] + r - T{1})/r; // B[0] prop (b+r-1)/r
        }

        if (b[1] < T{0}) {
            G = {0,0};
            I = {0,0};
            C = {0,0};
            B[1] = T{1};
        } else if (b[1] < r) {
            B[1] = T{1} - b[1]/r;
        } else if (b[1] > 1) {
            G = {0,0};
            I = {0,0};
            C = {0,0};
            B[1] = T{-1};
        } else if (b[1] > (1-r)) {
            B[1] = -(b[1] + r - T{1})/r; // B[1] prop (b+r-1)/r
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
                // Apply a slowdown based on distance toedge
                T toedge = distances.shortest();
                T speedlimit = T{0.1} + (T{2} / (T{1} + std::exp(-T{50} * std::abs(toedge)))) - T{1};
                //std::cout << "Not going in. To edge: " << toedge << " and speedlimit: " << speedlimit << std::endl;
                b += db*speedlimit;
                this->next = b;
            }

        } else {
            b += G * m[0] + C * m[1] + I * m[2] + B * m[3];
            this->next = b;
        }
    }
    // The location and all previous locations of this branch. Only populate for those
    // axons we're going to view. Also, only add to this when the axon moves an
    // appreciable distance, to keep it small. It can be dynamically managed in memory.
    morph::vVector<morph::Vector<T, 2>> path;
    // To allow use of vVector instead of deque, hold an iterator to the first member of
    // path that's 'in history', allowing path to become arbitrarily large in long
    // simulations (I accept the memory hit).
    //typename morph::vVector<morph::Vector<T, 2>>::iterator pathfront = path.begin();
    size_t pathfront = 0;
    // Pointer to path for current simulation step.
    //typename morph::vVector<morph::Vector<T, 2>>::iterator pathcur = path.begin();
    size_t pathcur = 0;

    // Place the next computed location for path in 'next' so that while computing, we
    // don't modify the numbers we're working from. After looping through all branches,
    // add this to path.
    morph::Vector<T, 2> next;
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
