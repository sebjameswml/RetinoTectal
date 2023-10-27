#pragma once

#include "branch_base.h"

#include <limits>
#include <vector>
#include <bitset>
#include <morph/vec.h>
#include <morph/MathAlgo.h>
#include "tissue.h"

// A branch which implements the koulakov model.
// template param T The type for numbers (float or double)
// template param N The number of receptor-ligand species
template<typename T, size_t N>
struct branch_koulakov : public branch_base<T,N>
{
protected:
    // protected vars
public:
    // Do any init here
    virtual void init() {  }

    // Estimate the ligand gradient, given the true ligand gradient
    virtual morph::vec<T, 2*N> estimate_ligand_gradient (morph::vec<T,2*N>& true_lgnd_grad,
                                                         morph::vec<T,N>& true_lgnd_exp)
    {
        // This is the non-stochastic implementation...
        return true_lgnd_grad;
    }

    static constexpr bool debug_koulakov = true;
    // Compute the next position for this branch, using information from all other
    // branches and the parameters vector, m. Will also need location on target tissue,
    // to get its ephrin values. For that we have this->current.
    // rns are a set of random numbers to multiply the ligand gradient with.f
    int koulakov_swap (std::vector<branch_koulakov<T, N>>& branches, // Not const
                       const guidingtissue<T, N>* source_tissue,
                       const guidingtissue<T, N>* tissue,
                       const morph::vec<T, 5>& m,
                       const T gr, // Copy of gr from parent
                       morph::RandUniform<T, std::mt19937>* rng)
    {
        // Current location is named b
        morph::vec<T, 2> b = this->current;

        // Find the branches that are next door by evaluating distance and randomly
        // choose the axis to do the switching for.
        unsigned int axis = (rng->get() < T{0.5}) ? 1 : 0;

        int bmax = static_cast<int>(branches.size());
        int k = 0;
        bool got_neighbour = false;
        for (; k < bmax; ++k) {
            auto bk = branches[k].current - this->current;
            if constexpr (debug_koulakov) {
                if (bk.length() < 2*gr) {
                    std::cout << "Branch at "<< bk << ": bk[axis=" << axis <<"] = " << bk[axis] << "    bk[1-axis=" << (1-axis)<< "] = " << bk[1-axis] << "     gr=" << gr << std::endl;
                }
            }
            if ((bk[axis]-gr) < std::numeric_limits<T>::epsilon()
                && bk[1-axis] < std::numeric_limits<T>::epsilon()) {
                // This is the neighbour
                if constexpr (debug_koulakov) { std::cout << "Got neighbour\n"; }
                got_neighbour = true;
                break;
            }
        }
        if (!got_neighbour) {
            if constexpr (debug_koulakov) { std::cout << "No neighbour?\n"; }
            return -1;
        }

        // Do probabilities here...
        // ra[1]: receptor expression of this
        T ra1 = this->rcpt[axis]; // or 1-axis
        // ra[2]: receptor expression of branch[k]
        T ra2 = branches[k].rcpt[axis]; // or 1-axis
        // la[1] is Ligand expression at this->current
        unsigned int cur_t_idx = static_cast<unsigned int>(std::round(this->current[0]/gr) + std::round(this->current[1]/gr) * tissue->w);
        unsigned int br_t_idx = static_cast<unsigned int>(std::round(branches[k].current[0]/gr) + std::round(branches[k].current[1]/gr) * tissue->w);
        // or br_t_idx = cur_t_idx + axis ? 1 : tissue->w;

        T la1 = tissue->lgnd[axis][cur_t_idx];
        T la2 = tissue->lgnd[axis][br_t_idx];


        // la[2] is Ligand expression at branches[k].current
        // Re-use m_g for alpha in Koulakov's expression
        T p_ex = T{0.5} - m[0] * (ra1-ra2) * (la1-la2);
        if constexpr (debug_koulakov) {
            std::cout << "ra1/2: " << ra1 << "/"<< ra2 << " and la1/2: " << la1 << "/" << la2 << std::endl;
            std::cout << "p_ex = 0.5 - " << m[0] << " * (" << ra1 << "-" << ra2 << ") * (" << la1 << "-" << la2 << ")\n";
            std::cout << "p_ex = 0.5 - " << m[0] << " * (" << ra1 - ra2 << ") * (" << la1 - la2 << ")\n";
            std::cout << "p_ex = " << p_ex << std::endl;
        }

        // Get a random number
        if (rng->get() < p_ex) {
            // exchange! need to swap this and k
            this->current = branches[k].current;
            branches[k].current = b;

            // For debug
            branches[k].setr(0.02);
            this->setr(0.02);

            return k; // >-1 means exchange happened
        }

        return -1; // -1 means no exchange happened
    }
};
