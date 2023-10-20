#pragma once

#include "branch_base.h"

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

    // Compute the next position for this branch, using information from all other
    // branches and the parameters vector, m. Will also need location on target tissue,
    // to get its ephrin values.
    // rns are a set of random numbers to multiply the ligand gradient with.
    void compute_next (const std::vector<branch_koulakov<T, N>>& branches,
                       const guidingtissue<T, N>* source_tissue,
                       const guidingtissue<T, N>* tissue,
                       const morph::vec<T, 5>& m,
                       const morph::vec<T, 2*N>& rns)
    {
        // Current location is named b
        morph::vec<T, 2> b = this->current;

        // The change in location computed from the model
        morph::vec<T, 2> db = {std::numeric_limits<T>::epsilon(), 0.001f}; // placeholder

        // Finally add b and db to get next
        this->next = b + db;
    }
};
