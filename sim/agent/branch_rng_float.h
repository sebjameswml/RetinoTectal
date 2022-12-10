/*
 * A singleton RNG function used by branch.h to add a small stochastic movement to each
 * branch on each timestep.
 */
#pragma once
#include <morph/Random.h>
#include <morph/vec.h>
// brng: "branch random number generator"
class brng
{
private:
    brng() { this->rng = new morph::RandUniform<float>(-1.0f, 1.0f); }
    ~brng() { delete this->rng; }
    static brng* pInstance_brng;
public:
    static brng* i() // The instance public function.
    {
        if (brng::pInstance_brng == 0) { brng::pInstance_brng = new brng; }
        return brng::pInstance_brng;
    }
    float get() { return this->rng->get(); }

    template<size_t arraysz = 2>
    void get(morph::vec<float,arraysz>& ar) { return this->rng->get(ar); }

    morph::RandUniform<float>* rng;
};
// Globally initialise brng instance pointer to NULL. Wherever you know what T is
brng* brng::pInstance_brng = 0;
