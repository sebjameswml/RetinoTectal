#pragma once

#include "branch_base.h"

#include <vector>
#include <morph/vvec.h>
#include <morph/vec.h>
#include <morph/Random.h>
#include <morph/mathconst.h>
#include "tissue.h"

// A branch model inspired by Gebhardt et al 2012 (Balancing of ephrin/Eph forward and
// reverse signalling as the driving force of adaptive topographic mapping)
template<typename T, size_t N>
struct branch_geb : public branch_base<T,N>
{
    // The cone radius, hard coded
    void setr (T _r) { this->r = _r; }
    T getr() { return this->r; }
    T getr_c() { return this->r_c; }
    T getr_j() { return this->r_j; }
    T getr_i() { return this->r_i; }
    void setr_c(T _r) {}
    void setr_j(T _r) {}
    void setr_i(T _r) {}
protected:
    T r = T{0.05};
    T r_c = T{0.0};
    T r_j = T{0.0};
    T r_i = T{0.0};

public:
    // For random angle choice
    morph::RandUniform<T> rng;
    T s; // unused, dummy.
    morph::vec<T, N> minses; // unused, dummy.
    morph::vec<T, N> maxes; // unused, dummy.

    // No init for branch_geb
    virtual void init() {}

    // To achieve rotation of tectum wrt retina, as in biology then, if N==4:
    // receptor 0 interacts primarily with ligand 1
    // receptor 1 interacts primarily with ligand 0
    // receptor 2 interacts primarily with ligand 3 // opposite to receptor 0
    // receptor 3 interacts primarily with ligand 2 // opposite to receptor 1
    void compute_next (const std::vector<branch_geb<T, N>>& branches,
                       const guidingtissue<T, N>* source_tissue,
                       const guidingtissue<T, N>* tissue,
                       const morph::vec<T, 5>& m,
                       const morph::vec<T, 2*N>& rns)
    {
        // Sums
        morph::vec<T, N> LFRT; // reverse fibre-target signal (for N receptors, N ligands)
        morph::vec<T, N> LFRF; // reverse cis (fibre-self) signal
        morph::vec<T, N> LFRf; // reverse fibre - other fibres signal

        morph::vec<T, N> RFLT; // forward fibre-target signal
        morph::vec<T, N> RFLF; // forward cis (fibre-self) signal
        morph::vec<T, N> RFLf; // forward fibre - other fibres signal

        // Will need to run through other branches to determine the R_f/L_f
        morph::vvec<morph::vec<T, N>> Rf(tissue->posn.size());
        morph::vvec<morph::vec<T, N>> Lf(tissue->posn.size());
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
        morph::vvec<morph::vec<T,2>> candposns;
        candposns.push_back (this->current);
        // Choose a random direction to move a distance r (2r?) and evaluate that (G' in Gebhardt et al)
        T angle = this->rng.get() * morph::mathconst<T>::two_pi;
        morph::vec<T,2> randvec = { this->r * std::sin(angle), this->r * std::cos(angle) };
        candposns.push_back (this->current+randvec);

        //std::cout << "Current and candidate positions: " << candposns[0] << ", " << candposns[1] << std::endl;

        // Candidate probabilities
        morph::vvec<T> pdG (2, T{0});;

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
                morph::vec<T,2> pvec = tissue->posn[i] - candposns[j];
                T Gaussian = std::exp (-(a*(pvec[0]*pvec[0]) + a*(pvec[1]*pvec[1])));
                //std::cout << "Gaussian = " << Gaussian << std::endl;
                //if (pvec.length() <= this->r) { // Actually, in the Gebhardt work, r is defined by the Gaussian getting to 0.01
                if (Gaussian >= T{0.01}) {
                    // then posn[i] is under our cone!
                    morph::vec<T, N> LgndGauss = this->lgnd * Gaussian;
                    morph::vec<T, N> RcptGauss = this->rcpt * Gaussian;

                    // My scheme is that receptor 0 interacts with ligand 1 etc. to achieve
                    // rotation of tectum wrt retina.
                    morph::vec<T, N> LgndGaussR = LgndGauss;
                    morph::vec<T, N> RcptGaussR = RcptGauss;
                    morph::vec<T, N> tissueRcptR = tissue->rcpt[i];
                    morph::vec<T, N> tissueLgndR = tissue->lgnd[i];
                    morph::vec<T, N> RfR = Rf[i];
                    morph::vec<T, N> LfR = Lf[i];
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
            morph::vec<T, N> G_num = LFRT + LFRF + LFRf;
            morph::vec<T, N> G_denom = RFLT + RFLF + RFLf;
            morph::vec<T, N> G;
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
            pdG[j] = (T{1}/(sigma*morph::mathconst<T>::two_pi)) * std::exp ( - GFi * GFi / (2*sigma*sigma));
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
