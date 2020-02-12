/*
 * 2D Karbowski system with *divisive* normalization of a_i, deriving
 * from RD_James base class.
 */

#include "rd_james.h"

template <class Flt>
class RD_James_divnorm : public RD_James<Flt>
{
public:

    /*!
     * An N element vector holding the sum of a_i for each TC type.
     */
    alignas(vector<Flt>) vector<Flt> sum_a;

    /*!
     * An N element vector holding the initial sum of a_i for each TC type.
     */
    alignas(vector<Flt>) vector<Flt> sum_a_init;

    /*!
     * Simple constructor; no arguments. Just calls base constructor.
     */
    RD_James_divnorm (void)
        : RD_James<Flt>() {
    }

    virtual void allocate (void) {
        RD_James<Flt>::allocate();
        this->resize_vector_param (this->sum_a, this->N);
        this->resize_vector_param (this->sum_a_init, this->N);
    }

    virtual void init (void) {
        RD_James<Flt>::init();
        // Now compute sum of a and record this as sum_a_init.
        for (unsigned int i = 0; i < this->N; ++i) {
            this->sum_a_computation(i);
        }
        this->sum_a_init = this->sum_a;
#if 0
        // Print these on screen
        for (unsigned int ii = 0; ii < this->sum_a_init.size(); ++ii) {
            cout << "sum_a_init["<<ii<<"] = " << this->sum_a_init[ii] << endl;
        }
#endif
    }

    /*!
     * Computation methods
     */
    //@{

    /*!
     * A possibly normalization-function specific task to carry out
     * once after the sum of a has been computed.
     */
    virtual void sum_a_computation (const unsigned int _i) {
        // Compute the sum of a[i] across the sheet.
        this->sum_a[_i] = 0.0;
        Flt sum_tmp = 0.0;
#pragma omp parallel for reduction(+:sum_tmp)
        for (unsigned int h=0; h<this->nhex; ++h) {
            sum_tmp += this->a[_i][h];
        }
        this->sum_a[_i] = sum_tmp;
    }

    /*!
     * The normalization/transfer function.
     */
    virtual inline Flt transfer_a (const Flt& _a, const unsigned int _i) {
#ifdef NORMALIZE_TO_ONE
        // Divisive normalization to one
        Flt a_rtn = this->nhex * _a / this->sum_a[_i];
#else
        // Divisive norm with initial sum multiplier
        Flt a_rtn = this->sum_a_init[_i] * _a / this->sum_a[_i];
#endif
        // Prevent a from becoming negative, necessary only when competition is implemented:
        return (a_rtn < 0.0) ? 0.0 : a_rtn;
    }

}; // RD_James_norm
