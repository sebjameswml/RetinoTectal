/*
 * Retino-tectal system. This specialization inits all the RT projections.
 */

#include "rd_james_dncomp.h"

template <class Flt>
class RD_RetTec : public RD_James_dncomp<Flt>
{
public:
    RD_RetTec (void)
        : RD_James_dncomp<Flt>() {
    }

    virtual void allocate (void) {
        RD_James_dncomp<Flt>::allocate();
    }

    virtual void init (void) {
        RD_James_dncomp<Flt>::init();
#if 0
        for (unsigned int i = 0; i < this->N; ++i) {
            // Initialize according to config
        }
#endif
    }

}; // RD_RetTec
