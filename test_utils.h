#include "firmware/data.h"

bool had_equals(const HadCaloObj &out_ref, const HadCaloObj &out, const char *what, int idx) ;
bool pf_equals(const PFChargedObj &out_ref, const PFChargedObj &out, const char *what, int idx) ;
bool pf_equals(const PFNeutralObj &out_ref, const PFNeutralObj &out, const char *what, int idx) ;

template<typename T> 
unsigned int count_nonzero(T objs[], unsigned int NMAX) {
    unsigned int ret = 0;
    for (unsigned int i = 0; i < NMAX; ++i) ret += (objs[i].hwPt > 0);
    return ret;
}
