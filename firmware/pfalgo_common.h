#ifndef FIRMWARE_PFALGO_COMMON_H
#define FIRMWARE_PFALGO_COMMON_H

#include "data.h"

inline int dr2_int(eta_t eta1, phi_t phi1, eta_t eta2, phi_t phi2) {
    ap_int<eta_t::width+1> deta = (eta1-eta2);
    ap_int<phi_t::width+1> dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}

#ifndef CMSSW_GIT_HASH
#define PFALGO_DR2MAX_TK_MU 2101
#endif

#endif
