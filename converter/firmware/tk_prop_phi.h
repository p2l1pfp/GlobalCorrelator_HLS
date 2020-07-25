//
#include "tk_input_converter.h"

template<class phi_T>
void init_dphi_table(phi_T table_out[(1<<DPHI_TAB_SIZE)]) {
    // index is a uint from 0 to 111...11=2^PT_INV_TAB_SIZE-1 that encodes 0.000 to 0.4999999
    // resulting value phi_T is a uint
    table_out[0] = 0;
    for (unsigned int i = 1; i < (1<<DPHI_TAB_SIZE); i++) {
        //float invpt = float(i)/(1<<DPHI_TAB_SIZE) * 0.5; // in 1/GeV
        //float invpt = float((2*i+1)/2.)/(1<<DPHI_TAB_SIZE) * 0.5; // in 1/GeV
        float invpt = float((2*i+4)/2.)/(1<<DPHI_TAB_SIZE) * 0.5; // in 1/GeV (works better somehow?)
        float rCurv = (1/invpt) * (129./2)/(0.735); // curv is pt * (looper radius / 2)/(looper pt)
        float x = (DETR)/(2*rCurv);
        float dPhi = atan(x / sqrt(1-x*x));
        table_out[i] = dPhi * PF_ETAPHI_SCALE;
        // overflow guard. shouldn't happen
        if( dPhi * PF_ETAPHI_SCALE >  (1<<(phi_T::width-1))-1)  table_out[i] = (1<<(phi_T::width-1))-1;
    }
    return;
}

template<class pt_inv_T, class phi_T> 
void convert_dphi_BRAM(pt_inv_T inv, phi_T &dphi){

    // Initialize the lookup tables
#ifdef __HLS_SYN__
    bool initialized = false;
    phi_T dphi_table[(1<<DPHI_TAB_SIZE)];
#else 
    static bool initialized = false;
    static phi_T dphi_table[(1<<DPHI_TAB_SIZE)];
#endif
    if (!initialized) {
        init_dphi_table<phi_T>(dphi_table);
        initialized = true;
    }

    if(inv<0) inv = -inv;
    urinv_t uinv = inv;

    // cutoffs at high and low pt
    if(uinv >= urinv_t(0.5)){
        dphi = 0.395 * PF_ETAPHI_SCALE; // low-pt: curv @ 2 GeV is 0.395
        return;
    } else if (uinv <= urinv_t(1./(1<<PT_INV_MAX_BITS))){
        dphi=0; // high-pt = straight track
        return;
    }

    ap_uint<DPHI_TAB_SIZE> index;
    const int offset = 1; // ignore the first bit since 0b0.01111.. = 0.499.. is largest value
    #pragma unroll
    for(int i=0; i<PT_INV_TAB_SIZE; i++){
        index[DPHI_TAB_SIZE-1-i] = uinv[urinv_t::width-1-i-offset]; //msb down to lowest
    }

    dphi = dphi_table[index];
}

template<class pt_inv_T, class phi_T> 
void convert_dphi_DSP(pt_inv_T inv, phi_T &dphi){
    // dPhi = r
    // (DETR)/(2*rCurv);
    //dphi = inv * DETR*(1/2.) * (0.735)/(129./2); // (735 MeV / 129/2 cm) pt -> radius conversion
    /* dphi = bigfix_t(DETR*(1/2.) * (0.735)/(129./2)) * inv; // (735 MeV / 129/2 cm) pt -> radius conversion */
    /* dphi = 0.4 * 0.5 * inv; //bigfix_t(DETR*(1/2.) * (0.735)/(129./2)) * inv; // (735 MeV / 129/2 cm) pt -> radius conversion */


    // kSynchrotron = (1.0/(0.3*3.8)); // in meter/GeV
    // rCurv = pt * kSynchrotron * 100 // in cm
    // dPhi = DETR / (2*rCurv)
    //      = DETR / (2 * pt * kSynchrotron * 100)
    //      = 1/pt * (DETR / (2 * kSynchrotron * 100))
    //dphi = inv * pt_inv_T(DETR / (2 * l1tk::kSynchrotron * 100.));
    //dphi = bigfix_t(inv) * bigfix_t(DETR / (2. * (1.0/(0.3*3.8)) * 100.) * PF_ETAPHI_SCALE);
    dphi = inv * bigfix_t(DETR / (2. * (1.0/(0.3*3.8)) * 100.) * PF_ETAPHI_SCALE);
    //std::cout << inv << " " <<  bigfix_t(DETR / (2. * (1.0/(0.3*3.8)) * 100.)) << " " << bigfix_t(inv) * bigfix_t(DETR / (2. * (1.0/(0.3*3.8)) * 100.)) << " " << dphi << std::endl;
    //std::cout << inv << " " <<  dphi << std::endl;
    //         dphi = inv * (DETR / (2 * l1tk::kSynchrotron * 100.));
}

template<class pt_inv_T, class phi_T> 
void convert_dphi(pt_inv_T inv, phi_T &dphi){
    //convert_dphi_BRAM(inv, dphi);
    convert_dphi_DSP(inv, dphi);
}
