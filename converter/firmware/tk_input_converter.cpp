/*
HLS implementation of input conversion wrapper
*/
#include "tk_input_converter.h"
#include "tk_pack.h"
#include "tk_prop_tanlambda.h"
//#include "tk_pt_inversion.h"
#include "tk_resolution.h"
#ifndef __SYNTHESIS__
#include <cstdio>
using std::cout;
using std::endl;
#endif

void pf_input_track_conv_hw(input_t in, output_t& out, numlink_t nlink){
    #pragma HLS pipeline II=1

    // unpack L1Track format
    rinv_t     rinv    ;
    tkphi_t    tkphi   ;
    tanlam_t   tanlam  ;
    tkz0_t     tkz0    ;
    tkd0_t     tkd0    ;
    chi2rphi_t chi2rphi;
    chi2rz_t   chi2rz  ;
    bendChi2_t bendChi2;
    hit_t      hit     ;
    trackMVA_t trackMVA;
    extraMVA_t extraMVA;
    valid_t    valid   ;

    unpack_L1T_track(in, rinv, tkphi, tanlam, tkz0, tkd0, chi2rphi, chi2rz, bendChi2, hit, trackMVA, extraMVA, valid);

    // selection
    if(!valid){
        out=0;
        return;
    }

    // track propagation to calo surface in eta
    tanlam_t tanlam_at_calo;
    propagate_tanlam(tkz0, tanlam, tanlam_at_calo);

    // track propagation to calo surface in phi
    rinv_t pt_inv = rinv;
    etaphi_t dphi;
    if (rinv<0) pt_inv = -rinv;
    convert_dphi(pt_inv, dphi);

    etaphi_t pf_phi_at_calo;
    if (rinv<0) pf_phi_at_calo = tkphi - dphi;
    else pf_phi_at_calo = tkphi + dphi;


    // converters
    pt_t conv_pt;
    //std::cout << "inverse pt " << rinv << std::endl;
    convert_pt(pt_inv, conv_pt);
    //std::cout << "        pt " << conv_pt << std::endl;

    etaphi_t pf_eta_at_calo;
    convert_eta(tanlam_at_calo, pf_eta_at_calo);

    //     // initialize conversion constants
    // #ifdef __HLS_SYN__
    //     bool initialized = false;
    //     pt_t phi_offsets[TT_NPHI_SECTORS];
    // #else 
    //     static bool initialized = false;
    //     static pt_t phi_offsets[TT_NPHI_SECTORS];
    // #endif
    //     if (!initialized) {
    //         for(int i=0;i<TT_NPHI_SECTORS;i++)
    //             phi_offsets[i] = (PF_ETAPHI_SCALE * M_PI / TT_NPHI_SECTORS)*(2*i-1);
    //         initialized = true;
    //     }
    //  etaphi_t pf_phi = bigfix_t(tkphi)*bigfix_t(PF_ETAPHI_SCALE) + phi_offsets[nlink];
    //etaphi_t pf_phi = tkphi;

    // z0_t pf_z0 = bigfix_t(tkz0)*bigfix_t(PF_Z0_SCALE); // scale=20
    // for now, copy z0 values without explicitly converting... must be replaced
    z0_t pf_z0 = zdet_t(tkz0) * zdet_t(PF_Z0_SCALE);

    // pack in PF format
    pt_t pf_pt = conv_pt;
    pt_t pf_pterr;
    reso_calo(pf_pt, pf_eta_at_calo, pf_pterr); // needs eta calo
    bool pf_TightQuality = 1;
    pack_pf_track(out, pf_pt, pf_pterr, pf_eta_at_calo, pf_phi_at_calo, pf_z0, pf_TightQuality);
    //std::cout << out.to_string(16) << std::endl;
}


template<class pt_T>
void init_pt_inv_table(pt_T table_out[(1<<PT_INV_TAB_SIZE)]) {
    // index is a uint from 0 to 111...11=2^PT_INV_TAB_SIZE-1 that encodes 0.000 to 0.4999999
    // resulting value pt_T is a uint from 0 to 2^16-1
    table_out[0] = (1<<PT_INV_MAX_BITS);
    for (unsigned int i = 1; i < (1<<PT_INV_TAB_SIZE); i++) {
        float invpt = float(i)/(1<<PT_INV_TAB_SIZE) * 0.5; // in 1/GeV
        table_out[i] = PF_PT_SCALE / invpt;
    }
    return;
}

template<class pt_inv_T, class pt_T>
void convert_pt(pt_inv_T inv, pt_T &pt){
    //pt_inv_T is ap_fixed<> (signed)
    //pt_T is ap_uint<> !

    // Initialize the lookup tables
#ifdef __HLS_SYN__
    bool initialized = false;
    pt_t inv_table[(1<<PT_INV_TAB_SIZE)];
#else 
    static bool initialized = false;
    static pt_t inv_table[(1<<PT_INV_TAB_SIZE)];
#endif
    if (!initialized) {
        init_pt_inv_table<pt_T>(inv_table);
        initialized = true;
    }

    // if(inv<0) inv = -inv; assume input is positive (upstream)
    urinv_t uinv = inv;

    // cutoffs at high and low pt
    if(uinv >= urinv_t(0.5)){
        pt = 2. * PF_PT_SCALE;
        return;
    } else if (uinv <= urinv_t(1./(1<<PT_INV_MAX_BITS))){
        pt=(1<<PT_INV_MAX_BITS) * PF_PT_SCALE;
        return;
    }

    ap_uint<PT_INV_TAB_SIZE> index;
    const int offset = 1; // ignore the first bit since 0b0.01111.. = 0.499.. is largest value
    #pragma unroll
    for(int i=0; i<PT_INV_TAB_SIZE; i++){
        index[PT_INV_TAB_SIZE-1-i] = uinv[urinv_t::width-1-i-offset]; //msb down to lowest
    }

    pt = inv_table[index];
}


template<class eta_T>
void init_eta_table(eta_T table_out[(1<<ETA_TAB_SIZE)]) {
    // index is a uint from 0 to 111...11=2^ETA_TAB_SIZE-1 that encodes 0.000 to 8
    // resulting value eta_T is a uint

    for (unsigned int i = 0; i < (1<<ETA_TAB_SIZE); i++) {
        // eta =  -ln(tan((pi/2 - arctan(TANLAM))/2))
        float tanlam = float(i)/(1<<ETA_TAB_SIZE) * 8.;
        float eta = -log(tan((M_PI/2 - atan(tanlam))/2));
        table_out[i] = eta * PF_ETAPHI_SCALE;
        // phi in -511,512 can only hold eta up to 2.23. else saturate for now
        if (eta * PF_ETAPHI_SCALE > (1<<(eta_T::width-1))-1) table_out[i] = (1<<(eta_T::width-1))-1;
    }
    return;
}

template<class tanlam_T, class eta_T>
void convert_eta(tanlam_T tanlam, eta_T &eta){
    // tanlam_T is ap_fixed<16,3>
    // eta_T is ap_int<10>

    // Initialize the lookup tables
#ifdef __HLS_SYN__
    bool initialized = false;
    etaphi_t eta_table[(1<<ETA_TAB_SIZE)];
#else 
    static bool initialized = false;
    static etaphi_t eta_table[(1<<ETA_TAB_SIZE)];
#endif
    if (!initialized) {
        init_eta_table<eta_T>(eta_table);
        initialized = true;
    }
    bool flip = false;
    if(tanlam<0){
        tanlam = -tanlam;
        flip=true;
    }
    utanlam_t utanlam = tanlam;

    ap_uint<ETA_TAB_SIZE> index;
    #pragma unroll
    for(int i=0; i<ETA_TAB_SIZE; i++){
        index[ETA_TAB_SIZE-1-i] = utanlam[utanlam_t::width-1-i]; //msb down to lowest
    }

    eta = eta_table[index];
    if(flip) eta = -eta;
}


template<class phi_T>
void init_dphi_table(phi_T table_out[(1<<DPHI_TAB_SIZE)]) {
    // index is a uint from 0 to 111...11=2^PT_INV_TAB_SIZE-1 that encodes 0.000 to 0.4999999
    // resulting value phi_T is a uint
    table_out[0] = 0;
    for (unsigned int i = 1; i < (1<<DPHI_TAB_SIZE); i++) {
        float invpt = float(i)/(1<<DPHI_TAB_SIZE) * 0.5; // in 1/GeV
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
void convert_dphi(pt_inv_T inv, phi_T &dphi){

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
