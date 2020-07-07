/*
HLS implementation of input conversion wrapper
*/
#include "wrapper.h"
#ifndef __SYNTHESIS__
#include <cstdio>
using std::cout;
using std::endl;
#endif

// #include "../../submodules/GTT_MET_HLS/eta/src/eta_top.cc"
// #include "../../submodules/GTT_MET_HLS/p_T/src/p_T_top.cc"

void pf_input_track_conv_hw(input_t in, output_t& out){
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

    // converters
    pt_t conv_pt;
    convert_pt(rinv,conv_pt);

    etaphi_t conv_eta;
    convert_eta(tanlam,conv_eta);

    // in_pt_t conv_in_pt = rinv;
    // out_pt_t conv_out_pt;
    // p_T(conv_in_pt, &conv_out_pt);

    // in_eta_t conv_in_eta = tanlam;
    // out_eta_t conv_out_eta;
    // bool parity=1;
    // bool erase;
    // eta(conv_in_eta, &conv_out_eta, parity, &erase);
    // tanlam_t tan_lambda = (M_PI/2.)-(2.*atan(exp(-1.*track_in[i].floatEta())));

    // pack in PF format
    pt_t pf_pt = conv_pt;
    pt_t pf_pterr = conv_pt; // TODO
    etaphi_t pf_eta = conv_eta;
    etaphi_t pf_phi = bigfix_t(tkphi)*bigfix_t(PF_ETAPHI_SCALE);
    z0_t pf_z0 = bigfix_t(tkz0)*bigfix_t(PF_Z0_SCALE);
    bool pf_TightQuality = 1;
    pack_pf_track(out, pf_pt, pf_pterr, pf_eta, pf_phi, pf_z0, pf_TightQuality);
}

void pack_L1T_track(ap_uint<kTrackWordSize> &tk,
                    rinv_t     rinv    ,
                    tkphi_t    tkphi   ,
                    tanlam_t   tanlam  ,
                    tkz0_t     tkz0    ,
                    tkd0_t     tkd0    ,
                    chi2rphi_t chi2rphi,
                    chi2rz_t   chi2rz  ,
                    bendChi2_t bendChi2,
                    hit_t      hit     ,
                    trackMVA_t trackMVA,
                    extraMVA_t extraMVA,
                    valid_t    valid   ){

    ap_uint<rinv_t    ::width> word_rinv      = rinv    .range(rinv_t    ::width-1,0);
    ap_uint<tkphi_t   ::width> word_tkphi     = tkphi   .range(tkphi_t   ::width-1,0);
    ap_uint<tanlam_t  ::width> word_tanlam    = tanlam  .range(tanlam_t  ::width-1,0);
    ap_uint<tkz0_t    ::width> word_tkz0      = tkz0    .range(tkz0_t    ::width-1,0);
    ap_uint<tkd0_t    ::width> word_tkd0      = tkd0    .range(tkd0_t    ::width-1,0);
    ap_uint<chi2rphi_t::width> word_chi2rphi  = chi2rphi.range(chi2rphi_t::width-1,0);
    ap_uint<chi2rz_t  ::width> word_chi2rz    = chi2rz  .range(chi2rz_t  ::width-1,0);
    ap_uint<bendChi2_t::width> word_bendChi2  = bendChi2.range(bendChi2_t::width-1,0);
    ap_uint<hit_t     ::width> word_hit       = hit     .range(hit_t     ::width-1,0);
    ap_uint<trackMVA_t::width> word_trackMVA  = trackMVA.range(trackMVA_t::width-1,0);
    ap_uint<extraMVA_t::width> word_extraMVA  = extraMVA.range(extraMVA_t::width-1,0);
    ap_uint<valid_t   ::width> word_valid     = valid   .range(valid_t   ::width-1,0);

    tk = (word_valid, word_extraMVA, word_trackMVA, word_hit, word_bendChi2, word_chi2rz, 
          word_chi2rphi, word_tkd0, word_tkz0, word_tanlam, word_tkphi, word_rinv);    
}

void unpack_L1T_track(ap_uint<kTrackWordSize> in,
                    rinv_t     &rinv    ,
                    tkphi_t    &tkphi   ,
                    tanlam_t   &tanlam  ,
                    tkz0_t     &tkz0    ,
                    tkd0_t     &tkd0    ,
                    chi2rphi_t &chi2rphi,
                    chi2rz_t   &chi2rz  ,
                    bendChi2_t &bendChi2,
                    hit_t      &hit     ,
                    trackMVA_t &trackMVA,
                    extraMVA_t &extraMVA,
                    valid_t    &valid   ){
    unsigned int lo = 0;
    unsigned int len = 0;
    len=rinv_t    ::width; bit_copy(in, rinv    , lo); lo += len;
    len=tkphi_t   ::width; bit_copy(in, tkphi   , lo); lo += len;
    len=tanlam_t  ::width; bit_copy(in, tanlam  , lo); lo += len;
    len=tkz0_t    ::width; bit_copy(in, tkz0    , lo); lo += len;
    len=tkd0_t    ::width; bit_copy(in, tkd0    , lo); lo += len;
    len=chi2rphi_t::width; bit_copy(in, chi2rphi, lo); lo += len;
    len=chi2rz_t  ::width; bit_copy(in, chi2rz  , lo); lo += len;
    len=bendChi2_t::width; bit_copy(in, bendChi2, lo); lo += len;
    len=hit_t     ::width; bit_copy(in, hit     , lo); lo += len;
    len=trackMVA_t::width; bit_copy(in, trackMVA, lo); lo += len;
    len=extraMVA_t::width; bit_copy(in, extraMVA, lo); lo += len;
    len=valid_t   ::width; bit_copy(in, valid   , lo); lo += len;
}

void pack_pf_track(ap_uint<64> &tk,
                   pt_t     pf_pt   ,
                   pt_t     pf_pterr,
                   etaphi_t pf_eta  ,
                   etaphi_t pf_phi  ,
                   z0_t     pf_z0   ,
                   bool     pf_TightQuality){
    tk = (pf_TightQuality, pf_z0, pf_phi, pf_eta, pf_pterr, pf_pt);    
}

void unpack_pf_track(ap_uint<64> in,
                   pt_t     &pf_pt   ,
                   pt_t     &pf_pterr,
                   etaphi_t &pf_eta  ,
                   etaphi_t &pf_phi  ,
                   z0_t     &pf_z0   ,
                   bool     &pf_TightQuality){
    unsigned int lo = 0;
    unsigned int len = 0;
    len=pt_t    ::width; bit_copy(in, pf_pt          , lo); lo += len;
    len=pt_t    ::width; bit_copy(in, pf_pterr       , lo); lo += len;
    len=etaphi_t::width; bit_copy(in, pf_eta         , lo); lo += len;
    len=etaphi_t::width; bit_copy(in, pf_phi         , lo); lo += len;
    len=z0_t    ::width; bit_copy(in, pf_z0          , lo); lo += len;
    pf_TightQuality = in[lo];
}

template<class in_t, class out_t> 
void bit_copy(in_t in, out_t &out, int offset){
    for(int i  =out_t::width-1; i>=0; i--){
        out[i] = in[i+offset];
    }    
}

template<class pt_inv_T, class pt_T>
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
        init_pt_inv_table<pt_inv_T,pt_T>(inv_table);
        initialized = true;
    }

    if(inv<0) inv = -inv;
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
