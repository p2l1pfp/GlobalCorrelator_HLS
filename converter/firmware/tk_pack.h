#include "tk_input_converter.h"

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
    bool debug=false;
    if(debug){
        
        std::cout << "pack_pf_track" << std::endl;
        std::cout << "  pf_pt           " << pf_pt          .to_string(16) << std::endl;
        std::cout << "  pf_pterr        " << pf_pterr       .to_string(16) << std::endl;
        std::cout << "  pf_eta          " << pf_eta         .to_string(16) << std::endl;
        std::cout << "  pf_phi          " << pf_phi         .to_string(16) << std::endl;
        std::cout << "  pf_z0           " << pf_z0          .to_string(16) << std::endl;
        std::cout << "  pf_TightQuality " << pf_TightQuality << std::endl;
        
    }
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
