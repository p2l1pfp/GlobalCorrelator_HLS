/*
HLS implementation of input conversion wrapper
*/
#include "tk_input_converter.h"
#include "tk_pack.h"
#include "tk_prop_tanlambda.h"
#include "tk_prop_phi.h"
#include "tk_pt_inversion.h"
#include "tk_resolution.h"
#include "tk_tanlambda_to_eta.h"
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


