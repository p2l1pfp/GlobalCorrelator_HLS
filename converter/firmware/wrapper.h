#ifndef WRAPPER_H
#define WRAPPER_H

#include <iostream>
#include <cmath>
#include "ap_int.h"
#include "ap_fixed.h"

/* #include "types.h" */
/* #include "convert_pt.h" */

#define PT_INV_TAB_SIZE 8
#define PT_INV_MAX_BITS 8
#define ETA_TAB_SIZE 8

typedef ap_uint<96> input_t;
typedef ap_uint<64> output_t;

void pf_input_track_conv_hw(input_t in, output_t& out);



enum TrackBitWidths { 
//MSB
    kValid              = 1,                                    // Valid bit

    kMVAExtraSize       = 6,                                    // Space for two specialized MVA selections
    kTrackMVASize       = 3,                                    // Width of track quality MVA
    kHitPatternSize     = 7,                                    // Width of the hit pattern for stubs
    kBendChi2Size       = 3,                                    // Width of the Bend-Chi2
    kChi2RZSize         = 4,                                    // Width of Chi2 for r-z
    kChi2RPhiSize       = 4,                                    // Width of Chi2 for r-phi
    kTrackQualitySize   = kMVAExtraSize + kTrackMVASize \
                        + kHitPatternSize + kBendChi2Size \
                        + kChi2RZSize + kChi2RPhiSize,          // Width of track quality

    kD0Size             = 13,                                   // Width of D0
    kD0MagSize          = 5,                                    // Width of D0 magnitude (signed)
    kZ0Size             = 12,                                   // Width of z-position
    kZ0MagSize          = 5,                                    // Width of z-position magnitude (signed)
    kEtaSize            = 16,                                   // Width of eta
    kEtaMagSize         = 3+1, // ch twiki seems wrong (sign)   // Width of eta magnitude (signed)
    kPhiSize            = 12,                                   // Width of phi
    kPhiMagSize         = 0,                                    // Width of phi magnitude (signed)
    kChargeSize         = 1,                                    // Charge of pT
    kChargeMagSize      = 1,                                    // Width of charge magnitude (signed)
    kPtSize             = 14,                                   // Width of pT
    kPtMagSize          = 0,                                    // Width of pT magnitude (unsigned)
//LSB
    kTrackParamsSize    = kD0Size + kZ0Size + kEtaSize \
                        + kPhiSize + kChargeSize + kPtSize,     // Width of track parameters

    kTrackWordSize      = kValid \
                        + kTrackQualitySize \
                        + kTrackParamsSize,                     // Width of the Track Finding Word in bits
};
// spare
typedef ap_uint<kValid>                                             valid_t;        // valid bit

// track quality types
typedef ap_uint<kTrackQualitySize>                                  trackQuality_t; // All track quality bits
typedef ap_uint<kMVAExtraSize>                                      extraMVA_t;      // Specialized MVA selection
typedef ap_uint<kTrackMVASize>                                      trackMVA_t;     // Track quality MVA
typedef ap_uint<kHitPatternSize>                                    hit_t;          // Hit mask bits
typedef ap_uint<kBendChi2Size>                                      bendChi2_t;     // Bend-Chi2
typedef ap_uint<kChi2RZSize>                                        chi2rz_t;       // Chi2 r-z
typedef ap_uint<kChi2RPhiSize>                                      chi2rphi_t;     // Chi2 r-phi

// track parameters types
typedef ap_fixed<kD0Size, kD0MagSize, AP_RND_CONV, AP_SAT>          tkd0_t;           // D0
typedef ap_fixed<kZ0Size, kZ0MagSize, AP_RND_CONV, AP_SAT>          tkz0_t;           // Track z
typedef ap_fixed<kEtaSize, kEtaMagSize, AP_RND_CONV, AP_SAT>        tanlam_t;         // Track eta
typedef ap_fixed<kPhiSize, kPhiMagSize, AP_RND_CONV, AP_SAT>        tkphi_t;          // Track phi
typedef bool                                                        q_t;              // Charge of track PT
typedef ap_ufixed<kPtSize, kPtMagSize, AP_RND_CONV, AP_SAT>         tkpt_t;           // Track PT
typedef ap_fixed<kPtSize+kChargeSize,0, AP_RND_CONV, AP_SAT>        rinv_t;           // 1/RT
typedef ap_ufixed<kPtSize,0, AP_RND_CONV, AP_SAT>                   urinv_t;           // 1/RT unsigned
typedef ap_ufixed<kEtaSize-1, kEtaMagSize-1, AP_RND_CONV, AP_SAT>   utanlam_t;         // Track eta

typedef ap_fixed<24,12, AP_RND_CONV, AP_SAT> bigfix_t; // helper type

//#include "../../submodules/GlobalCorrelator_HLS/firmware/data.h"
#define PF_PT_SCALE (4.0)
#define PF_ETAPHI_SCALE (4*180./3.1415)
#define PF_Z0_SCALE (20)
typedef ap_int<16> pt_t;
typedef ap_int<10> etaphi_t;
typedef ap_int<5>  vtx_t;
typedef ap_uint<3> particleid_t;
typedef ap_int<10> z0_t;

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
                    valid_t    valid   ); 

void unpack_L1T_track(ap_uint<kTrackWordSize> tk,
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
                      valid_t    &valid   );

void pack_pf_track(ap_uint<64> &tk,
                   pt_t     pf_pt   ,
                   pt_t     pf_pterr,
                   etaphi_t pf_eta  ,
                   etaphi_t pf_phi  ,
                   z0_t     pf_z0   ,
                   bool     pf_TightQuality);

void unpack_pf_track(ap_uint<64> tk,
                     pt_t     &pf_pt   ,
                     pt_t     &pf_pterr,
                     etaphi_t &pf_eta  ,
                     etaphi_t &pf_phi  ,
                     z0_t     &pf_z0   ,
                     bool     &pf_TightQuality);

template<class in_t, class out_t> void bit_copy(in_t in, out_t &out, int offset=0);

template<class pt_inv_T, class pt_T> void init_pt_inv_table(pt_T table_out[(1<<PT_INV_TAB_SIZE)]);
template<class pt_inv_T, class pt_T> void convert_pt(pt_inv_T inv, pt_T &pt);

template<class eta_T> void init_eta_table(eta_T table_out[(1<<ETA_TAB_SIZE)]);
template<class tanlam_T, class eta_T> void convert_eta(tanlam_T tanlam, eta_T &eta);

inline float tanlam_to_eta(float tanlam){return -log(tan((3.1415/2 - atan(tanlam))/2));}


#endif
