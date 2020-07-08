#ifndef LONETRACK_TYPES_H
#define LONETRACK_TYPES_H
// From https://gitlab.cern.ch/GTT/common/-/blob/master/DataFormats/interface/Track.h
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/HybridDataFormat#Fitted_Tracks_written_by_KalmanF

namespace l1tk {
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
        kTrackParamsSize    = kD0Size + kZ0Size + kEtaSize              \
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

    const double kSynchrotron = (1.0/(0.3*3.8));
    const double kRmax = (128.*100.0*kSynchrotron); //128 is max pt
    const double kRmin = (2.*100.0*kSynchrotron); //2 is min pt
    const float rzphiChi2Bins[16] = {0, 0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 40, 100, 200, 500, 1000, 3000};
}
#endif
