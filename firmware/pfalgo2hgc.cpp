#include "pfalgo2hgc.h"

#include "pfalgo_common.icc"

void tk2calo_sumtk_hgc(const TkObj track[NTRACK], const bool isMu[NTRACK], const int tkerr2[NTRACK], const ap_uint<NCALO> calo_track_link_bit[NTRACK], pt_t sumtk[NCALO], int sumtkerr2[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t sum = 0;
        int sumerr = 0;
        for (int it = 0; it < NTRACK; ++it) {
            if (calo_track_link_bit[it][icalo] && !isMu[it]) { sum += track[it].hwPt; sumerr += tkerr2[it]; }
        }
        sumtk[icalo] = sum;
        sumtkerr2[icalo] = sumerr;
    }
}

void tk2calo_elealgo_hgc(const TkObj track[NTRACK], const HadCaloObj calo[NCALO], const ap_uint<NCALO> calo_track_link_bit[NTRACK], bool isEle[NTRACK]) {
    for (int it = 0; it < NTRACK; ++it) {
        bool ele = false;
        for (int icalo = 0; icalo < NCALO; ++icalo) {
            if (calo_track_link_bit[it][icalo] && calo[icalo].hwIsEM) ele = true;
        }
        isEle[it] = ele;
    }
} 

void tk2calo_tkalgo_hgc(const TkObj track[NTRACK], const bool isEle[NTRACK], const bool isMu[NTRACK], const ap_uint<NCALO> calo_track_link_bit[NTRACK], PFChargedObj pfout[NTRACK]) {
    const pt_t TKPT_MAX_LOOSE = PFALGO_TK_MAXINVPT_LOOSE; // 20 * PT_SCALE;
    const pt_t TKPT_MAX_TIGHT = PFALGO_TK_MAXINVPT_TIGHT; // 20 * PT_SCALE;
    for (int it = 0; it < NTRACK; ++it) {
        bool goodByPt = track[it].hwPt < (track[it].hwTightQuality ? TKPT_MAX_TIGHT : TKPT_MAX_LOOSE);
        bool good = isMu[it] || isEle[it] || goodByPt || calo_track_link_bit[it].or_reduce();
        if (good) {
            pfout[it].hwPt  = track[it].hwPt;
            pfout[it].hwEta = track[it].hwEta;
            pfout[it].hwPhi = track[it].hwPhi;
            pfout[it].hwId  = isMu[it] ? PID_Muon : ( isEle[it] ? PID_Electron : PID_Charged );
            pfout[it].hwZ0  = track[it].hwZ0;
        } else {
            pfout[it].hwPt  = 0;
            pfout[it].hwEta = 0;
            pfout[it].hwPhi = 0;
            pfout[it].hwId  = 0;
            pfout[it].hwZ0  = 0;
        }
    }
}


void tk2calo_caloalgo_hgc(const HadCaloObj calo[NCALO], const pt_t sumtk[NCALO], const int sumtkerr2[NCALO], PFNeutralObj pfout[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t calopt;
        if (sumtk[icalo] == 0) {
            calopt = calo[icalo].hwPt;
        } else {
            pt_t ptdiff = calo[icalo].hwPt - sumtk[icalo];
            if (ptdiff > 0 && (ptdiff*ptdiff) > (sumtkerr2[icalo] + (sumtkerr2[icalo] >> 1))) {
                calopt = ptdiff;
            } else {
                calopt = 0;
            }
        }
        pfout[icalo].hwPt  = calopt;
        pfout[icalo].hwEta = calopt ? calo[icalo].hwEta : etaphi_t(0);
        pfout[icalo].hwPhi = calopt ? calo[icalo].hwPhi : etaphi_t(0);
        pfout[icalo].hwId  = calopt ? (calo[icalo].hwIsEM ? PID_Photon : PID_Neutral) : 0;
    }
}



void pfalgo2hgc(const HadCaloObj calo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) {
    #pragma HLS ARRAY_PARTITION variable=calo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete    
    #pragma HLS ARRAY_PARTITION variable=outch complete
    #pragma HLS ARRAY_PARTITION variable=outne complete
    #pragma HLS ARRAY_PARTITION variable=outmu complete

    #pragma HLS pipeline II=2

    // ---------------------------------------------------------------
    // TK-MU Linking
    ap_uint<NMU> mu_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=mu_track_link_bit complete
    bool isMu[NTRACK];
    for (int it = 0; it < NTRACK; ++it) { isMu[it] = 0; }

    mutrk_link(mu, track, mu_track_link_bit);
    pfmualgo(mu, track, mu_track_link_bit, outmu, isMu);

    bool isEle[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=isEle complete

    // ---------------------------------------------------------------
    // TK-HAD Linking
    ap_uint<NCALO> calo_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=calo_track_link_bit complete
    tk2calo_link_drdpt(calo, track, calo_track_link_bit);
    //tk2calo_link_dronly(hadcalo_sub, track, calo_track_link_bit);

    int tkerr2[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=tkerr2 complete
    tk2calo_tkerr2(track, tkerr2);

    pt_t sumtk[NCALO]; int sumtkerr2[NCALO];
    #pragma HLS ARRAY_PARTITION variable=sumtk complete
    #pragma HLS ARRAY_PARTITION variable=sumtkerr2 complete

    tk2calo_elealgo_hgc(track, calo, calo_track_link_bit, isEle);

    tk2calo_tkalgo_hgc(track, isEle, isMu, calo_track_link_bit, outch);
    tk2calo_sumtk_hgc(track,  isMu, tkerr2, calo_track_link_bit, sumtk, sumtkerr2);

#if NCALO == NSELCALO
    tk2calo_caloalgo_hgc(calo, sumtk, sumtkerr2, outne);
#else
    PFNeutralObj outne_all[NCALO];
    #pragma HLS ARRAY_PARTITION variable=outne_all complete
    tk2calo_caloalgo_hgc(calo, sumtk, sumtkerr2, outne_all);
    ptsort_hwopt<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);
#endif

}

