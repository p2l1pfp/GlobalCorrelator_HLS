#include "pfalgo3.h"

#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
bool gdebug_;
void pfalgo3_set_debug(bool debug) { gdebug_ = debug; }
#else
void pfalgo3_set_debug(bool debug) { }
#endif

#include "pfalgo_common.icc"

template<int DR2MAX>
void tk2em_drvals(const EmCaloObj calo[NEMCALO], const TkObj track[NTRACK], const bool isMu[NTRACK], tk2em_dr_t calo_track_drval[NTRACK][NEMCALO]) {
    const tk2em_dr_t eDR2MAX = DR2MAX;
    for (int it = 0; it < NTRACK; ++it) {
        for (int icalo = 0; icalo < NEMCALO; ++icalo) {
            if (isMu[it] || track[it].hwPt == 0 || calo[icalo].hwPt == 0) {calo_track_drval[it][icalo] = eDR2MAX; } // set to DR max if the track is a muon or null
            else { calo_track_drval[it][icalo] = dr2_int_cap(track[it].hwEta, track[it].hwPhi, calo[icalo].hwEta, calo[icalo].hwPhi, eDR2MAX); } 
        }
    }
}


template<int DR2MAX>
void em2calo_drvals(const HadCaloObj hadcalo[NCALO], const EmCaloObj emcalo[NEMCALO], em2calo_dr_t hadcalo_emcalo_drval[NEMCALO][NCALO]) {
    const em2calo_dr_t eDR2MAX = DR2MAX;
    for (int it = 0; it < NEMCALO; ++it) {
        pt_t hadcaloEmPtMin = (emcalo[it].hwPt >> 1);
        for (int icalo = 0; icalo < NCALO; ++icalo) {
            if (hadcalo[icalo].hwEmPt > hadcaloEmPtMin) {
                hadcalo_emcalo_drval[it][icalo] = dr2_int_cap(emcalo[it].hwEta, emcalo[it].hwPhi, hadcalo[icalo].hwEta, hadcalo[icalo].hwPhi, eDR2MAX);
            } else {
                hadcalo_emcalo_drval[it][icalo] = eDR2MAX;
            }
        }
    }
}


void tk2em_link(const EmCaloObj calo[NEMCALO], const TkObj track[NTRACK], const bool isMu[NTRACK], ap_uint<NEMCALO> calo_track_link_bit[NTRACK]) {
    const int DR2MAX = PFALGO_DR2MAX_TK_EM;

    tk2em_dr_t drvals[NTRACK][NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0

    tk2em_drvals<DR2MAX>(calo, track, isMu, drvals);
    pick_closest<DR2MAX,NTRACK,NEMCALO,tk2em_dr_t>(drvals, calo_track_link_bit);
}

void em2calo_link(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], ap_uint<NCALO> em_calo_link_bit[NCALO]) {
    const int DR2MAX = PFALGO_DR2MAX_EM_CALO;
    em2calo_dr_t drvals[NEMCALO][NCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0

    em2calo_drvals<DR2MAX>(hadcalo, emcalo, drvals);
    pick_closest<DR2MAX,NEMCALO,NCALO,em2calo_dr_t>(drvals, em_calo_link_bit);
}

void tk2calo_sumtk(const TkObj track[NTRACK], const bool isEle[NTRACK], const bool isMu[NTRACK], const int tkerr2[NTRACK], const ap_uint<NCALO> calo_track_link_bit[NTRACK], pt_t sumtk[NCALO], int sumtkerr2[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t sum = 0;
        int sumerr = 0;
        for (int it = 0; it < NTRACK; ++it) {
            if (!isEle[it] && calo_track_link_bit[it][icalo] && !isMu[it]) { sum += track[it].hwPt; sumerr += tkerr2[it]; }
        }
        sumtk[icalo] = sum;
        sumtkerr2[icalo] = sumerr;
    }
}

void tk2em_sumtk(const TkObj track[NTRACK], const ap_uint<NEMCALO> calo_track_link_bit[NTRACK], pt_t sumtk[NEMCALO]) {
    for (int icalo = 0; icalo < NEMCALO; ++icalo) {
        pt_t sum = 0;
        for (int it = 0; it < NTRACK; ++it) {
            if (calo_track_link_bit[it][icalo]) { sum += track[it].hwPt; }
        }
        sumtk[icalo] = sum;
    }
}

void em2calo_sumem(const EmCaloObj emcalo[NEMCALO], const bool isEM[NEMCALO], const ap_uint<NCALO> em_had_link_bit[NTRACK], pt_t sumem[NEMCALO], bool keepcalo[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t sum = 0; bool keep = false;
        for (int iem = 0; iem < NEMCALO; ++iem) {
            if (em_had_link_bit[iem][icalo]) { 
                if (isEM[iem]) sum += emcalo[iem].hwPt; 
                else keep = true;
            }
        }
        sumem[icalo] = sum;
        keepcalo[icalo] = keep;
    }
}

void tk2calo_tkalgo(const TkObj track[NTRACK], const bool isEle[NTRACK], const bool isMu[NTRACK], const ap_uint<NCALO> calo_track_link_bit[NTRACK], PFChargedObj pfout[NTRACK]) {
    const pt_t TKPT_MAX_LOOSE = PFALGO_TK_MAXINVPT_LOOSE; // 20 * PT_SCALE;
    const pt_t TKPT_MAX_TIGHT = PFALGO_TK_MAXINVPT_TIGHT; // 20 * PT_SCALE;
    for (int it = 0; it < NTRACK; ++it) {
        bool goodByPt = track[it].hwPt < (track[it].hwTightQuality ? TKPT_MAX_TIGHT : TKPT_MAX_LOOSE);
        bool good = isMu[it] || isEle[it] || goodByPt || calo_track_link_bit[it].or_reduce();
        if (good) {
            pfout[it].hwPt  = track[it].hwPt;
            pfout[it].hwEta = track[it].hwEta;
            pfout[it].hwPhi = track[it].hwPhi;
            pfout[it].hwId  = isEle[it] ? PID_Electron : ( isMu[it] ? PID_Muon : PID_Charged );
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


void tk2calo_caloalgo(const HadCaloObj calo[NCALO], const pt_t sumtk[NCALO], const int sumtkerr2[NCALO], PFNeutralObj pfout[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t calopt;
        if (sumtk[icalo] == 0) {
            calopt = calo[icalo].hwPt;
            #ifndef __SYNTHESIS__
            if (gdebug_ && calo[icalo].hwPt) printf("HW hadcalo %2d pt %7d sumtk %7d                                -->  calopt %7d \n", icalo, int(calo[icalo].hwPt), 
                   int(sumtk[icalo]), int(calopt));
            #endif
        } else {
            pt_t ptdiff = calo[icalo].hwPt - sumtk[icalo];
            if (ptdiff > 0 && (ptdiff*ptdiff > sumtkerr2[icalo])) {
                calopt = ptdiff;
            } else {
                calopt = 0;
            }
            #ifndef __SYNTHESIS__
            if (gdebug_ && calo[icalo].hwPt) printf("HW hadcalo %2d pt %7d sumtk %7d +- %7.1f --> ptdiff %7d  -->  calopt %7d \n", icalo, int(calo[icalo].hwPt), 
                   int(sumtk[icalo]), std::sqrt(float(sumtkerr2[icalo])), int(ptdiff), int(calopt));
            #endif
        }
        pfout[icalo].hwPt  = calopt;
        pfout[icalo].hwEta = calopt ? calo[icalo].hwEta : eta_t(0);
        pfout[icalo].hwPhi = calopt ? calo[icalo].hwPhi : phi_t(0);
        pfout[icalo].hwId  = calopt ? PID_Neutral : 0;
    }
}

void tk2em_emalgo(const EmCaloObj calo[NEMCALO], const pt_t sumtk[NEMCALO], bool isEM[NEMCALO], pt_t photonPt[NEMCALO]) {
    for (int icalo = 0; icalo < NEMCALO; ++icalo) {
        if (sumtk[icalo] == 0) {
            isEM[icalo] = true;
            photonPt[icalo] = calo[icalo].hwPt;
        } else {
            pt_t ptdiff = calo[icalo].hwPt - sumtk[icalo];
            int ptdiff2 = ptdiff*ptdiff;
            int sigma2 = (calo[icalo].hwPtErr*calo[icalo].hwPtErr);
            int sigma2Lo = (sigma2 << 2), sigma2Hi = sigma2; // + (sigma2 >> 1);
            if ((ptdiff >= 0 && ptdiff2 <= sigma2Hi) || (ptdiff < 0 && ptdiff2 < sigma2Lo)) {
                photonPt[icalo] = 0;    
                isEM[icalo] = true;
            } else if (ptdiff > 0) {
                photonPt[icalo] = ptdiff;    
                isEM[icalo] = true;
            } else {
                photonPt[icalo] = 0;    
                isEM[icalo] = false;
            }
        }
    }
}
void tk2em_elealgo(const ap_uint<NEMCALO> em_track_link_bit[NTRACK], const bool isEM[NEMCALO], bool isEle[NTRACK]) {
    for (int it = 0; it < NTRACK; ++it) {
        bool ele = false;
        for (int icalo = 0; icalo < NEMCALO; ++icalo) { 
            if (isEM[icalo] && em_track_link_bit[it][icalo]) ele = true;
        }
        isEle[it] = ele;
        #ifndef __SYNTHESIS__
        if (gdebug_ && ele) printf("HW track %2d promoted to an electron\n", it);
        #endif
    }
}
void tk2em_photons(const EmCaloObj calo[NEMCALO], const pt_t photonPt[NEMCALO], PFNeutralObj pfout[NSELCALO]) {
    for (int icalo = 0; icalo < NEMCALO; ++icalo) {
        pfout[icalo].hwPt  = photonPt[icalo];
        pfout[icalo].hwEta = photonPt[icalo] ? calo[icalo].hwEta : eta_t(0);
        pfout[icalo].hwPhi = photonPt[icalo] ? calo[icalo].hwPhi : phi_t(0);
        pfout[icalo].hwId  = photonPt[icalo] ? PID_Photon : 0;
        #ifndef __SYNTHESIS__
        if (gdebug_ && photonPt[icalo]) printf("HW emcalo %2d pt %7d promoted to a photon with pt %7d\n", icalo, int(calo[icalo].hwPt), int(photonPt[icalo]));
        #endif
    }
}

void em2calo_sub(const HadCaloObj calo[NCALO], const pt_t sumem[NCALO], const bool keepcalo[NCALO], HadCaloObj calo_out[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t ptsub = calo[icalo].hwPt   - sumem[icalo];
        pt_t emsub = calo[icalo].hwEmPt - sumem[icalo];
        if ((ptsub <= (calo[icalo].hwPt >> 4)) || 
                (calo[icalo].hwIsEM && (emsub <= (calo[icalo].hwEmPt>>3)) && !keepcalo[icalo])) {
#ifndef __SYNTHESIS__
            if (gdebug_ && calo[icalo].hwPt) printf("HW hadcalo %2d pt %7d empt %7d sumem %7d keepcalo %1d  --> discarded\n", icalo, int(calo[icalo].hwPt), int(calo[icalo].hwEmPt), int(sumem[icalo]), int(keepcalo[icalo]));
#endif
            calo_out[icalo].hwPt   = 0;
            calo_out[icalo].hwEmPt = 0;
            calo_out[icalo].hwEta  = 0;
            calo_out[icalo].hwPhi  = 0;
            calo_out[icalo].hwIsEM = 0;
        } else {
#ifndef __SYNTHESIS__
            if (gdebug_ && calo[icalo].hwPt) printf("HW hadcalo %2d pt %7d empt %7d sumem %7d keepcalo %1d  --> kept with pt %7d empt %7d\n", icalo, int(calo[icalo].hwPt), int(calo[icalo].hwEmPt), int(sumem[icalo]), int(keepcalo[icalo]), int(ptsub), int((emsub > 0 ? emsub : pt_t(0))));
#endif
            calo_out[icalo].hwPt   = ptsub;
            calo_out[icalo].hwEmPt = (emsub > 0 ? emsub : pt_t(0));
            calo_out[icalo].hwEta  = calo[icalo].hwEta;
            calo_out[icalo].hwPhi  = calo[icalo].hwPhi;
            calo_out[icalo].hwIsEM = calo[icalo].hwIsEM;
        }
    }
}




void pfalgo3(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) {
    #pragma HLS ARRAY_PARTITION variable=emcalo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete    
    #pragma HLS ARRAY_PARTITION variable=outch complete
    #pragma HLS ARRAY_PARTITION variable=outpho complete
    #pragma HLS ARRAY_PARTITION variable=outne complete
    #pragma HLS ARRAY_PARTITION variable=outmu complete

#ifdef HLS_pipeline_II
 #if HLS_pipeline_II == 1
    #pragma HLS pipeline II=1
 #elif HLS_pipeline_II == 2
    #pragma HLS pipeline II=2
 #elif HLS_pipeline_II == 3
    #pragma HLS pipeline II=3
 #elif HLS_pipeline_II == 4
    #pragma HLS pipeline II=4
 #elif HLS_pipeline_II == 6
    #pragma HLS pipeline II=6
 #endif
#else
    #pragma HLS pipeline II=2
#endif

    // ---------------------------------------------------------------
    // TK-MU Linking
    ap_uint<NMU> mu_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=mu_track_link_bit complete
    bool isMu[NTRACK];
    for (int it = 0; it < NTRACK; ++it) { isMu[it] = 0; }

    mutrk_link(mu, track, mu_track_link_bit);
    pfmualgo(mu, track, mu_track_link_bit, outmu, isMu);
    #ifndef __SYNTHESIS__
    for (int it = 0; it < NTRACK; ++it) { for (int imu = 0; imu < NMU; ++imu) {
        if (gdebug_ && isMu[it] && mu_track_link_bit[it][imu]) printf("HW track %2d is linked to muon %d\n", it, imu);
    } }
    #endif

    // ---------------------------------------------------------------
    // TK-EM Linking
    ap_uint<NEMCALO> em_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=em_track_link_bit complete

    tk2em_link(emcalo, track, isMu, em_track_link_bit);
    #ifndef __SYNTHESIS__
    for (int it = 0; it < NTRACK; ++it) { for (int ic = 0; ic < NEMCALO; ++ic) {
        if (gdebug_ && track[it].hwPt && em_track_link_bit[it][ic]) printf("HW track %2d pt %7d is linked to em calo %2d pt %7d\n", it, int(track[it].hwPt), ic, int(emcalo[ic].hwPt));
    } }
    #endif
    
    pt_t sumtk2em[NEMCALO]; 
    #pragma HLS ARRAY_PARTITION variable=sumtk2em complete

    pt_t photonPt[NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=photonPt complete

    bool isEM[NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=isEM complete

    tk2em_sumtk(track, em_track_link_bit, sumtk2em);
    tk2em_emalgo(emcalo, sumtk2em, isEM, photonPt);
    tk2em_photons(emcalo, photonPt, outpho);

    bool isEle[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=isEle complete

   tk2em_elealgo(em_track_link_bit, isEM, isEle);

    ap_uint<NCALO> em_calo_link_bit[NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=em_calo_link_bit complete
    em2calo_link(emcalo, hadcalo, em_calo_link_bit);
    #ifndef __SYNTHESIS__
    for (int ie = 0; ie < NEMCALO; ++ie) { for (int ic = 0; ic < NCALO; ++ic) {
        if (gdebug_ && emcalo[ie].hwPt && em_calo_link_bit[ie][ic]) printf("HW em calo %2d pt %7d is linked to had calo %2d pt %7d empt %7d\n", ie, int(emcalo[ie].hwPt), ic, int(hadcalo[ic].hwPt), int(hadcalo[ic].hwEmPt));
    } }
    #endif

    bool keepcalo[NCALO];
    pt_t sumem[NCALO]; 
    #pragma HLS ARRAY_PARTITION variable=sumem complete
    em2calo_sumem(emcalo, isEM, em_calo_link_bit, sumem, keepcalo);

    HadCaloObj hadcalo_sub[NCALO];
    #pragma HLS ARRAY_PARTITION variable=hadcalo_sub complete
    em2calo_sub(hadcalo, sumem, keepcalo, hadcalo_sub);

#if 0    
    HadCaloObj hadcalo_sub_out[NCALO];
    #pragma HLS ARRAY_PARTITION variable=hadcalo_sub_out complete
    buffer_ff<HadCaloObj,NCALO>(hadcalo_sub,hadcalo_sub_out);
#endif

    // ---------------------------------------------------------------
    // TK-HAD Linking
    ap_uint<NCALO> calo_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=calo_track_link_bit complete

    tk2calo_link_drdpt(hadcalo_sub, track, calo_track_link_bit);
    #ifndef __SYNTHESIS__
    for (int it = 0; it < NTRACK; ++it) { for (int ic = 0; ic < NCALO; ++ic) {
        if (gdebug_ && track[it].hwPt && calo_track_link_bit[it][ic]) printf("HW track %2d pt %7d is linked to had calo %2d pt %7d\n", it, int(track[it].hwPt), ic, int(hadcalo_sub[ic].hwPt));
    } }
    #endif

    int tkerr2[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=tkerr2 complete
    tk2calo_tkerr2(track, tkerr2);

    pt_t sumtk[NCALO]; int sumtkerr2[NCALO];
    #pragma HLS ARRAY_PARTITION variable=sumtk complete
    #pragma HLS ARRAY_PARTITION variable=sumtkerr2 complete

    tk2calo_tkalgo(track, isEle, isMu, calo_track_link_bit, outch);
    tk2calo_sumtk(track, isEle, isMu, tkerr2, calo_track_link_bit, sumtk, sumtkerr2);

    PFNeutralObj outne_all[NCALO];
    #pragma HLS ARRAY_PARTITION variable=outne_all complete
    tk2calo_caloalgo(hadcalo_sub, sumtk, sumtkerr2, outne_all);

    ptsort_hwopt<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);
}

#if defined(PACKING_DATA_SIZE) && defined(PACKING_NCHANN)
#include "l1pf_encoding.h"

void packed_pfalgo3(const ap_uint<PACKING_DATA_SIZE> input[PACKING_NCHANN], ap_uint<PACKING_DATA_SIZE> output[PACKING_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=output complete
    #pragma HLS pipeline II=2

    EmCaloObj emcalo[NEMCALO]; HadCaloObj hadcalo[NCALO]; TkObj track[NTRACK]; MuObj mu[NMU]; 
    PFChargedObj outch[NTRACK]; PFNeutralObj outpho[NPHOTON], outne[NSELCALO]; PFChargedObj outmu[NMU];
    #pragma HLS ARRAY_PARTITION variable=emcalo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete    
    #pragma HLS ARRAY_PARTITION variable=outch complete
    #pragma HLS ARRAY_PARTITION variable=outpho complete
    #pragma HLS ARRAY_PARTITION variable=outne complete
    #pragma HLS ARRAY_PARTITION variable=outmu complete

    pfalgo3_unpack_in(input, emcalo, hadcalo, track, mu);
    pfalgo3(emcalo, hadcalo, track, mu, outch, outpho, outne, outmu);
    pfalgo3_pack_out(outch, outpho, outne, outmu, output);
}

void pfalgo3_pack_in(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], ap_uint<PACKING_DATA_SIZE> input[PACKING_NCHANN]) {
    assert(NTRACK+NEMCALO+NCALO+NMU <= PACKING_NCHANN);
    const int EMCALO_OFFS = 0, HADCALO_OFFS = NEMCALO, TRACK_OFFS = HADCALO_OFFS + NCALO, MU_OFFS = TRACK_OFFS + NTRACK;
    l1pf_pattern_pack<NEMCALO, EMCALO_OFFS >(emcalo,  input);
    l1pf_pattern_pack<NCALO,   HADCALO_OFFS>(hadcalo, input);
    l1pf_pattern_pack<NTRACK,  TRACK_OFFS  >(track,   input);
    l1pf_pattern_pack<NMU,     MU_OFFS     >(mu,      input);
}

void pfalgo3_unpack_in(const ap_uint<PACKING_DATA_SIZE> input[PACKING_NCHANN], EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU]) {
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=emcalo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete    
    #pragma HLS inline
    assert(NTRACK+NPHOTON+NSELCALO+NMU <= PACKING_NCHANN);
    const int EMCALO_OFFS = 0, HADCALO_OFFS = NEMCALO, TRACK_OFFS = HADCALO_OFFS + NCALO, MU_OFFS = TRACK_OFFS + NTRACK;
    l1pf_pattern_unpack<NEMCALO, EMCALO_OFFS >(input, emcalo);
    l1pf_pattern_unpack<NCALO,   HADCALO_OFFS>(input, hadcalo);
    l1pf_pattern_unpack<NTRACK,  TRACK_OFFS  >(input, track);
    l1pf_pattern_unpack<NMU,     MU_OFFS     >(input, mu);
}

void pfalgo3_pack_out(const PFChargedObj outch[NTRACK], const PFNeutralObj outpho[NPHOTON], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU], ap_uint<PACKING_DATA_SIZE> output[PACKING_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=output complete
    #pragma HLS ARRAY_PARTITION variable=outch complete
    #pragma HLS ARRAY_PARTITION variable=outpho complete
    #pragma HLS ARRAY_PARTITION variable=outne complete
    #pragma HLS ARRAY_PARTITION variable=outmu complete
    #pragma HLS inline

    assert(NTRACK+NPHOTON+NSELCALO+NMU <= PACKING_NCHANN);
    const int PFCH_OFFS = 0, PFPH_OFFS = PFCH_OFFS + NTRACK, PFNE_OFFS = PFPH_OFFS + NPHOTON, PFMU_OFFS = PFNE_OFFS + NSELCALO;
    l1pf_pattern_pack<NTRACK,  PFCH_OFFS>(outch,  output);
    l1pf_pattern_pack<NPHOTON, PFPH_OFFS>(outpho, output);
    l1pf_pattern_pack<NSELCALO,PFNE_OFFS>(outne,  output);
    l1pf_pattern_pack<NMU,     PFMU_OFFS>(outmu,  output);
}

void pfalgo3_unpack_out(const ap_uint<PACKING_DATA_SIZE> output[PACKING_NCHANN], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) {
    assert(NTRACK+NPHOTON+NSELCALO+NMU <= PACKING_NCHANN);
    const int PFCH_OFFS = 0, PFPH_OFFS = PFCH_OFFS + NTRACK, PFNE_OFFS = PFPH_OFFS + NPHOTON, PFMU_OFFS = PFNE_OFFS + NSELCALO;
    l1pf_pattern_unpack<NTRACK,  PFCH_OFFS>(output, outch);
    l1pf_pattern_unpack<NPHOTON, PFPH_OFFS>(output, outpho);
    l1pf_pattern_unpack<NSELCALO,PFNE_OFFS>(output, outne);
    l1pf_pattern_unpack<NMU,     PFMU_OFFS>(output, outmu);
}

#endif

