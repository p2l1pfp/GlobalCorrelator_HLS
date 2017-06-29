#include "simple_pfalgo3.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

typedef ap_uint<7> tk2em_dr_t;
typedef ap_uint<10> tk2calo_dr_t;
typedef ap_uint<10> em2calo_dr_t;
typedef ap_uint<12> tk2calo_dq_t;

int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}
template<int NB>
ap_uint<NB> dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    int dr2 = deta*deta + dphi*dphi;
    return (dr2 < int(max) ? ap_uint<NB>(dr2) : max);
}
template<int NB, typename PTS_t>
ap_uint<NB> dr2_dpt_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, pt_t pt1, pt_t pt2, PTS_t ptscale, ap_uint<NB> dr2max, ap_uint<NB> max) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    int dr2 = deta*deta + dphi*dphi;
    pt_t dpt = pt1 - pt2;
    if (dpt < 0) dpt = 0;
    int dq = dr2 + ((dpt*dpt*ptscale) >> 8);
    return ((dr2 < int(dr2max)) && (dq < int(max))) ? ap_uint<NB>(dq) : max;
}


template<int DR2MAX>
void tk2em_drvals(EmCaloObj calo[NEMCALO], TkObj track[NTRACK], tk2em_dr_t calo_track_drval[NTRACK][NCALO]) {
    const tk2em_dr_t eDR2MAX = DR2MAX;
    for (int it = 0; it < NTRACK; ++it) {
        for (int icalo = 0; icalo < NEMCALO; ++icalo) {
            calo_track_drval[it][icalo] = dr2_int_cap(track[it].hwEta, track[it].hwPhi, calo[icalo].hwEta, calo[icalo].hwPhi, eDR2MAX);
        }
    }
}

template<int DR2MAX>
void tk2calo_drvals(HadCaloObj calo[NCALO], TkObj track[NTRACK], tk2calo_dr_t calo_track_drval[NTRACK][NCALO]) {
    const tk2calo_dr_t eDR2MAX = DR2MAX;
    for (int it = 0; it < NTRACK; ++it) {
        pt_t caloPtMin = track[it].hwPt - 2*(track[it].hwPtErr);
        if (caloPtMin < 0) caloPtMin = 0;
        for (int icalo = 0; icalo < NCALO; ++icalo) {
            if (calo[icalo].hwPt > caloPtMin) {
                calo_track_drval[it][icalo] = dr2_int_cap(track[it].hwEta, track[it].hwPhi, calo[icalo].hwEta, calo[icalo].hwPhi, eDR2MAX);
            } else {
                calo_track_drval[it][icalo] = eDR2MAX;
            }
        }
    }
}

template<int DR2MAX>
void init_dr2max_times_pterr2_inv(int vals[512]) {
    for (int i = 0; i < 512; ++i) {
        vals[i] = (DR2MAX<<8)/(i?i*i:1);
    }
}
template<int DR2MAX>
int calc_dptscale(pt_t trackHwPtErr) {
    // LUT for 1/ptErr2
    int _dr2max_times_pterr2_inv_vals[512];
    init_dr2max_times_pterr2_inv<DR2MAX>(_dr2max_times_pterr2_inv_vals);
    if (trackHwPtErr < 512) {
        return _dr2max_times_pterr2_inv_vals[trackHwPtErr];
    } else {
        return 0;
    }
}

template<int DR2MAX>
void tk2calo_drdptvals(HadCaloObj calo[NCALO], TkObj track[NTRACK], tk2calo_dq_t calo_track_drval[NTRACK][NCALO]) {
    const tk2calo_dq_t eDR2MAX = DR2MAX;
    const tk2calo_dq_t eDQMAX  = 5*DR2MAX; // at most we're 2 sigma away in pt, so that's a factor 4
    for (int it = 0; it < NTRACK; ++it) {
        pt_t caloPtMin = track[it].hwPt - 2*(track[it].hwPtErr);
        int  dptscale  = calc_dptscale<DR2MAX>(track[it].hwPtErr);
        if (caloPtMin < 0) caloPtMin = 0;
        for (int icalo = 0; icalo < NCALO; ++icalo) {
            if (calo[icalo].hwPt > caloPtMin) {
                calo_track_drval[it][icalo] = dr2_dpt_int_cap(track[it].hwEta, track[it].hwPhi, calo[icalo].hwEta, calo[icalo].hwPhi, track[it].hwPt, calo[icalo].hwPt, dptscale, eDR2MAX, eDQMAX);
                //if (calo_track_drval[it][icalo] < eDQMAX) printf("HWO DQ(track %+7d %+7d  calo %3d) = %12d\n", int(track[it].hwEta), int(track[it].hwPhi), icalo, int(calo_track_drval[it][icalo]));
            } else {
                calo_track_drval[it][icalo] = eDQMAX;
            }
        }
    }
}

template<int DR2MAX>
void em2calo_drvals(HadCaloObj hadcalo[NCALO], EmCaloObj emcalo[NEMCALO], em2calo_dr_t hadcalo_emcalo_drval[NEMCALO][NCALO]) {
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

template<int DR2MAX, int NTK, int NCA, typename DR_T>
void pick_closest(DR_T calo_track_drval[NTK][NCA], ap_uint<NCA> calo_track_link_bit[NTK]) {
    const DR_T eDR2MAX = DR2MAX;
    for (int it = 0; it < NTK; ++it) {
        for (int icalo = 0; icalo < NCA; ++icalo) {
            DR_T mydr = calo_track_drval[it][icalo];
            bool link = (mydr != eDR2MAX);
            for (int j = 0; j < NCA; ++j) {
                if (icalo <= j) link = link && (calo_track_drval[it][j] >= mydr);
                else            link = link && (calo_track_drval[it][j] >  mydr);
            }
            calo_track_link_bit[it][icalo] = link;
        }
    }
}
void tk2calo_link_dronly(HadCaloObj calo[NCALO], TkObj track[NTRACK], ap_uint<NCALO> calo_track_link_bit[NTRACK]) {
    const int DR2MAX = PFALGO3_DR2MAX_TK_CALO;
    tk2calo_dr_t drvals[NTRACK][NCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0

    tk2calo_drvals<DR2MAX>(calo, track, drvals);
    pick_closest<DR2MAX,NTRACK,NCALO,tk2calo_dr_t>(drvals, calo_track_link_bit);
}
void tk2calo_link_drdpt(HadCaloObj calo[NCALO], TkObj track[NTRACK], ap_uint<NCALO> calo_track_link_bit[NTRACK]) {
    const int DR2MAX = PFALGO3_DR2MAX_TK_CALO;
    const int DQMAX = 5*DR2MAX;
    tk2calo_dq_t drvals[NTRACK][NCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0

    tk2calo_drdptvals<DR2MAX>(calo, track, drvals);
    pick_closest<DQMAX,NTRACK,NCALO,tk2calo_dq_t>(drvals, calo_track_link_bit);
}

void tk2em_link(EmCaloObj calo[NEMCALO], TkObj track[NTRACK], ap_uint<NEMCALO> calo_track_link_bit[NTRACK]) {
    const int DR2MAX = PFALGO3_DR2MAX_TK_EM;
    tk2em_dr_t drvals[NTRACK][NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0

    tk2em_drvals<DR2MAX>(calo, track, drvals);
    pick_closest<DR2MAX,NTRACK,NEMCALO,tk2em_dr_t>(drvals, calo_track_link_bit);
}
void em2calo_link(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], ap_uint<NCALO> em_calo_link_bit[NCALO]) {
    const int DR2MAX = PFALGO3_DR2MAX_EM_CALO;
    em2calo_dr_t drvals[NEMCALO][NCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0

    em2calo_drvals<DR2MAX>(hadcalo, emcalo, drvals);
    pick_closest<DR2MAX,NEMCALO,NCALO,em2calo_dr_t>(drvals, em_calo_link_bit);
}



void tk2calo_tkerr2(TkObj track[NTRACK], int tkerr2[NTRACK]) {
    for (int it = 0; it < NTRACK; ++it) {
        tkerr2[it] = (track[it].hwPtErr * track[it].hwPtErr) << 2; // we will want (2*error)^2
    }
}
void tk2calo_sumtk(TkObj track[NTRACK], int tkerr2[NTRACK], ap_uint<NCALO> calo_track_link_bit[NTRACK], pt_t sumtk[NCALO], int sumtkerr2[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t sum = 0;
        int sumerr = 0;
        for (int it = 0; it < NTRACK; ++it) {
            if (calo_track_link_bit[it][icalo]) { sum += track[it].hwPt; sumerr += tkerr2[it]; }
        }
        sumtk[icalo] = sum;
        sumtkerr2[icalo] = sumerr;
    }
}
void tk2calo_sumtk(TkObj track[NTRACK], bool isEle[NTRACK], int tkerr2[NTRACK], ap_uint<NCALO> calo_track_link_bit[NTRACK], pt_t sumtk[NCALO], int sumtkerr2[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t sum = 0;
        int sumerr = 0;
        for (int it = 0; it < NTRACK; ++it) {
            if (!isEle[it] && calo_track_link_bit[it][icalo]) { sum += track[it].hwPt; sumerr += tkerr2[it]; }
        }
        sumtk[icalo] = sum;
        sumtkerr2[icalo] = sumerr;
    }
}

void tk2em_sumtk(TkObj track[NTRACK], ap_uint<NEMCALO> calo_track_link_bit[NTRACK], pt_t sumtk[NEMCALO]) {
    for (int icalo = 0; icalo < NEMCALO; ++icalo) {
        pt_t sum = 0;
        for (int it = 0; it < NTRACK; ++it) {
            if (calo_track_link_bit[it][icalo]) { sum += track[it].hwPt; }
        }
        sumtk[icalo] = sum;
    }
}

void em2calo_sumem(EmCaloObj emcalo[NEMCALO], bool isEM[NEMCALO], ap_uint<NCALO> em_had_link_bit[NTRACK], pt_t sumem[NEMCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t sum = 0;
        for (int iem = 0; iem < NEMCALO; ++iem) {
            if (isEM[iem] && em_had_link_bit[iem][icalo]) { sum += emcalo[iem].hwPt; }
        }
        sumem[icalo] = sum;
    }
}



void tk2calo_tkalgo(TkObj track[NTRACK], ap_uint<NCALO> calo_track_link_bit[NTRACK], PFChargedObj pfout[NTRACK]) {
    const pt_t TKPT_MAX = PFALGO3_TK_MAXINVPT; // 20 * PT_SCALE;
    for (int it = 0; it < NTRACK; ++it) {
        bool good = (track[it].hwPt < TKPT_MAX) || calo_track_link_bit[it].or_reduce();
        if (good) {
            pfout[it].hwPt  = track[it].hwPt;
            pfout[it].hwEta = track[it].hwEta;
            pfout[it].hwPhi = track[it].hwPhi;
            pfout[it].hwId  = PID_Charged;
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
void tk2calo_tkalgo(TkObj track[NTRACK], bool isEle[NTRACK], ap_uint<NCALO> calo_track_link_bit[NTRACK], PFChargedObj pfout[NTRACK]) {
    const pt_t TKPT_MAX = PFALGO3_TK_MAXINVPT; // 20 * PT_SCALE;
    for (int it = 0; it < NTRACK; ++it) {
        bool good = isEle[it] || (track[it].hwPt < TKPT_MAX) || calo_track_link_bit[it].or_reduce();
        if (good) {
            pfout[it].hwPt  = track[it].hwPt;
            pfout[it].hwEta = track[it].hwEta;
            pfout[it].hwPhi = track[it].hwPhi;
            pfout[it].hwId  = isEle[it] ? PID_Electron : PID_Charged;
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


void tk2calo_caloalgo(HadCaloObj calo[NCALO], pt_t sumtk[NCALO], int sumtkerr2[NCALO], PFNeutralObj pfout[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t calopt;
        if (sumtk[icalo] == 0) {
            calopt = calo[icalo].hwPt;
        } else {
            pt_t ptdiff = calo[icalo].hwPt - sumtk[icalo];
            if (ptdiff > 0 && (ptdiff*ptdiff) > sumtkerr2[icalo]) {
                calopt = ptdiff;
            } else {
                calopt = 0;
            }
        }
        pfout[icalo].hwPt  = calopt;
        pfout[icalo].hwEta = calopt ? calo[icalo].hwEta : etaphi_t(0);
        pfout[icalo].hwPhi = calopt ? calo[icalo].hwPhi : etaphi_t(0);
        pfout[icalo].hwId  = calopt ? PID_Neutral : 0;
    }
}
void tk2em_emalgo(EmCaloObj calo[NEMCALO], pt_t sumtk[NEMCALO], bool isEM[NEMCALO], pt_t photonPt[NEMCALO]) {
    for (int icalo = 0; icalo < NEMCALO; ++icalo) {
        if (sumtk[icalo] == 0) {
            isEM[icalo] = true;
            photonPt[icalo] = calo[icalo].hwPt;
        } else {
            pt_t ptdiff = calo[icalo].hwPt - sumtk[icalo];
            if (ptdiff*ptdiff <= ((calo[icalo].hwPtErr*calo[icalo].hwPtErr)<<2)) {
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
void tk2em_elealgo(ap_uint<NEMCALO> em_track_link_bit[NTRACK], bool isEM[NEMCALO], bool isEle[NTRACK]) {
    for (int it = 0; it < NTRACK; ++it) {
        bool ele = false;
        for (int icalo = 0; icalo < NEMCALO; ++icalo) { 
            if (isEM[icalo] && em_track_link_bit[it][icalo]) ele = true;
        }
        isEle[it] = ele;
    }
}
void tk2em_photons(EmCaloObj calo[NEMCALO], pt_t photonPt[NEMCALO], PFNeutralObj pfout[NSELCALO]) {
    for (int icalo = 0; icalo < NEMCALO; ++icalo) {
        pfout[icalo].hwPt  = photonPt[icalo];
        pfout[icalo].hwEta = photonPt[icalo] ? calo[icalo].hwEta : etaphi_t(0);
        pfout[icalo].hwPhi = photonPt[icalo] ? calo[icalo].hwPhi : etaphi_t(0);
        pfout[icalo].hwId  = photonPt[icalo] ? PID_Photon : 0;
    }
}

void em2calo_sub(HadCaloObj calo[NCALO], pt_t sumem[NCALO], HadCaloObj calo_out[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t ptsub = calo[icalo].hwPt   - sumem[icalo];
        pt_t emsub = calo[icalo].hwEmPt - sumem[icalo];
        if ((ptsub < (calo[icalo].hwPt >> 4)) || 
                (calo[icalo].hwIsEM && (emsub < (calo[icalo].hwEmPt>>3)))) {
            calo_out[icalo].hwPt   = 0;
            calo_out[icalo].hwEmPt = 0;
            calo_out[icalo].hwEta  = 0;
            calo_out[icalo].hwPhi  = 0;
            calo_out[icalo].hwIsEM = 0;
        } else {
            calo_out[icalo].hwPt   = ptsub;
            calo_out[icalo].hwEmPt = (emsub > 0 ? emsub : pt_t(0));
            calo_out[icalo].hwEta  = calo[icalo].hwEta;
            calo_out[icalo].hwPhi  = calo[icalo].hwPhi;
            calo_out[icalo].hwIsEM = calo[icalo].hwIsEM;
        }
    }
}

void pfalgo3_calo(HadCaloObj calo[NCALO], TkObj track[NTRACK], PFChargedObj outch[NTRACK], PFNeutralObj outne[NSELCALO]) {
    #pragma HLS ARRAY_PARTITION variable=calo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=outch complete
    #pragma HLS ARRAY_PARTITION variable=outne complete

    #pragma HLS pipeline II=5

    ap_uint<NCALO> calo_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=calo_track_link_bit complete
    tk2calo_link_drdpt(calo, track, calo_track_link_bit);
    //tk2calo_link_dronly(calo, track, calo_track_link_bit);

    int tkerr2[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=tkerr2 complete
    tk2calo_tkerr2(track, tkerr2);

    pt_t sumtk[NCALO]; int sumtkerr2[NCALO];
    #pragma HLS ARRAY_PARTITION variable=sumtk complete
    #pragma HLS ARRAY_PARTITION variable=sumtkerr2 complete

    tk2calo_tkalgo(track, calo_track_link_bit, outch);
    tk2calo_sumtk(track, tkerr2, calo_track_link_bit, sumtk, sumtkerr2);

    PFNeutralObj outne_all[NCALO];
    #pragma HLS ARRAY_PARTITION variable=outne_all complete
    tk2calo_caloalgo(calo, sumtk, sumtkerr2, outne_all);
    ptsort_hwopt<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);
}

void pfalgo3_em(EmCaloObj calo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], bool isEle[NTRACK], PFNeutralObj outpho[NPHOTON], HadCaloObj hadcalo_out[NCALO]) {
    #pragma HLS ARRAY_PARTITION variable=calo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=isEle complete
    #pragma HLS ARRAY_PARTITION variable=outpho complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo_out complete

    #pragma HLS pipeline II=5

    ap_uint<NEMCALO> em_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=em_track_link_bit complete
    tk2em_link(calo, track, em_track_link_bit);

    pt_t sumtk[NEMCALO]; 
    #pragma HLS ARRAY_PARTITION variable=sumtk complete

    pt_t photonPt[NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=photonPt complete

    bool isEM[NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=isEM complete

    tk2em_sumtk(track, em_track_link_bit, sumtk);
    tk2em_emalgo(calo, sumtk, isEM, photonPt);
    tk2em_photons(calo, photonPt, outpho);
    tk2em_elealgo(em_track_link_bit, isEM, isEle);

    ap_uint<NCALO> em_calo_link_bit[NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=em_calo_link_bit complete
    em2calo_link(calo, hadcalo, em_calo_link_bit);

    pt_t sumem[NCALO]; 
    #pragma HLS ARRAY_PARTITION variable=sumem complete
    em2calo_sumem(calo, isEM, em_calo_link_bit, sumem);
    em2calo_sub(hadcalo, sumem, hadcalo_out);
}

void pfalgo3_full(EmCaloObj calo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO]) {
    #pragma HLS ARRAY_PARTITION variable=calo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=outch complete
    #pragma HLS ARRAY_PARTITION variable=outpho complete
    #pragma HLS ARRAY_PARTITION variable=outne complete

    #pragma HLS pipeline II=5

    ap_uint<NEMCALO> em_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=em_track_link_bit complete
    tk2em_link(calo, track, em_track_link_bit);

    pt_t sumtk2em[NEMCALO]; 
    #pragma HLS ARRAY_PARTITION variable=sumtk2em complete

    pt_t photonPt[NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=photonPt complete

    bool isEM[NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=isEM complete

    tk2em_sumtk(track, em_track_link_bit, sumtk2em);
    tk2em_emalgo(calo, sumtk2em, isEM, photonPt);
    tk2em_photons(calo, photonPt, outpho);

    bool isEle[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=isEle complete
    tk2em_elealgo(em_track_link_bit, isEM, isEle);

    ap_uint<NCALO> em_calo_link_bit[NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=em_calo_link_bit complete
    em2calo_link(calo, hadcalo, em_calo_link_bit);

    pt_t sumem[NCALO]; 
    #pragma HLS ARRAY_PARTITION variable=sumem complete
    em2calo_sumem(calo, isEM, em_calo_link_bit, sumem);

    HadCaloObj hadcalo_sub[NCALO];
    #pragma HLS ARRAY_PARTITION variable=hadcalo_sub complete

    em2calo_sub(hadcalo, sumem, hadcalo_sub);

    ap_uint<NCALO> calo_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=calo_track_link_bit complete
    tk2calo_link_drdpt(hadcalo_sub, track, calo_track_link_bit);
    //tk2calo_link_dronly(hadcalo_sub, track, calo_track_link_bit);

    int tkerr2[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=tkerr2 complete
    tk2calo_tkerr2(track, tkerr2);

    pt_t sumtk[NCALO]; int sumtkerr2[NCALO];
    #pragma HLS ARRAY_PARTITION variable=sumtk complete
    #pragma HLS ARRAY_PARTITION variable=sumtkerr2 complete

    tk2calo_tkalgo(track, isEle, calo_track_link_bit, outch);
    tk2calo_sumtk(track, isEle, tkerr2, calo_track_link_bit, sumtk, sumtkerr2);

    PFNeutralObj outne_all[NCALO];
    #pragma HLS ARRAY_PARTITION variable=outne_all complete
    tk2calo_caloalgo(hadcalo_sub, sumtk, sumtkerr2, outne_all);
    ptsort_hwopt<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);
}

void mp7wrapped_pack_in(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MP7DataWord data[MP7_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=emcalo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    // pack inputs
    assert(2*NEMCALO + 2*NTRACK + 2*NCALO <= MP7_NCHANN);
    #define HADOFFS 2*NEMCALO
    #define TKOFFS 2*NCALO+HADOFFS
    for (unsigned int i = 0; i < NEMCALO; ++i) {
        data[2*i+0] = ( emcalo[i].hwPtErr, emcalo[i].hwPt );
        data[2*i+1] = ( emcalo[i].hwPhi,   emcalo[i].hwEta );
    }
    for (unsigned int i = 0; i < NCALO; ++i) {
        data[2*i+0+HADOFFS] = ( hadcalo[i].hwEmPt, hadcalo[i].hwPt );
        data[2*i+1+HADOFFS] = ( hadcalo[i].hwIsEM, hadcalo[i].hwPhi, hadcalo[i].hwEta );
    }
    for (unsigned int i = 0; i < NTRACK; ++i) {
        data[2*i+0+TKOFFS] = ( track[i].hwPtErr, track[i].hwPt );
        data[2*i+1+TKOFFS] = ( track[i].hwZ0, track[i].hwPhi, track[i].hwEta );
   }
}
void mp7wrapped_unpack_in(MP7DataWord data[MP7_NCHANN], EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=emcalo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    // unpack inputs
    assert(2*NEMCALO + 2*NTRACK + 2*NCALO <= MP7_NCHANN);
    #define HADOFFS 2*NEMCALO
    #define TKOFFS 2*NCALO+HADOFFS
    for (unsigned int i = 0; i < NEMCALO; ++i) {
        emcalo[i].hwPt    = data[2*i+0](15, 0);
        emcalo[i].hwPtErr = data[2*i+0](31,16);
        emcalo[i].hwEta   = data[2*i+1](8, 0);
        emcalo[i].hwPhi   = data[2*i+1](17,9);
    }
    for (unsigned int i = 0; i < NCALO; ++i) {
        hadcalo[i].hwPt   = data[2*i+0+HADOFFS](15, 0);
        hadcalo[i].hwEmPt = data[2*i+0+HADOFFS](31,16);
        hadcalo[i].hwEta  = data[2*i+1+HADOFFS](8, 0);
        hadcalo[i].hwPhi  = data[2*i+1+HADOFFS](17,9);
        hadcalo[i].hwIsEM = data[2*i+1+HADOFFS][18];
    }
    for (unsigned int i = 0; i < NTRACK; ++i) {
        track[i].hwPt    = data[2*i+0+TKOFFS](15, 0);
        track[i].hwPtErr = data[2*i+0+TKOFFS](31,16);
        track[i].hwEta   = data[2*i+1+TKOFFS](8, 0);
        track[i].hwPhi   = data[2*i+1+TKOFFS](17,9);
        track[i].hwZ0    = data[2*i+1+TKOFFS](28,18);
    }

}
void mp7wrapped_pack_out( PFChargedObj pfch[NTRACK], PFNeutralObj pfpho[NPHOTON], PFNeutralObj pfne[NSELCALO], MP7DataWord data[MP7_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfpho complete
    #pragma HLS ARRAY_PARTITION variable=pfne complete
    // pack outputs
    assert(2*NTRACK + 2*NPHOTON + 2*NSELCALO <= MP7_NCHANN);
    #define PHOOFFS 2*NTRACK
    #define NHOFFS 2*NPHOTON+PHOOFFS
    for (unsigned int i = 0; i < NTRACK; ++i) {
        data[2*i+0] = ( pfch[i].hwId,  pfch[i].hwPt );
        data[2*i+1] = ( pfch[i].hwZ0, pfch[i].hwPhi, pfch[i].hwEta );
    }
    for (unsigned int i = 0; i < NPHOTON; ++i) {
        data[2*i+0+PHOOFFS] = ( pfpho[i].hwId, pfpho[i].hwPt );
        data[2*i+1+PHOOFFS] = ( pfpho[i].hwPhi, pfpho[i].hwEta );
    }
    for (unsigned int i = 0; i < NSELCALO; ++i) {
        data[2*i+0+NHOFFS] = ( pfne[i].hwId, pfne[i].hwPt );
        data[2*i+1+NHOFFS] = ( pfne[i].hwPhi, pfne[i].hwEta );
    }
}
void mp7wrapped_unpack_out( MP7DataWord data[MP7_NCHANN], PFChargedObj pfch[NTRACK], PFNeutralObj pfpho[NPHOTON], PFNeutralObj pfne[NSELCALO]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfpho complete
    #pragma HLS ARRAY_PARTITION variable=pfne complete
    // unpack outputs
    assert(2*NTRACK + 2*NPHOTON + 2*NSELCALO <= MP7_NCHANN);
    #define PHOOFFS 2*NTRACK
    #define NHOFFS 2*NPHOTON+PHOOFFS
    for (unsigned int i = 0; i < NTRACK; ++i) {
        pfch[i].hwPt  = data[2*i+0](15, 0);
        pfch[i].hwId  = data[2*i+0](17,16);
        pfch[i].hwEta = data[2*i+1](8, 0);
        pfch[i].hwPhi = data[2*i+1](17,9);
        pfch[i].hwZ0  = data[2*i+1](28,18);
    }
    for (unsigned int i = 0; i < NPHOTON; ++i) {
        pfpho[i].hwPt  = data[2*i+0+PHOOFFS](15, 0);
        pfpho[i].hwId  = data[2*i+0+PHOOFFS](17,16);
        pfpho[i].hwEta = data[2*i+1+PHOOFFS](8, 0);
        pfpho[i].hwPhi = data[2*i+1+PHOOFFS](17,9);
    }
    for (unsigned int i = 0; i < NSELCALO; ++i) {
        pfne[i].hwPt  = data[2*i+0+NHOFFS](15, 0);
        pfne[i].hwId  = data[2*i+0+NHOFFS](17,16);
        pfne[i].hwEta = data[2*i+1+NHOFFS](8, 0);
        pfne[i].hwPhi = data[2*i+1+NHOFFS](17,9);
    }
}

void mp7wrapped_pfalgo3_full(MP7DataWord input[MP7_NCHANN], MP7DataWord output[MP7_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=output complete

    #pragma HLS pipeline II=5

    EmCaloObj emcalo[NEMCALO]; HadCaloObj hadcalo[NCALO]; TkObj track[NTRACK]; 
    #pragma HLS ARRAY_PARTITION variable=emcalo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete

    PFChargedObj pfch[NTRACK]; PFNeutralObj pfpho[NPHOTON]; PFNeutralObj pfne[NSELCALO];
    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfpho complete
    #pragma HLS ARRAY_PARTITION variable=pfne complete

    mp7wrapped_unpack_in(input, emcalo, hadcalo, track);
    pfalgo3_full(emcalo, hadcalo, track, pfch, pfpho, pfne);
    mp7wrapped_pack_out(pfch, pfpho, pfne, output);
}

