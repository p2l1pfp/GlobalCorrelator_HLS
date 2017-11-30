#include "simple_fullpfalgo.h"
#include "mp7pf_encoding.h"
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
    ap_int<26> dpt2 = (dpt > 5792) ? ap_int<26>((1<<25)-1) : ap_int<26>(dpt*dpt);
    int dq = dr2 + (dpt2*ptscale >> 8);
    return ((dr2 < int(dr2max)) && (dq < int(max))) ? ap_uint<NB>(dq) : max;
}

template<typename T, int NIn, int NOut>
void ptsort_hwopt(T in[NIn], T out[NOut]) {
    T tmp[NOut];
    #pragma HLS ARRAY_PARTITION variable=tmp complete

    for (int iout = 0; iout < NOut; ++iout) {
        #pragma HLS unroll
        tmp[iout].hwPt = 0;
    }

    for (int it = 0; it < NIn; ++it) {
        for (int iout = NOut-1; iout >= 0; --iout) {
            if (tmp[iout].hwPt <= in[it].hwPt) {
                if (iout == 0 || tmp[iout-1].hwPt > in[it].hwPt) {
                    tmp[iout] = in[it];
                } else {
                    tmp[iout] = tmp[iout-1];
                }
            }
        }

    }
    for (int iout = 0; iout < NOut; ++iout) {
        out[iout] = tmp[iout];
    }

}



template<int DR2MAX>
void tk2em_drvals(EmCaloObj calo[NEMCALO], TkObj track[NTRACK], tk2em_dr_t calo_track_drval[NTRACK][NEMCALO], bool isMu[NTRACK]) {
    const tk2em_dr_t eDR2MAX = DR2MAX;
    for (int it = 0; it < NTRACK; ++it) {
        for (int icalo = 0; icalo < NEMCALO; ++icalo) {
            if (isMu[it]){calo_track_drval[it][icalo] = eDR2MAX; } // set to DR max if the track is a muon
            else { calo_track_drval[it][icalo] = dr2_int_cap(track[it].hwEta, track[it].hwPhi, calo[icalo].hwEta, calo[icalo].hwPhi, eDR2MAX); } 
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
    	int tmp = (DR2MAX<<8)/(i?i*i:1), int18_max = (1<<17)-1;
        vals[i] = (tmp > int18_max ? int18_max : tmp);
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
    // now, DR2MAX is 10 bits, so dptscale max is at most 10+8 bits = 18 bits
    for (int it = 0; it < NTRACK; ++it) {
        pt_t caloPtMin = track[it].hwPt - 2*(track[it].hwPtErr);
        ap_int<18> dptscale  = calc_dptscale<DR2MAX>(track[it].hwPtErr);
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

void tk2em_link(EmCaloObj calo[NEMCALO], TkObj track[NTRACK], ap_uint<NEMCALO> calo_track_link_bit[NTRACK], bool isMu[NTRACK]) {
    const int DR2MAX = PFALGO3_DR2MAX_TK_EM;
    tk2em_dr_t drvals[NTRACK][NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0

    tk2em_drvals<DR2MAX>(calo, track, drvals, isMu);
    pick_closest<DR2MAX,NTRACK,NEMCALO,tk2em_dr_t>(drvals, calo_track_link_bit);

    // for (int it = 0; it < NTRACK; it++){
    //     printf("HLS: tk index = %i and match = %i \n", it, int(calo_track_link_bit[it]));        
    // }

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
        tkerr2[it] = (track[it].hwPtErr * track[it].hwPtErr);
    }
}
void tk2calo_sumtk(TkObj track[NTRACK], bool isEle[NTRACK], bool isMu[NTRACK], int tkerr2[NTRACK], ap_uint<NCALO> calo_track_link_bit[NTRACK], pt_t sumtk[NCALO], int sumtkerr2[NCALO]) {
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

void tk2em_sumtk(TkObj track[NTRACK], ap_uint<NEMCALO> calo_track_link_bit[NTRACK], pt_t sumtk[NEMCALO]) {
    for (int icalo = 0; icalo < NEMCALO; ++icalo) {
        pt_t sum = 0;
        for (int it = 0; it < NTRACK; ++it) {
            if (calo_track_link_bit[it][icalo]) { sum += track[it].hwPt; }
        }
        sumtk[icalo] = sum;
        // printf("HLS: emcalo index = %i and sumtk = %i \n", icalo, int(sumtk[icalo]));                
    }
}

void em2calo_sumem(EmCaloObj emcalo[NEMCALO], bool isEM[NEMCALO], ap_uint<NCALO> em_had_link_bit[NTRACK], pt_t sumem[NEMCALO], bool keepcalo[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t sum = 0; bool keep = false;
        for (int iem = 0; iem < NEMCALO; ++iem) {
            if (em_had_link_bit[iem][icalo]) { 
                if (isEM[iem]) sum += emcalo[iem].hwPt; 
                else keep = false;
            }
        }
        sumem[icalo] = sum;
        keepcalo[icalo] = keep;
    }
}

void tk2calo_tkalgo(TkObj track[NTRACK], bool isEle[NTRACK], bool isMu[NTRACK], ap_uint<NCALO> calo_track_link_bit[NTRACK], PFChargedObj pfout[NTRACK]) {
    const pt_t TKPT_MAX_LOOSE = PFALGO3_TK_MAXINVPT_LOOSE; // 20 * PT_SCALE;
    const pt_t TKPT_MAX_TIGHT = PFALGO3_TK_MAXINVPT_TIGHT; // 20 * PT_SCALE;
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


void tk2calo_caloalgo(HadCaloObj calo[NCALO], pt_t sumtk[NCALO], int sumtkerr2[NCALO], PFNeutralObj pfout[NCALO]) {
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
            int ptdiff2 = ptdiff*ptdiff;
            int sigma2 = (calo[icalo].hwPtErr*calo[icalo].hwPtErr);
            int sigma2Lo = (sigma2 << 2), sigma2Hi = sigma2 + (sigma2 >> 1);
            if ((ptdiff > 0 && ptdiff2 <= sigma2Hi) || (ptdiff < 0 && ptdiff2 < sigma2Lo)) {
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

void em2calo_sub(HadCaloObj calo[NCALO], pt_t sumem[NCALO], bool keepcalo[NCALO], HadCaloObj calo_out[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t ptsub = calo[icalo].hwPt   - sumem[icalo];
        pt_t emsub = calo[icalo].hwEmPt - sumem[icalo];
        if ((ptsub <= (calo[icalo].hwPt >> 4)) || 
                (calo[icalo].hwIsEM && (emsub <= (calo[icalo].hwEmPt>>3)) && !keepcalo[icalo])) {
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

//-------------------------------------------------------
// TK-MU Algos
//-------------------------------------------------------

void spfph_mu2trk_dptvals(MuObj mu[NMU], TkObj track[NTRACK], pt_t mu_track_dptval[NMU][NTRACK]) {
    const ap_uint<12> DR2MAX = PFALGO3_DR2MAX_TK_MU;
    for (int im = 0; im < NMU; ++im) {
        for (int it = 0; it < NTRACK; ++it) {
            pt_t dpt = mu[im].hwPt - track[it].hwPt;
            if (dr2_int_cap<12>(mu[im].hwEta, mu[im].hwPhi, track[it].hwEta, track[it].hwPhi, DR2MAX) < DR2MAX) {
                mu_track_dptval[im][it] = (dpt > 0 ? dpt : pt_t(-dpt));
            } else {
                mu_track_dptval[im][it] = mu[im].hwPt >> 1;
            }
        }
    }
}

void spfph_mu2trk_linkstep(MuObj mu[NMU], pt_t mu_track_dptval[NMU][NTRACK], ap_uint<NMU> mu_track_link_bit[NTRACK]) {
    for (int im = 0; im < NMU; ++im) {
        for (int it = 0; it < NTRACK; ++it) {
            pt_t mydpt = mu_track_dptval[im][it];
            bool link = (mydpt < (mu[im].hwPt >> 1));
            for (int j = 0; j < NTRACK; ++j) {
                if (it <= j) link = link && (mu_track_dptval[im][j] >= mydpt);
                else         link = link && (mu_track_dptval[im][j] >  mydpt);
            }   
            mu_track_link_bit[it][im] = link;
        }
    }
}

void spfph_mutrk_link(MuObj mu[NMU], TkObj track[NTRACK], ap_uint<NMU> mu_track_link_bit[NTRACK]) {
    
    #pragma HLS ARRAY_PARTITION variable=mu complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu_track_link_bit complete dim=0

    pt_t dptvals[NMU][NTRACK];
    #pragma HLS ARRAY_PARTITION variable=dptvals complete dim=0

    spfph_mu2trk_dptvals(mu, track, dptvals);
    spfph_mu2trk_linkstep(mu, dptvals, mu_track_link_bit);
}

void spfph_mualgo(MuObj mu[NMU], TkObj track[NTRACK], ap_uint<NMU> mu_track_link_bit[NTRACK], PFChargedObj pfmuout[NMU], bool isMu[NTRACK]) {
    #pragma HLS ARRAY_PARTITION variable=isMu complete

    for (int im = 0; im < NMU; ++im) {
        bool good = false;
        int ibest = -1;
        for (int it = 0; it < NTRACK; ++it) {
            if (mu_track_link_bit[it][im]){ good = true; ibest = it; }
        }
        if (mu[im].hwPt > 0 && good && ibest != -1) {
            pfmuout[im].hwPt  = track[ibest].hwPt;
            pfmuout[im].hwEta = track[ibest].hwEta;
            pfmuout[im].hwPhi = track[ibest].hwPhi;
            pfmuout[im].hwId  = PID_Muon;
            pfmuout[im].hwZ0  = track[ibest].hwZ0;
            isMu[ibest] = 1;
        } else {
            pfmuout[im].hwPt  = 0;
            pfmuout[im].hwEta = 0;
            pfmuout[im].hwPhi = 0;
            pfmuout[im].hwId  = 0;
            pfmuout[im].hwZ0  = 0;
        }
    }
}

//-------------------------------------------------------
// PF Algos
//-------------------------------------------------------

void pfalgo3_full(EmCaloObj calo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) {
    
    #pragma HLS ARRAY_PARTITION variable=calo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete    
    #pragma HLS ARRAY_PARTITION variable=outch complete
    #pragma HLS ARRAY_PARTITION variable=outpho complete
    #pragma HLS ARRAY_PARTITION variable=outne complete
    #pragma HLS ARRAY_PARTITION variable=outmu complete

    #pragma HLS pipeline II=HLS_pipeline_II

    // ---------------------------------------------------------------
    // TK-MU Linking
    ap_uint<NMU> mu_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=mu_track_link_bit complete
    bool isMu[NTRACK];
    for (int it = 0; it < NTRACK; ++it) { isMu[it] = 0; }

    spfph_mutrk_link(mu, track, mu_track_link_bit);
    spfph_mualgo(mu, track, mu_track_link_bit, outmu, isMu);

    // ---------------------------------------------------------------
    // TK-EM Linking
    ap_uint<NEMCALO> em_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=em_track_link_bit complete
    tk2em_link(calo, track, em_track_link_bit, isMu);

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

    bool keepcalo[NCALO];
    pt_t sumem[NCALO]; 
    #pragma HLS ARRAY_PARTITION variable=sumem complete
    em2calo_sumem(calo, isEM, em_calo_link_bit, sumem, keepcalo);

    HadCaloObj hadcalo_sub[NCALO];
    #pragma HLS ARRAY_PARTITION variable=hadcalo_sub complete

    em2calo_sub(hadcalo, sumem, keepcalo, hadcalo_sub);

    // ---------------------------------------------------------------
    // TK-HAD Linking
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

    tk2calo_tkalgo(track, isEle, isMu, calo_track_link_bit, outch);
    tk2calo_sumtk(track, isEle, isMu, tkerr2, calo_track_link_bit, sumtk, sumtkerr2);

    PFNeutralObj outne_all[NCALO];
    #pragma HLS ARRAY_PARTITION variable=outne_all complete
    tk2calo_caloalgo(hadcalo_sub, sumtk, sumtkerr2, outne_all);
    ptsort_hwopt<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);

}


void mp7wrapped_pack_in(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], MP7DataWord data[MP7_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=emcalo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete
    // pack inputs
    assert(2*NEMCALO + 2*NTRACK + 2*NCALO + 2*NMU <= MP7_NCHANN);
    #define HADOFFS 2*NEMCALO
    #define TKOFFS 2*NCALO+HADOFFS
    #define MUOFFS 2*NTRACK+TKOFFS
    mp7_pack<NEMCALO,0>(emcalo, data);
    mp7_pack<NCALO,HADOFFS>(hadcalo, data);
    mp7_pack<NTRACK,TKOFFS>(track, data);
    mp7_pack<NMU,MUOFFS>(mu, data);
}

void mp7wrapped_unpack_in(MP7DataWord data[MP7_NCHANN], EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=emcalo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete
    // unpack inputs
    assert(2*NEMCALO + 2*NTRACK + 2*NCALO + 2*NMU <= MP7_NCHANN);
    #define HADOFFS 2*NEMCALO
    #define TKOFFS 2*NCALO+HADOFFS
    #define MUOFFS 2*NTRACK+TKOFFS
    mp7_unpack<NEMCALO,0>(data, emcalo);
    mp7_unpack<NCALO,HADOFFS>(data, hadcalo);
    mp7_unpack<NTRACK,TKOFFS>(data, track);
    mp7_unpack<NMU,MUOFFS>(data, mu);
}

void mp7wrapped_pack_out( PFChargedObj pfch[NTRACK], PFNeutralObj pfpho[NPHOTON], PFNeutralObj pfne[NSELCALO], PFChargedObj pfmu[NMU], MP7DataWord data[MP7_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfpho complete
    #pragma HLS ARRAY_PARTITION variable=pfne complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete
    // pack outputs
    assert(2*NTRACK + 2*NPHOTON + 2*NSELCALO + 2*NMU <= MP7_NCHANN);
    #define PHOOFFS 2*NTRACK
    #define NHOFFS 2*NPHOTON+PHOOFFS
    #define PFMUOFFS 2*NSELCALO+NHOFFS
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
    for (unsigned int i = 0; i < NMU; ++i) {
        data[2*i+0+PFMUOFFS] = ( pfmu[i].hwId, pfmu[i].hwPt );
        data[2*i+1+PFMUOFFS] = ( pfmu[i].hwZ0, pfmu[i].hwPhi, pfmu[i].hwEta );
    }

}
void mp7wrapped_unpack_out( MP7DataWord data[MP7_NCHANN], PFChargedObj pfch[NTRACK], PFNeutralObj pfpho[NPHOTON], PFNeutralObj pfne[NSELCALO], PFChargedObj pfmu[NMU]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfpho complete
    #pragma HLS ARRAY_PARTITION variable=pfne complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete
    // unpack outputs
    assert(2*NTRACK + 2*NPHOTON + 2*NSELCALO + 2*NMU <= MP7_NCHANN);
    #define PHOOFFS 2*NTRACK
    #define NHOFFS 2*NPHOTON+PHOOFFS
    #define PFMUOFFS 2*NSELCALO+NHOFFS
    for (unsigned int i = 0; i < NTRACK; ++i) {
        pfch[i].hwPt  = data[2*i+0](15, 0);
        pfch[i].hwId  = data[2*i+0](18,16);
        pfch[i].hwEta = data[2*i+1](9, 0);
        pfch[i].hwPhi = data[2*i+1](19,10);
        pfch[i].hwZ0  = data[2*i+1](29,20);
    }
    for (unsigned int i = 0; i < NPHOTON; ++i) {
        pfpho[i].hwPt  = data[2*i+0+PHOOFFS](15, 0);
        pfpho[i].hwId  = data[2*i+0+PHOOFFS](18,16);
        pfpho[i].hwEta = data[2*i+1+PHOOFFS](9, 0);
        pfpho[i].hwPhi = data[2*i+1+PHOOFFS](19,10);
    }
    for (unsigned int i = 0; i < NSELCALO; ++i) {
        pfne[i].hwPt  = data[2*i+0+NHOFFS](15, 0);
        pfne[i].hwId  = data[2*i+0+NHOFFS](18,16);
        pfne[i].hwEta = data[2*i+1+NHOFFS](9, 0);
        pfne[i].hwPhi = data[2*i+1+NHOFFS](19,10);
    }
    for (unsigned int i = 0; i < NMU; ++i) {
        pfmu[i].hwPt  = data[2*i+0+PFMUOFFS](15, 0);
        pfmu[i].hwId  = data[2*i+0+PFMUOFFS](18,16);
        pfmu[i].hwEta = data[2*i+1+PFMUOFFS](9, 0);
        pfmu[i].hwPhi = data[2*i+1+PFMUOFFS](19,10);
        pfmu[i].hwZ0  = data[2*i+1+PFMUOFFS](29,20);
    }

}

void mp7wrapped_pfalgo3_full(MP7DataWord input[MP7_NCHANN], MP7DataWord output[MP7_NCHANN]) {
    
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=output complete
    #pragma HLS INTERFACE ap_none port=output

    #pragma HLS pipeline II=HLS_pipeline_II

    EmCaloObj emcalo[NEMCALO]; HadCaloObj hadcalo[NCALO]; TkObj track[NTRACK]; MuObj mu[NMU];
    #pragma HLS ARRAY_PARTITION variable=emcalo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete

    PFChargedObj pfch[NTRACK]; PFNeutralObj pfpho[NPHOTON]; PFNeutralObj pfne[NSELCALO]; PFChargedObj pfmu[NMU];
    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfpho complete
    #pragma HLS ARRAY_PARTITION variable=pfne complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete

    mp7wrapped_unpack_in(input, emcalo, hadcalo, track, mu);
    pfalgo3_full(emcalo, hadcalo, track, mu, pfch, pfpho, pfne, pfmu);
    mp7wrapped_pack_out(pfch, pfpho, pfne, pfmu, output);

}

