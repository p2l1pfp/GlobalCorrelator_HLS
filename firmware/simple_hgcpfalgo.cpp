#include "simple_hgcpfalgo.h"
#include "../puppi/firmware/simple_puppi.h"
#include "mp7pf_encoding.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

//typedef ap_uint<7> tk2em_dr_t;
//typedef ap_uint<10> tk2calo_dr_t;
typedef ap_uint<10> em2calo_dr_t;
typedef ap_uint<12> tk2calo_dq_t;

int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}
template<int NB>
ap_uint<NB> dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) {
    //hardcode for etaphi size
    int tmpe = eta1-eta2;
    ap_uint<NB> deta = (tmpe > 0 ? tmpe : -tmpe);
    int tmpp = phi1-phi2;
    ap_uint<NB> dphi = (tmpp > 0 ? tmpp : -tmpp);
    int dr2 = max;
    if ((deta >> (NB/2))==0 && (dphi >> (NB/2))==0) {
        ap_uint<NB> deta2 = deta*deta;
        ap_uint<NB> dphi2 = dphi*dphi;
        dr2 = deta2 + dphi2;
    }
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
template<int NB, typename PTS_t>
ap_uint<NB> dr2_plus_dpt_int_cap(int dr2, pt_t pt1, pt_t pt2, PTS_t ptscale, ap_uint<NB> dr2max, ap_uint<NB> max) {
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

template<typename T, int NIn, int NOut>
void ptsort_hwopt_ind(T in[NIn], T out[NOut], ap_uint<6> indout[NOut]) {
    T tmp[NOut];
    #pragma HLS ARRAY_PARTITION variable=tmp complete
    ap_uint<6> tmpind[NOut];
    #pragma HLS ARRAY_PARTITION variable=tmpind complete

    for (int iout = 0; iout < NOut; ++iout) {
        #pragma HLS unroll
        tmp[iout].hwPt = 0;
        tmpind[iout] = ap_uint<6>(iout);
    }

    for (int it = 0; it < NIn; ++it) {
        for (int iout = NOut-1; iout >= 0; --iout) {
            if (tmp[iout].hwPt <= in[it].hwPt) {
                if (iout == 0 || tmp[iout-1].hwPt > in[it].hwPt) {
                    tmp[iout] = in[it];
                    tmpind[iout] = ap_uint<6>(it);
                } else {
                    tmp[iout] = tmp[iout-1];
                    tmpind[iout] = tmpind[iout-1];
                }
            }
        }

    }
    for (int iout = 0; iout < NOut; ++iout) {
        out[iout] = tmp[iout];
        indout[iout] = tmpind[iout];
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
void tk2calo_drvals_cut(HadCaloObj calo[NCALO], TkObj track[NTRACK], tk2calo_dr_t calo_track_drval[NTRACK][NCALO], tk2calo_dr_t calo_track_drval_cut[NTRACK][NCALO]) {
    const tk2calo_dr_t eDR2MAX = DR2MAX;
    for (int it = 0; it < NTRACK; ++it) {
        pt_t caloPtMin = track[it].hwPt - 2*(track[it].hwPtErr);
        if (caloPtMin < 0) caloPtMin = 0;
        for (int icalo = 0; icalo < NCALO; ++icalo) {
            calo_track_drval[it][icalo] = dr2_int_cap(track[it].hwEta, track[it].hwPhi, calo[icalo].hwEta, calo[icalo].hwPhi, eDR2MAX);
            if (calo[icalo].hwPt > caloPtMin) {
                calo_track_drval_cut[it][icalo] = calo_track_drval[it][icalo];
            } else {
                calo_track_drval_cut[it][icalo] = eDR2MAX;
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

template<int DR2MAX, int DR2MAXcut>
void tk2calo_drdptvals_cut(HadCaloObj calo[NCALO], TkObj track[NTRACK], tk2calo_dr_t calo_track_drval[NTRACK][NCALO], tk2calo_dq_t calo_track_drval_cut[NTRACK][NCALO]) {
    const tk2calo_dr_t eDR2MAX = DR2MAX;
    const tk2calo_dq_t eDR2MAXcut = DR2MAXcut;
    const tk2calo_dq_t eDQMAX  = 5*DR2MAXcut; // at most we're 2 sigma away in pt, so that's a factor 4
    // now, DR2MAX is 10 bits, so dptscale max is at most 10+8 bits = 18 bits
    for (int it = 0; it < NTRACK; ++it) {
        pt_t caloPtMin = track[it].hwPt - 2*(track[it].hwPtErr);
        ap_int<18> dptscale  = calc_dptscale<DR2MAXcut>(track[it].hwPtErr);
        if (caloPtMin < 0) caloPtMin = 0;
        for (int icalo = 0; icalo < NCALO; ++icalo) {
            int dr2 = dr2_int_cap(track[it].hwEta, track[it].hwPhi, calo[icalo].hwEta, calo[icalo].hwPhi, eDR2MAX);
            if (calo[icalo].hwPt > caloPtMin) {
                //calo_track_drval_cut[it][icalo] = dr2_dpt_int_cap(track[it].hwEta, track[it].hwPhi, calo[icalo].hwEta, calo[icalo].hwPhi, track[it].hwPt, calo[icalo].hwPt, dptscale, eDR2MAX, eDQMAX);
                calo_track_drval_cut[it][icalo] = dr2_plus_dpt_int_cap(dr2, track[it].hwPt, calo[icalo].hwPt, dptscale, eDR2MAXcut, eDQMAX);
                //if (calo_track_drval_cut[it][icalo] < eDQMAX) printf("HWO DQ(track %+7d %+7d  calo %3d) = %12d\n", int(track[it].hwEta), int(track[it].hwPhi), icalo, int(calo_track_drval_cut[it][icalo]));
            } else {
                calo_track_drval_cut[it][icalo] = eDQMAX;
            }
            calo_track_drval[it][icalo] = tk2calo_dr_t(dr2);
        }
    }
}

template<int DR2MAX, int NTK, int NCA, typename DR_T>
void pick_closest(DR_T calo_track_drval[NTK][NCA], ap_uint<NCA> calo_track_link_bit[NTK]) {
    const DR_T eDR2MAX = DR2MAX;
    for (int it = 0; it < NTK; ++it) {
        /*for (int icalo = 0; icalo < NCA; ++icalo) {
            DR_T mydr = calo_track_drval[it][icalo];
            bool link = (mydr != eDR2MAX);
            for (int j = 0; j < NCA; ++j) {
                if (icalo <= j) link = link && (calo_track_drval[it][j] >= mydr);
                else            link = link && (calo_track_drval[it][j] >  mydr);
            }
            calo_track_link_bit[it][icalo] = link;
        }*/
//causes loop unroll issue for large number of calo/tracks (happens at 25^3)
//calo_track_link_bit[it][icalo] is true if all preceeding calo are >= and all succeeding are >
        DR_T mydr = calo_track_drval[it][0];
        int index = 0;
        for (int icalo = 1; icalo < NCA; ++icalo) {
            if (mydr >= calo_track_drval[it][icalo]) {
                mydr = calo_track_drval[it][icalo];
                index = icalo;
            }
        }
        calo_track_link_bit[it] = 0;
        calo_track_link_bit[it][index] = (mydr >= eDR2MAX) ? 0 : 1;
    }
}

void tk2calo_link_dronly(HadCaloObj calo[NCALO], TkObj track[NTRACK], ap_uint<NCALO> calo_track_link_bit[NTRACK], tk2calo_dr_t drvals[NTRACK][NCALO]) {
    const int DR2MAX = PFALGO2_HGC_DR2MAX_TK_CALO;
    //tk2calo_dr_t drvals[NTRACK][NCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0
    tk2calo_dr_t drvals_cut[NTRACK][NCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals_cut complete dim=0

    //tk2calo_drvals<DR2MAX>(calo, track, drvals);
    tk2calo_drvals_cut<PFPUPPI_DR2MAX>(calo, track, drvals, drvals_cut);
    pick_closest<DR2MAX,NTRACK,NCALO,tk2calo_dr_t>(drvals_cut, calo_track_link_bit);
}
void tk2calo_link_drdpt(HadCaloObj calo[NCALO], TkObj track[NTRACK], ap_uint<NCALO> calo_track_link_bit[NTRACK], tk2calo_dr_t drvals[NTRACK][NCALO]) {
    const int DR2MAX = PFALGO2_HGC_DR2MAX_TK_CALO;
    const int DQMAX = 5*DR2MAX;
    //tk2calo_dq_t drvals[NTRACK][NCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0
    tk2calo_dq_t drvals_cut[NTRACK][NCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals_cut complete dim=0

    //tk2calo_drdptvals<DR2MAX>(calo, track, drvals);
    tk2calo_drdptvals_cut<PFPUPPI_DR2MAX,DR2MAX>(calo, track, drvals, drvals_cut);
    pick_closest<DQMAX,NTRACK,NCALO,tk2calo_dq_t>(drvals_cut, calo_track_link_bit);
}

void tk2calo_em(HadCaloObj calo[NCALO], ap_uint<NCALO> calo_track_link_bit[NTRACK], bool isEM[NTRACK]) {
    for (int it = 0; it < NTRACK; it++) {
        for (int ic = 0; ic < NCALO; ic++) {
            if (calo_track_link_bit[it][ic] == 1) isEM[it] = calo[ic].hwIsEM;
        }
    }
}

void tk2calo_tkerr2(TkObj track[NTRACK], int tkerr2[NTRACK]) {
    for (int it = 0; it < NTRACK; ++it) {
        tkerr2[it] = (track[it].hwPtErr * track[it].hwPtErr);
    }
}
void tk2calo_sumtk(TkObj track[NTRACK], bool isMu[NTRACK], int tkerr2[NTRACK], ap_uint<NCALO> calo_track_link_bit[NTRACK], pt_t sumtk[NCALO], int sumtkerr2[NCALO]) {
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

void tk2calo_tkalgo(TkObj track[NTRACK], bool isMu[NTRACK], bool isEM[NTRACK], ap_uint<NCALO> calo_track_link_bit[NTRACK], PFChargedObj pfout[NTRACK]) {
    const pt_t TKPT_MAX_LOOSE = PFALGO2_HGC_TK_MAXINVPT_LOOSE; // 20 * PT_SCALE;
    const pt_t TKPT_MAX_TIGHT = PFALGO2_HGC_TK_MAXINVPT_TIGHT; // 20 * PT_SCALE;
    for (int it = 0; it < NTRACK; ++it) {
        bool goodByPt = track[it].hwPt < (track[it].hwTightQuality ? TKPT_MAX_TIGHT : TKPT_MAX_LOOSE);
        bool good = isMu[it] || goodByPt || calo_track_link_bit[it].or_reduce();
        if (good) {
            pfout[it].hwPt  = track[it].hwPt;
            pfout[it].hwEta = track[it].hwEta;
            pfout[it].hwPhi = track[it].hwPhi;
            pfout[it].hwId  = isMu[it] ? PID_Muon : (isEM[it] ? PID_Electron : PID_Charged);
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


void tk2calo_caloalgo(HadCaloObj calo[NCALO], pt_t sumtk[NCALO], int sumtkerr2[NCALO], PFNeutralObj pfout[NCALO], PFNeutralObj phoout[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t calopt = 0;
        pt_t emcalopt = 0;
        pt_t ptdiff = calo[icalo].hwPt - sumtk[icalo];
        if (ptdiff > 0 && (ptdiff*ptdiff) > (sumtkerr2[icalo] + (sumtkerr2[icalo] >> 1))) {
            if (calo[icalo].hwEmPt > 1) {
                emcalopt = (calo[icalo].hwEmPt < ptdiff ? calo[icalo].hwEmPt : ptdiff);
                ptdiff -= emcalopt;
            }
            if (ptdiff > 1) {
                calopt = ptdiff;
            }
        }
        pfout[icalo].hwPt  = calopt;
        pfout[icalo].hwEta = calopt > 0 ? calo[icalo].hwEta : etaphi_t(0);
        pfout[icalo].hwPhi = calopt > 0 ? calo[icalo].hwPhi : etaphi_t(0);
        pfout[icalo].hwId  = calopt > 0 ? PID_Neutral : 0;
        phoout[icalo].hwPt  = emcalopt;
        phoout[icalo].hwEta = emcalopt > 0 ? calo[icalo].hwEta : etaphi_t(0);
        phoout[icalo].hwPhi = emcalopt > 0 ? calo[icalo].hwPhi : etaphi_t(0);
        phoout[icalo].hwId  = emcalopt > 0 ? PID_Photon : 0;
    }
}

//-------------------------------------------------------
// TK-MU Algos
//-------------------------------------------------------

void spfph_mu2trk_dptvals(MuObj mu[NMU], TkObj track[NTRACK], pt_t mu_track_dptval[NMU][NTRACK]) {
    const ap_uint<12> DR2MAX = PFALGO2_HGC_DR2MAX_TK_MU;
    //const etaphi_t DR2MAX = PFALGO2_HGC_DR2MAX_TK_MU;
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

void pfalgo2_hgc(HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU], tk2em_dr_t drvals_tk2em[NTRACK][NPHOTON], tk2calo_dr_t drvals_tk2calo[NTRACK][NSELCALO]) {
    
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete    
    #pragma HLS ARRAY_PARTITION variable=outch complete
    #pragma HLS ARRAY_PARTITION variable=outpho complete
    #pragma HLS ARRAY_PARTITION variable=outne complete
    #pragma HLS ARRAY_PARTITION variable=outmu complete

    #pragma HLS ARRAY_PARTITION variable=drvals_tk2em   complete dim=0
    #pragma HLS ARRAY_PARTITION variable=drvals_tk2calo complete dim=0

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
    // TK-HAD Linking
    ap_uint<NCALO> calo_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=calo_track_link_bit complete
    tk2calo_dr_t drvals_tk2calo_unsort[NTRACK][NCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals_tk2calo_unsort complete dim=0
    bool isEM[NTRACK];
    //initialize
    for (int it = 0; it < NTRACK; it++) {isEM[it] = 0;}
    tk2calo_link_drdpt(hadcalo, track, calo_track_link_bit, drvals_tk2calo_unsort);
    //tk2calo_link_dronly(hadcalo, track, calo_track_link_bit, drvals_tk2calo_unsort);

    #pragma HLS ARRAY_PARTITION variable=isEM complete
    tk2calo_em(hadcalo, calo_track_link_bit, isEM);

    int tkerr2[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=tkerr2 complete
    tk2calo_tkerr2(track, tkerr2);

    pt_t sumtk[NCALO]; int sumtkerr2[NCALO];
    #pragma HLS ARRAY_PARTITION variable=sumtk complete
    #pragma HLS ARRAY_PARTITION variable=sumtkerr2 complete

    tk2calo_tkalgo(track, isMu, isEM, calo_track_link_bit, outch);
    tk2calo_sumtk(track, isMu, tkerr2, calo_track_link_bit, sumtk, sumtkerr2);

    PFNeutralObj outne_all[NCALO];
    #pragma HLS ARRAY_PARTITION variable=outne_all complete
    tk2calo_caloalgo(hadcalo, sumtk, sumtkerr2, outne_all, outpho);
    //ptsort_hwopt<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);
    ap_uint<6> calind[NSELCALO];
    #pragma HLS ARRAY_PARTITION variable=calind complete
    ptsort_hwopt_ind<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne, calind);
    for (int ic=0; ic<NSELCALO; ic++) {
        #pragma HLS LOOP UNROLL
        for (int it=0; it<NTRACK; it++) {
            #pragma HLS LOOP UNROLL
            drvals_tk2calo[it][ic] = drvals_tk2calo_unsort[it][calind[ic]];
        }
    }
    for (int ic=0; ic<NCALO; ic++) {
        #pragma HLS LOOP UNROLL
        for (int it=0; it<NTRACK; it++) {
            if (!isMu[it]) drvals_tk2em[it][ic] = drvals_tk2calo_unsort[it][ic];
            else drvals_tk2em[it][ic] = PFALGO2_HGC_DR2MAX_TK_CALO;
        }
    }

}


void mp7wrapped_pack_in_hgc(HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], MP7DataWord data[MP7_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete
    // pack inputs
    assert(2*NTRACK + 2*NCALO + 2*NMU <= MP7_NCHANN);
    #define HADOFFS 0
    #define TKOFFS 2*NCALO+HADOFFS
    #define MUOFFS 2*NTRACK+TKOFFS
    mp7_pack<NCALO,HADOFFS>(hadcalo, data);
    mp7_pack<NTRACK,TKOFFS>(track, data);
    mp7_pack<NMU,MUOFFS>(mu, data);
}

void mp7wrapped_unpack_in_hgc(MP7DataWord data[MP7_NCHANN], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete
    // unpack inputs
    assert(2*NTRACK + 2*NCALO + 2*NMU <= MP7_NCHANN);
    #define HADOFFS 0
    #define TKOFFS 2*NCALO+HADOFFS
    #define MUOFFS 2*NTRACK+TKOFFS
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
void mp7wrapped_pack_out_necomb( PFChargedObj pfch[NTRACK], PFNeutralObj pfne_all[NNEUTRALS], PFChargedObj pfmu[NMU], MP7DataWord data[MP7_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfne_all complete
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
        data[2*i+0+PHOOFFS] = ( pfne_all[i].hwPtPuppi, pfne_all[i].hwPt );
        data[2*i+1+PHOOFFS] = ( pfne_all[i].hwId, pfne_all[i].hwPhi, pfne_all[i].hwEta );
    }
    for (unsigned int i = 0; i < NSELCALO; ++i) {
        data[2*i+0+NHOFFS] = ( pfne_all[i+NPHOTON].hwPtPuppi, pfne_all[i+NPHOTON].hwPt );
        data[2*i+1+NHOFFS] = ( pfne_all[i+NPHOTON].hwId, pfne_all[i+NPHOTON].hwPhi, pfne_all[i+NPHOTON].hwEta );
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
void mp7wrapped_unpack_out_necomb( MP7DataWord data[MP7_NCHANN], PFChargedObj pfch[NTRACK], PFNeutralObj pfpho[NPHOTON], PFNeutralObj pfne[NSELCALO], PFChargedObj pfmu[NMU]) {
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
        pfpho[i].hwPtPuppi  = data[2*i+0+PHOOFFS](31,16);
        pfpho[i].hwEta = data[2*i+1+PHOOFFS](9, 0);
        pfpho[i].hwPhi = data[2*i+1+PHOOFFS](19,10);
        pfpho[i].hwId  = data[2*i+1+PHOOFFS](22,20);
    }
    for (unsigned int i = 0; i < NSELCALO; ++i) {
        pfne[i].hwPt  = data[2*i+0+NHOFFS](15, 0);
        pfne[i].hwPtPuppi  = data[2*i+0+NHOFFS](31, 16);
        pfne[i].hwEta = data[2*i+1+NHOFFS](9, 0);
        pfne[i].hwPhi = data[2*i+1+NHOFFS](19,10);
        pfne[i].hwId  = data[2*i+1+NHOFFS](22,20);
    }
    for (unsigned int i = 0; i < NMU; ++i) {
        pfmu[i].hwPt  = data[2*i+0+PFMUOFFS](15, 0);
        pfmu[i].hwId  = data[2*i+0+PFMUOFFS](18,16);
        pfmu[i].hwEta = data[2*i+1+PFMUOFFS](9, 0);
        pfmu[i].hwPhi = data[2*i+1+PFMUOFFS](19,10);
        pfmu[i].hwZ0  = data[2*i+1+PFMUOFFS](29,20);
    }

}

void mp7wrapped_pfalgo2_hgc(MP7DataWord input[MP7_NCHANN], MP7DataWord output[MP7_NCHANN], z0_t Z0) {
    
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=output complete
    #pragma HLS INTERFACE ap_none port=output

    #pragma HLS pipeline II=HLS_pipeline_II

    HadCaloObj hadcalo[NCALO]; TkObj track[NTRACK]; MuObj mu[NMU];
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete

    PFChargedObj pfch[NTRACK]; PFNeutralObj pfpho[NPHOTON]; PFNeutralObj pfne[NSELCALO]; PFChargedObj pfmu[NMU];
    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfpho complete
    #pragma HLS ARRAY_PARTITION variable=pfne complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete

    mp7wrapped_unpack_in_hgc(input, hadcalo, track, mu);
    tk2em_dr_t drvals_tk2em[NTRACK][NPHOTON];
    tk2calo_dr_t drvals_tk2calo[NTRACK][NSELCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals_tk2em complete dim=0
    #pragma HLS ARRAY_PARTITION variable=drvals_tk2calo complete dim=0
    pfalgo2_hgc(hadcalo, track, mu, pfch, pfpho, pfne, pfmu, drvals_tk2em, drvals_tk2calo);
    //concat drvals, ne
    tk2calo_dr_t drvals[NTRACK][NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0
    PFNeutralObj pfne_all[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=pfne_all complete dim=0
    for (int i=0; i<NPHOTON; i++) {
        #pragma HLS LOOP UNROLL
        for (int j=0; j<NTRACK; j++) {
            #pragma HLS LOOP UNROLL
            drvals[j][i] = pfpho[i].hwPt == 0 ? tk2calo_dr_t(PFPUPPI_DR2MAX) : tk2calo_dr_t(drvals_tk2em[j][i]);
        }
        pfne_all[i] = pfpho[i];
    }
    for (int i=0; i<NSELCALO; i++) {
        #pragma HLS LOOP UNROLL
        for (int j=0; j<NTRACK; j++) {
            #pragma HLS LOOP UNROLL
            drvals[j][i+NPHOTON] = pfne[i].hwPt == 0 ? tk2calo_dr_t(PFPUPPI_DR2MAX) : drvals_tk2calo[j][i];
        }
        pfne_all[i+NPHOTON] = pfne[i];
    }
    simple_puppi_hw(pfch, pfne_all, drvals, Z0);
    mp7wrapped_pack_out_necomb(pfch, pfne_all, pfmu, output);
    //mp7wrapped_pack_out(pfch, pfpho, pfne, pfmu, output);

}

