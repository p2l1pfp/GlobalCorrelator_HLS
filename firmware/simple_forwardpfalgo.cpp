#include "simple_forwardpfalgo.h"
#include "../puppi/firmware/simple_puppi_forward.h"
#include "mp7pf_encoding.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

//typedef ap_uint<10> em2calo_dr_t;

bool match_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, etaphi_t boxSize) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    return ((deta<boxSize && deta>-boxSize) && (dphi<boxSize && dphi>-boxSize));
}
int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}
/*template<int NB>
ap_uint<NB> dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) {
    etaphi_t deta = (eta1-eta2);
    //#pragma HLS RESOURCE variable=deta core=AddSub_DSP
    etaphi_t dphi = (phi1-phi2);
    //#pragma HLS RESOURCE variable=dphi core=AddSub_DSP
    int deta2 = deta*deta;
    //#pragma HLS RESOURCE variable=deta2 core=Mul
    int dphi2 = dphi*dphi;
    //#pragma HLS RESOURCE variable=dphi2 core=Mul
    int dr2 = deta2 + dphi2;
    //#pragma HLS RESOURCE variable=dr2 core=AddSub_DSP
    //int dr2 = deta*deta + dphi*dphi;
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
}*/
//new
template<int NB>
ap_uint<NB> dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) {
    //hardcode for etaphi size
    int tmpe = eta1-eta2;
    ap_uint<NB> deta = (tmpe > 0 ? tmpe : -tmpe);
    //#pragma HLS RESOURCE variable=deta core=AddSub_DSP
    int tmpp = phi1-phi2;
    ap_uint<NB> dphi = (tmpp > 0 ? tmpp : -tmpp);
    //#pragma HLS RESOURCE variable=dphi core=AddSub_DSP
    int dr2 = max;
    if ((deta >> (NB/2))==0 && (dphi >> (NB/2))==0) {
        ap_uint<NB> deta2 = deta*deta;
        //#pragma HLS RESOURCE variable=deta2 core=Mul
        ap_uint<NB> dphi2 = dphi*dphi;
        //#pragma HLS RESOURCE variable=dphi2 core=Mul
        dr2 = deta2 + dphi2;
        //#pragma HLS RESOURCE variable=dr2 core=AddSub_DSP
    }
    //int dr2 = deta*deta + dphi*dphi;
    return (dr2 < int(max) ? ap_uint<NB>(dr2) : max);
}
//
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
void init_dr2max_times_pterr2_inv(int vals[512]) {
    for (int i = 0; i < 512; ++i) {
    	int tmp = (DR2MAX<<8)/(i?i*i:1), int18_max = (1<<17)-1;
        vals[i] = (tmp > int18_max ? int18_max : tmp);
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

template<int DR2MAX>
void em2calo_drvals_cut(HadCaloObj hadcalo[NCALO], EmCaloObj emcalo[NEMCALO], em2calo_dr_t hadcalo_emcalo_drval[NEMCALO][NCALO], em2calo_dr_t hadcalo_emcalo_drval_cut[NEMCALO][NCALO]) {
    const em2calo_dr_t eDR2MAX = DR2MAX;
    for (int it = 0; it < NEMCALO; ++it) {
        pt_t hadcaloEmPtMin = (emcalo[it].hwPt >> 1);
        for (int icalo = 0; icalo < NCALO; ++icalo) {
            hadcalo_emcalo_drval[it][icalo] = dr2_int_cap(emcalo[it].hwEta, emcalo[it].hwPhi, hadcalo[icalo].hwEta, hadcalo[icalo].hwPhi, eDR2MAX);
            if (hadcalo[icalo].hwEmPt > hadcaloEmPtMin) {
                hadcalo_emcalo_drval_cut[it][icalo] = hadcalo_emcalo_drval[it][icalo];
            } else {
                hadcalo_emcalo_drval_cut[it][icalo] = eDR2MAX;
            }
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

void em2calo_link(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], ap_uint<NCALO> em_calo_link_bit[NCALO], em2calo_dr_t drvals[NEMCALO][NCALO]) {
    const int DR2MAX = PFALGO3_DR2MAX_EM_CALO;
    em2calo_dr_t drvals_cut[NEMCALO][NCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0
    #pragma HLS ARRAY_PARTITION variable=drvals_cut complete dim=0

    em2calo_drvals_cut<PFPUPPI_DR2MAX>(hadcalo, emcalo, drvals, drvals_cut);
    pick_closest<DR2MAX,NEMCALO,NCALO,em2calo_dr_t>(drvals_cut, em_calo_link_bit);
}

void em2calo_sumem(EmCaloObj emcalo[NEMCALO], ap_uint<NCALO> em_had_link_bit[NEMCALO], pt_t sumem[NEMCALO], bool keepcalo[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t sum = 0; bool keep = false;
        for (int iem = 0; iem < NEMCALO; ++iem) {
            if (em_had_link_bit[iem][icalo]) { 
                sum += emcalo[iem].hwPt; 
            }
        }
        sumem[icalo] = sum;
        keepcalo[icalo] = keep;
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

void em_photons(EmCaloObj calo[NEMCALO], PFNeutralObj pfout[NPHOTON]) {
    for (int icalo = 0; icalo < NEMCALO; ++icalo) {
        pt_t phopt = calo[icalo].hwPt;
        pfout[icalo].hwPt  = phopt;
        pfout[icalo].hwEta = phopt ? calo[icalo].hwEta : etaphi_t(0);
        pfout[icalo].hwPhi = phopt ? calo[icalo].hwPhi : etaphi_t(0);
        pfout[icalo].hwId  = phopt ? PID_Photon : 0;
    }
}

void had_caloalgo(HadCaloObj calo[NCALO], PFNeutralObj pfout[NCALO]) {
    for (int icalo = 0; icalo < NCALO; ++icalo) {
        pt_t calopt = calo[icalo].hwPt;
        pfout[icalo].hwPt  = calopt;
        pfout[icalo].hwEta = calopt ? calo[icalo].hwEta : etaphi_t(0);
        pfout[icalo].hwPhi = calopt ? calo[icalo].hwPhi : etaphi_t(0);
        pfout[icalo].hwId  = calopt ? PID_Neutral : 0;
    }
}

void spfph_mualgo_forward(MuObj mu[NMU], PFChargedObj pfmuout[NMU]) {

    for (int im = 0; im < NMU; ++im) {
        if (mu[im].hwPt > 0) {
            pfmuout[im].hwPt  = mu[im].hwPt;
            pfmuout[im].hwEta = mu[im].hwEta;
            pfmuout[im].hwPhi = mu[im].hwPhi;
            pfmuout[im].hwId  = PID_Muon;
            pfmuout[im].hwZ0  = 0;
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

void pfalgo3_forward(EmCaloObj calo[NEMCALO], HadCaloObj hadcalo[NCALO], MuObj mu[NMU], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU], em2calo_dr_t drvals_em2calo[NPHOTON][NSELCALO]) {
//void pfalgo3_forward(EmCaloObj calo[NEMCALO], HadCaloObj hadcalo[NCALO], MuObj mu[NMU], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) {
    
    #pragma HLS ARRAY_PARTITION variable=calo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=mu complete    
    #pragma HLS ARRAY_PARTITION variable=outpho complete
    #pragma HLS ARRAY_PARTITION variable=outne complete
    #pragma HLS ARRAY_PARTITION variable=outmu complete

    #pragma HLS ARRAY_PARTITION variable=drvals_em2calo complete dim=0

    #pragma HLS pipeline II=HLS_pipeline_II

    // ---------------------------------------------------------------
    // MU
    spfph_mualgo_forward(mu, outmu);

    // ---------------------------------------------------------------
    // EM
    pt_t photonPt[NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=photonPt complete

    em_photons(calo, outpho);


    em2calo_dr_t drvals_unsort[NEMCALO][NCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals_unsort complete dim=0

    ap_uint<NCALO> em_calo_link_bit[NEMCALO];
    #pragma HLS ARRAY_PARTITION variable=em_calo_link_bit complete
    em2calo_link(calo, hadcalo, em_calo_link_bit, drvals_unsort);

    bool keepcalo[NCALO];
    pt_t sumem[NCALO]; 
    #pragma HLS ARRAY_PARTITION variable=sumem complete
    em2calo_sumem(calo, em_calo_link_bit, sumem, keepcalo);

    HadCaloObj hadcalo_sub[NCALO];
    #pragma HLS ARRAY_PARTITION variable=hadcalo_sub complete

    em2calo_sub(hadcalo, sumem, keepcalo, hadcalo_sub);

    // ---------------------------------------------------------------
    // HAD
    PFNeutralObj outne_all[NCALO];
    #pragma HLS ARRAY_PARTITION variable=outne_all complete
    had_caloalgo(hadcalo_sub, outne_all);
    ptsort_hwopt<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);

    ap_uint<6> calind[NSELCALO];
    #pragma HLS ARRAY_PARTITION variable=calind complete
    ptsort_hwopt_ind<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne, calind);
    for (int ic=0; ic<NSELCALO; ic++) {
        #pragma HLS LOOP UNROLL
        for (int it=0; it<NPHOTON; it++) {
            #pragma HLS LOOP UNROLL
            drvals_em2calo[it][ic] = drvals_unsort[it][calind[ic]];
        }
    }

}


void mp7wrapped_pack_in(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], MuObj mu[NMU], MP7DataWord data[MP7_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=emcalo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=mu complete
    // pack inputs
    assert(2*NEMCALO + 2*NCALO + 2*NMU <= MP7_NCHANN);
    #define HADOFFS 2*NEMCALO
    #define MUOFFS 2*NCALO+HADOFFS
    mp7_pack<NEMCALO,0>(emcalo, data);
    mp7_pack<NCALO,HADOFFS>(hadcalo, data);
    mp7_pack<NMU,MUOFFS>(mu, data);
}

void mp7wrapped_unpack_in(MP7DataWord data[MP7_NCHANN], EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], MuObj mu[NMU]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=emcalo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=mu complete
    // unpack inputs
    assert(2*NEMCALO + 2*NCALO + 2*NMU <= MP7_NCHANN);
    #define HADOFFS 2*NEMCALO
    #define MUOFFS 2*NCALO+HADOFFS
    mp7_unpack<NEMCALO,0>(data, emcalo);
    mp7_unpack<NCALO,HADOFFS>(data, hadcalo);
    mp7_unpack<NMU,MUOFFS>(data, mu);
}

void mp7wrapped_pack_out( PFNeutralObj pfpho[NPHOTON], PFNeutralObj pfne[NSELCALO], PFChargedObj pfmu[NMU], MP7DataWord data[MP7_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfpho complete
    #pragma HLS ARRAY_PARTITION variable=pfne complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete
    // pack outputs
    assert(2*NPHOTON + 2*NSELCALO + 2*NMU <= MP7_NCHANN);
    #define PHOOFFS 0
    #define NHOFFS 2*NPHOTON+PHOOFFS
    #define PFMUOFFS 2*NSELCALO+NHOFFS
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
void mp7wrapped_pack_out_necomb( PFNeutralObj pfne_all[NNEUTRALS], PFChargedObj pfmu[NMU], MP7DataWord data[MP7_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfne_all complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete
    // pack outputs
    assert(2*NPHOTON + 2*NSELCALO + 2*NMU <= MP7_NCHANN);
    #define PHOOFFS 0
    #define NHOFFS 2*NPHOTON+PHOOFFS
    #define PFMUOFFS 2*NSELCALO+NHOFFS
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
void mp7wrapped_unpack_out( MP7DataWord data[MP7_NCHANN], PFNeutralObj pfpho[NPHOTON], PFNeutralObj pfne[NSELCALO], PFChargedObj pfmu[NMU]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfpho complete
    #pragma HLS ARRAY_PARTITION variable=pfne complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete
    // unpack outputs
    assert(2*NPHOTON + 2*NSELCALO + 2*NMU <= MP7_NCHANN);
    #define PHOOFFS 0
    #define NHOFFS 2*NPHOTON+PHOOFFS
    #define PFMUOFFS 2*NSELCALO+NHOFFS
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
void mp7wrapped_unpack_out_necomb( MP7DataWord data[MP7_NCHANN], PFNeutralObj pfpho[NPHOTON], PFNeutralObj pfne[NSELCALO], PFChargedObj pfmu[NMU]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfpho complete
    #pragma HLS ARRAY_PARTITION variable=pfne complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete
    // unpack outputs
    assert(2*NPHOTON + 2*NSELCALO + 2*NMU <= MP7_NCHANN);
    #define PHOOFFS 0
    #define NHOFFS 2*NPHOTON+PHOOFFS
    #define PFMUOFFS 2*NSELCALO+NHOFFS
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

void mp7wrapped_pfalgo3_forward(MP7DataWord input[MP7_NCHANN], MP7DataWord output[MP7_NCHANN]) {
    
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=output complete
    #pragma HLS INTERFACE ap_none port=output

    #pragma HLS pipeline II=HLS_pipeline_II

    EmCaloObj emcalo[NEMCALO]; HadCaloObj hadcalo[NCALO]; MuObj mu[NMU];
    #pragma HLS ARRAY_PARTITION variable=emcalo complete
    #pragma HLS ARRAY_PARTITION variable=hadcalo complete
    #pragma HLS ARRAY_PARTITION variable=mu complete

    PFNeutralObj pfpho[NPHOTON]; PFNeutralObj pfne[NSELCALO]; PFChargedObj pfmu[NMU];
    #pragma HLS ARRAY_PARTITION variable=pfpho complete
    #pragma HLS ARRAY_PARTITION variable=pfne complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete

    mp7wrapped_unpack_in(input, emcalo, hadcalo, mu);
    em2calo_dr_t drvals_em2calo[NPHOTON][NSELCALO];
    #pragma HLS ARRAY_PARTITION variable=drvals_em2calo complete dim=0
    pfalgo3_forward(emcalo, hadcalo, mu, pfpho, pfne, pfmu, drvals_em2calo);
    //pfalgo3_forward(emcalo, hadcalo, mu, pfpho, pfne, pfmu);
    //concat drvals, ne
    PFNeutralObj pfne_all[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=pfne_all complete dim=0
    pt_t ptpuppi[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=ptpuppi complete
    for (int i=0; i<NPHOTON; i++) {
        #pragma HLS LOOP UNROLL
        pfne_all[i] = pfpho[i];
        ptpuppi[i] = 0;
    }
    for (int i=0; i<NSELCALO; i++) {
        #pragma HLS LOOP UNROLL
        pfne_all[i+NPHOTON] = pfne[i];
        ptpuppi[i+NPHOTON] = 0;
    }
    simple_puppi_forward_hw(pfne_all, ptpuppi, drvals_em2calo);
    for (int in = 0; in < NNEUTRALS; in++) {
        #pragma HLS LOOP UNROLL
        pfne_all[in].hwPtPuppi = ptpuppi[in];
    }
    mp7wrapped_pack_out_necomb(pfne_all, pfmu, output);
    //mp7wrapped_pack_out(pfpho, pfne, pfmu, output);

}

