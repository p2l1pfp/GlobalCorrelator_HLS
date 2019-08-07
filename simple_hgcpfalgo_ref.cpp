#include "firmware/data.h"
#include "firmware/simple_hgcpfalgo.h"
#include "DiscretePFInputs.h"
#include "utils/Firmware2DiscretePF.h"
#include <cmath>
#include <algorithm>

bool g_debug_ = 0;

void pfalgo2_hgc_ref_set_debug(bool debug) { g_debug_ = debug; }

template <typename T> int sqr(const T & t) { return t*t; }

template<int NCAL, int DR2MAX, bool doPtMin, typename CO_t>
int best_match_ref(CO_t calo[NCAL], const TkObj & track) {
    pt_t caloPtMin = track.hwPt - 2*(track.hwPtErr);
    if (caloPtMin < 0) caloPtMin = 0;
    int  drmin = DR2MAX, ibest = -1;
    for (int ic = 0; ic < NCAL; ++ic) {
            if (doPtMin && calo[ic].hwPt <= caloPtMin) continue;
            int dr = dr2_int(track.hwEta, track.hwPhi, calo[ic].hwEta, calo[ic].hwPhi);
            if (dr < drmin) { drmin = dr; ibest = ic; }
    }
    return ibest;
}
template<int NCAL, int DR2MAX>
int best_match_ref(HadCaloObj calo[NCAL], const EmCaloObj & em) {
    pt_t emPtMin = em.hwPt >> 1;
    int  drmin = DR2MAX, ibest = -1;
    for (int ic = 0; ic < NCAL; ++ic) {
            if (calo[ic].hwEmPt <= emPtMin) continue;
            int dr = dr2_int(em.hwEta, em.hwPhi, calo[ic].hwEta, calo[ic].hwPhi);
            if (dr < drmin) { drmin = dr; ibest = ic; }
    }
    return ibest;
}

template<int NCAL, int DR2MAX, typename CO_t>
int best_match_with_pt_ref(CO_t calo[NCAL], const TkObj & track) {
    pt_t caloPtMin = track.hwPt - 2*(track.hwPtErr);
    if (caloPtMin < 0) caloPtMin = 0;
    int dptscale = (DR2MAX<<8)/std::max<int>(1,sqr(track.hwPtErr));
    int drmin = 0, ibest = -1;
    for (int ic = 0; ic < NCAL; ++ic) {
            if (calo[ic].hwPt <= caloPtMin) continue;
            int dr = dr2_int(track.hwEta, track.hwPhi, calo[ic].hwEta, calo[ic].hwPhi);
            if (dr >= DR2MAX) continue;
            dr += (( sqr(std::max<int>(track.hwPt-calo[ic].hwPt,0))*dptscale ) >> 8);
            if (g_debug_) printf("REF DQ(track %+7d %+7d  calo %3d) = %12d\n", int(track.hwEta), int(track.hwPhi), ic, dr);
            if (ibest == -1 || dr <= drmin) { drmin = dr; ibest = ic; }
    }
    return ibest;
}


template<int NCAL, int DR2MAX, bool doPtMin, typename CO_t>
void link_ref(CO_t calo[NCAL], TkObj track[NTRACK], ap_uint<NCAL> calo_track_link_bit[NTRACK]) {
    for (int it = 0; it < NTRACK; ++it) {
        int ibest = best_match_ref<NCALO,DR2MAX,doPtMin,CO_t>(calo, track[it]);
        calo_track_link_bit[it] = 0;
        if (ibest != -1) calo_track_link_bit[it][ibest] = 1;
    }
}

template<typename T, int NIn, int NOut>
void ptsort_ref(T in[NIn], T out[NOut]) {
    for (int iout = 0; iout < NOut; ++iout) {
        out[iout].hwPt = 0;
    }
    for (int it = 0; it < NIn; ++it) {
        for (int iout = 0; iout < NOut; ++iout) {
            if (in[it].hwPt >= out[iout].hwPt) {
                for (int i2 = NOut-1; i2 > iout; --i2) {
                    out[i2] = out[i2-1];
                }
                out[iout] = in[it];
                break;
            }
        }
    }
}

void pfalgo2_hgc_ref(HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) {

    if (g_debug_) {
#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_MORE
        for (int i = 0; i < NTRACK; ++i) { if (track[i].hwPt == 0) continue;
            l1tpf_int::PropagatedTrack tk; fw2dpf::convert(track[i], tk); 
            printf("FW  \t track %3d: pt %8d [ %7.2f ]  calo eta %+7d [ %+5.2f ]  calo phi %+7d [ %+5.2f ]  calo ptErr %6d [ %7.2f ] \n", 
                                i, tk.hwPt, tk.floatPt(), tk.hwEta, tk.floatEta(), tk.hwPhi, tk.floatPhi(), tk.hwCaloPtErr, tk.floatCaloPtErr());
        }
        for (int i = 0; i < NCALO; ++i) { if (hadcalo[i].hwPt == 0) continue;
            l1tpf_int::CaloCluster calo; fw2dpf::convert(hadcalo[i], calo); 
            printf("FW  \t calo  %3d: pt %8d [ %7.2f ]  calo eta %+7d [ %+5.2f ]  calo phi %+7d [ %+5.2f ]  calo emPt %7d [ %7.2f ]   isEM %d \n", 
                                i, calo.hwPt, calo.floatPt(), calo.hwEta, calo.floatEta(), calo.hwPhi, calo.floatPhi(), calo.hwEmPt, calo.floatEmPt(), calo.isEM);
        } 
        for (int i = 0; i < NMU; ++i) { if (mu[i].hwPt == 0) continue;
            l1tpf_int::Muon muon; fw2dpf::convert(mu[i], muon); 
            printf("FW  \t muon  %3d: pt %8d [ %7.2f ]  muon eta %+7d [ %+5.2f ]  muon phi %+7d [ %+5.2f ]   \n", 
                                i, muon.hwPt, muon.floatPt(), muon.hwEta, muon.floatEta(), muon.hwPhi, muon.floatPhi());
        } 
#endif
    }

    // constants
    const pt_t     TKPT_MAX_LOOSE = PFALGO2_HGC_TK_MAXINVPT_LOOSE; // 20 * PT_SCALE;
    const pt_t     TKPT_MAX_TIGHT = PFALGO2_HGC_TK_MAXINVPT_TIGHT; // 20 * PT_SCALE;
    const int      DR2MAX   = PFALGO2_HGC_DR2MAX_TK_CALO;
    const int      DR2MAX_TM = PFALGO2_HGC_DR2MAX_TK_MU;

    ////////////////////////////////////////////////////
    // TK-MU Linking

    // initialize good track bit
    // bool mu_good[NMU];
    // for (int im = 0; im < NMU; ++im) { mu_good[im] = (mu[im].hwPt < TKPT_MAX); }

    // initialize output
    for (int ipf = 0; ipf < NMU; ++ipf) { outmu[ipf].hwPt = 0; outmu[ipf].hwEta = 0; outmu[ipf].hwPhi = 0; outmu[ipf].hwId  = 0; outmu[ipf].hwZ0  = 0; }

    bool isMu[NTRACK];
    for (int it = 0; it < NTRACK; ++it) { isMu[it] = 0; } // initialize
    // for each muon, find the closest track
    for (int im = 0; im < NMU; ++im) {
        if (mu[im].hwPt > 0) {
            int ibest = -1;
            int dptmin = mu[im].hwPt >> 1;
            for (int it = 0; it < NTRACK; ++it) {
                int dr = dr2_int(mu[im].hwEta, mu[im].hwPhi, track[it].hwEta, track[it].hwPhi);
                //printf("deltaR2(mu %d float pt %5.1f, tk %2d float pt %5.1f) = int %d  (float deltaR = %.3f); int cut at %d\n", im, 0.25*int(mu[im].hwPt), it, 0.25*int(track[it].hwPt), dr, std::sqrt(float(dr))/229.2, PFALGO2_HGC_DR2MAX_TK_MU);
                if (dr < DR2MAX_TM) { 
                    int dpt = std::abs(int(track[it].hwPt - mu[im].hwPt));
                    if (dpt < dptmin) {
                        dptmin = dpt; ibest = it; 
                    }
                }
            }
            if (ibest != -1) {
                outmu[im].hwPt = track[ibest].hwPt;
                outmu[im].hwEta = track[ibest].hwEta;
                outmu[im].hwPhi = track[ibest].hwPhi;
                outmu[im].hwId  = PID_Muon;
                outmu[im].hwZ0 = track[ibest].hwZ0;      
                isMu[ibest] = 1;
                if (g_debug_) printf("FW  \t muon %3d linked to track %3d \n", im, ibest);
            } else {
                if (g_debug_) printf("FW  \t muon %3d not linked to any track\n", im);
            }
        }
    }


    ////////////////////////////////////////////////////
    // TK-HAD Linking

    // initialize sum track pt
    pt_t calo_sumtk[NCALO], calo_pt[NCALO], emcalo_pt[NCALO];
    int  calo_sumtkErr2[NCALO];
    for (int ic = 0; ic < NCALO; ++ic) { calo_sumtk[ic] = 0;  calo_sumtkErr2[ic] = 0;}

    // initialize good track bit
    bool track_good[NTRACK];
    for (int it = 0; it < NTRACK; ++it) { 
        track_good[it] = (track[it].hwPt < (track[it].hwTightQuality ? TKPT_MAX_TIGHT : TKPT_MAX_LOOSE) || isMu[it]); 
    }

    // initialize output
    for (int ipf = 0; ipf < NTRACK; ++ipf) { outch[ipf].hwPt = 0; outch[ipf].hwEta = 0; outch[ipf].hwPhi = 0; outch[ipf].hwId = 0; outch[ipf].hwZ0 = 0; }
    for (int ipf = 0; ipf < NSELCALO; ++ipf) { outne[ipf].hwPt = 0; outne[ipf].hwEta = 0; outne[ipf].hwPhi = 0; outne[ipf].hwId = 0; }

    bool isEM[NTRACK];
    // for each track, find the closest calo
    for (int it = 0; it < NTRACK; ++it) {
        isEM[it] = 0; //initialize
        if (track[it].hwPt > 0 && !isMu[it]) {
            int  ibest = best_match_with_pt_ref<NCALO,DR2MAX,HadCaloObj>(hadcalo, track[it]);
            //int  ibest = best_match_ref<NCALO,DR2MAX,true,HadCaloObj>(calo, track[it]);
            if (ibest != -1) {
                if (g_debug_) printf("FW  \t track  %3d pt %7d matched to calo' %3d pt %7d\n", it, int(track[it].hwPt), ibest, int(hadcalo[ibest].hwPt));
                track_good[it] = 1;
                isEM[it] = hadcalo[ibest].hwIsEM;
                calo_sumtk[ibest]    += track[it].hwPt;
                calo_sumtkErr2[ibest] += sqr(track[it].hwPtErr);
            }
        }
    }

    for (int ic = 0; ic < NCALO; ++ic) {
        if (calo_sumtk[ic] > 0) {
            pt_t ptdiff = hadcalo[ic].hwPt - calo_sumtk[ic];
            int sigmamult = (calo_sumtkErr2[ic] + (calo_sumtkErr2[ic] >> 1)); // this multiplies by 1.5 = sqrt(1.5)^2 ~ (1.2)^2
            if (g_debug_ && (hadcalo[ic].hwPt > 0)) {
#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_MORE
                l1tpf_int::CaloCluster floatcalo; fw2dpf::convert(hadcalo[ic], floatcalo); 
                printf("FW  \t calo'  %3d pt %7d [ %7.2f ] eta %+7d [ %+5.2f ] has a sum track pt %7d, difference %7d +- %.2f \n",
                            ic, int(hadcalo[ic].hwPt), floatcalo.floatPt(), int(hadcalo[ic].hwEta), floatcalo.floatEta(), 
                                int(calo_sumtk[ic]), int(ptdiff), std::sqrt(float(int(calo_sumtkErr2[ic]))));
                printf("FW  \t \t\tEMpT %7d pT %7d \n", int(hadcalo[ic].hwEmPt), int(hadcalo[ic].hwPt));
#endif
                        
            }
            if (ptdiff > 0 && ptdiff*ptdiff > sigmamult) {
                if (hadcalo[ic].hwEmPt > 1) {
                    emcalo_pt[ic] = std::min(hadcalo[ic].hwEmPt,ptdiff);
                    ptdiff -= emcalo_pt[ic];
                }
                else emcalo_pt[ic] = 0;
                if (ptdiff > 1) {
                    calo_pt[ic] = ptdiff;
                }
                else calo_pt[ic] = 0;
            } else {
                calo_pt[ic] = 0;
                emcalo_pt[ic] = 0;
            }
        } else {
            emcalo_pt[ic] = hadcalo[ic].hwEmPt;
            calo_pt[ic] = hadcalo[ic].hwPt-hadcalo[ic].hwEmPt;
        }
        if (g_debug_ && (hadcalo[ic].hwPt > 0)) printf("FW  \t calo'  %3d pt %7d ---> %7d \n", ic, int(hadcalo[ic].hwPt), int(calo_pt[ic]));
    }

    // copy out charged hadrons
    for (int it = 0; it < NTRACK; ++it) {
        if (track_good[it]) {
            outch[it].hwPt = track[it].hwPt;
            outch[it].hwEta = track[it].hwEta;
            outch[it].hwPhi = track[it].hwPhi;
            outch[it].hwZ0 = track[it].hwZ0;
            outch[it].hwId  = isMu[it] ? PID_Muon : (isEM[it] ? PID_Electron : PID_Charged);
        }
    }

    // copy out neutral hadrons
    PFNeutralObj outne_all[NCALO];
    for (int ipf = 0; ipf < NCALO; ++ipf) { outne_all[ipf].hwPt = 0; outne_all[ipf].hwEta = 0; outne_all[ipf].hwPhi = 0; outne_all[ipf].hwId = 0; }
    for (int ic = 0; ic < NCALO; ++ic) {
        if (emcalo_pt[ic] > 0) {
            outpho[ic].hwPt  = emcalo_pt[ic];
            outpho[ic].hwEta = hadcalo[ic].hwEta;
            outpho[ic].hwPhi = hadcalo[ic].hwPhi;
            outpho[ic].hwId  = PID_Photon;
        } else {
            outpho[ic].hwPt  = 0;
            outpho[ic].hwEta = 0;
            outpho[ic].hwPhi = 0;
            outpho[ic].hwId  = 0;
        }
        if (calo_pt[ic] > 0) {
            outne_all[ic].hwPt  = calo_pt[ic];
            outne_all[ic].hwEta = hadcalo[ic].hwEta;
            outne_all[ic].hwPhi = hadcalo[ic].hwPhi;
            outne_all[ic].hwId  = PID_Neutral;
        } else {
            outne_all[ic].hwPt  = 0;
            outne_all[ic].hwEta = 0;
            outne_all[ic].hwPhi = 0;
            outne_all[ic].hwId  = 0;
        }
    }

    ptsort_ref<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);

}
