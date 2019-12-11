#include "pfalgo3_ref.h"
#include "pfalgo_common_ref.icc"

#include "DiscretePFInputs.h"
#include "utils/Firmware2DiscretePF.h"
#include <cmath>
#include <algorithm>


int g_pfalgo3_debug_ref_ = 0;

void pfalgo3_ref_set_debug(int debug) { g_pfalgo3_debug_ref_ = debug; }

template<int NCAL, int DR2MAX, bool doPtMin, typename CO_t>
int best_match_ref(const CO_t calo[NCAL], const TkObj & track) {
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
int best_match_ref(const HadCaloObj calo[NCAL], const EmCaloObj & em) {
    pt_t emPtMin = em.hwPt >> 1;
    int  drmin = DR2MAX, ibest = -1;
    for (int ic = 0; ic < NCAL; ++ic) {
            if (calo[ic].hwEmPt <= emPtMin) continue;
            int dr = dr2_int(em.hwEta, em.hwPhi, calo[ic].hwEta, calo[ic].hwPhi);
            if (dr < drmin) { drmin = dr; ibest = ic; }
    }
    return ibest;
}


void pfalgo3_em_ref(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const bool isMu[NTRACK], bool isEle[NTRACK], PFNeutralObj outpho[NPHOTON], HadCaloObj hadcalo_out[NCALO]) {
    // constants
    const int DR2MAX_TE = PFALGO_DR2MAX_TK_EM;
    const int DR2MAX_EH = PFALGO_DR2MAX_EM_CALO;

    // initialize sum track pt
    pt_t calo_sumtk[NEMCALO];
    for (int ic = 0; ic < NEMCALO; ++ic) {  calo_sumtk[ic] = 0; }
    int tk2em[NTRACK]; 
    bool isEM[NEMCALO];
    // for each track, find the closest calo
    for (int it = 0; it < NTRACK; ++it) {
        if (track[it].hwPt > 0 && !isMu[it]) {
            tk2em[it] = best_match_ref<NEMCALO,DR2MAX_TE,false,EmCaloObj>(emcalo, track[it]);
            if (tk2em[it] != -1) {
                if (g_pfalgo3_debug_ref_) printf("FW  \t track  %3d pt %7d matched to em calo %3d pt %7d\n", it, int(track[it].hwPt), tk2em[it], int(emcalo[tk2em[it]].hwPt));
                calo_sumtk[tk2em[it]] += track[it].hwPt;
            }
        } else {
        	tk2em[it] = -1;
        }
    }

    if (g_pfalgo3_debug_ref_) {
        for (int ic = 0; ic < NEMCALO; ++ic) {  if (emcalo[ic].hwPt > 0) printf("FW  \t emcalo %3d pt %7d has sumtk %7d\n", ic, int(emcalo[ic].hwPt), int(calo_sumtk[ic])); }
    }

    for (int ic = 0; ic < NEMCALO; ++ic) {
        pt_t photonPt;
        if (calo_sumtk[ic] > 0) {
            pt_t ptdiff = emcalo[ic].hwPt - calo_sumtk[ic];
            int sigma2 = sqr(emcalo[ic].hwPtErr);
            int sigma2Lo = 4*sigma2, sigma2Hi = sigma2 + (sigma2>>1);
            int ptdiff2 = ptdiff*ptdiff;
            if ((ptdiff > 0 && ptdiff2 <= sigma2Hi) || (ptdiff < 0 && ptdiff2 < sigma2Lo)) {
                // electron
                photonPt = 0; 
                isEM[ic] = true;
                if (g_pfalgo3_debug_ref_) printf("FW  \t emcalo %3d pt %7d ptdiff %7d [match window: -%.2f / +%.2f] flagged as electron\n", ic, int(emcalo[ic].hwPt), int(ptdiff), std::sqrt(float(sigma2Lo)), std::sqrt(float(sigma2Hi)));
            } else if (ptdiff > 0) {
                // electron + photon
                photonPt = ptdiff; 
                isEM[ic] = true;
                if (g_pfalgo3_debug_ref_) printf("FW  \t emcalo %3d pt %7d ptdiff %7d [match window: -%.2f / +%.2f] flagged as electron + photon of pt %7d\n", ic, int(emcalo[ic].hwPt), int(ptdiff), std::sqrt(float(sigma2Lo)), std::sqrt(float(sigma2Hi)), int(photonPt));
            } else {
                // pion
                photonPt = 0;
                isEM[ic] = false;
                if (g_pfalgo3_debug_ref_) printf("FW  \t emcalo %3d pt %7d ptdiff %7d [match window: -%.2f / +%.2f] flagged as pion\n", ic, int(emcalo[ic].hwPt), int(ptdiff), std::sqrt(float(sigma2Lo)), std::sqrt(float(sigma2Hi)));
            }
        } else {
            // photon
            isEM[ic] = true;
            photonPt = emcalo[ic].hwPt;
            if (g_pfalgo3_debug_ref_ && emcalo[ic].hwPt > 0) printf("FW  \t emcalo %3d pt %7d flagged as photon\n", ic, int(emcalo[ic].hwPt));
        }
        outpho[ic].hwPt  = photonPt;
        outpho[ic].hwEta = photonPt ? emcalo[ic].hwEta : etaphi_t(0);
        outpho[ic].hwPhi = photonPt ? emcalo[ic].hwPhi : etaphi_t(0);
        outpho[ic].hwId  = photonPt ? PID_Photon : particleid_t(0);

    }

    for (int it = 0; it < NTRACK; ++it) {
        isEle[it] = (tk2em[it] != -1) && isEM[tk2em[it]];
        if (g_pfalgo3_debug_ref_ && isEle[it]) printf("FW  \t track  %3d pt %7d flagged as electron.\n", it, int(track[it].hwPt));
    }

    int em2calo[NEMCALO];
    for (int ic = 0; ic < NEMCALO; ++ic) {
        em2calo[ic] = best_match_ref<NCALO,DR2MAX_EH>(hadcalo, emcalo[ic]);
        if (g_pfalgo3_debug_ref_ && (emcalo[ic].hwPt > 0)) {
             printf("FW  \t emcalo %3d pt %7d isEM %d matched to hadcalo %7d pt %7d emPt %7d isEM %d\n", 
                                ic, int(emcalo[ic].hwPt), isEM[ic], em2calo[ic], (em2calo[ic] >= 0 ? int(hadcalo[em2calo[ic]].hwPt) : -1), 
                                (em2calo[ic] >= 0 ? int(hadcalo[em2calo[ic]].hwEmPt) : -1), (em2calo[ic] >= 0 ? int(hadcalo[em2calo[ic]].hwIsEM) : 0));
        }
    }
    
    for (int ih = 0; ih < NCALO; ++ih) {
        hadcalo_out[ih] = hadcalo[ih];
        pt_t sub = 0; bool keep = false;
        for (int ic = 0; ic < NEMCALO; ++ic) {
            if (em2calo[ic] == ih) {
                if (isEM[ic]) sub += emcalo[ic].hwPt;
                else keep = true;
            }
        }
        pt_t emdiff  = hadcalo[ih].hwEmPt - sub;
        pt_t alldiff = hadcalo[ih].hwPt - sub;
        if (g_pfalgo3_debug_ref_ && (hadcalo[ih].hwPt > 0)) {
            printf("FW  \t calo   %3d pt %7d has a subtracted pt of %7d, empt %7d -> %7d   isem %d mustkeep %d \n",
                        ih, int(hadcalo[ih].hwPt), int(alldiff), int(hadcalo[ih].hwEmPt), int(emdiff), int(hadcalo[ih].hwIsEM), keep);
                    
        }
        if (alldiff <= ( hadcalo[ih].hwPt >>  4 ) ) {
            hadcalo_out[ih].hwPt = 0;   // kill
            hadcalo_out[ih].hwEmPt = 0; // kill
            if (g_pfalgo3_debug_ref_ && (hadcalo[ih].hwPt > 0)) printf("FW  \t calo   %3d pt %7d --> discarded (zero pt)\n", ih, int(hadcalo[ih].hwPt));
        } else if ((hadcalo[ih].hwIsEM && emdiff <= ( hadcalo[ih].hwEmPt >> 3 )) && !keep) {
            hadcalo_out[ih].hwPt = 0;   // kill
            hadcalo_out[ih].hwEmPt = 0; // kill
            if (g_pfalgo3_debug_ref_ && (hadcalo[ih].hwPt > 0)) printf("FW  \t calo   %3d pt %7d --> discarded (zero em)\n", ih, int(hadcalo[ih].hwPt));
        } else {
            hadcalo_out[ih].hwPt   = alldiff;   
            hadcalo_out[ih].hwEmPt = (emdiff > 0 ? emdiff : pt_t(0)); 
        }
    }
}

void pfalgo3_ref(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) {

    if (g_pfalgo3_debug_ref_) {
#ifdef L1Trigger_Phase2L1ParticleFlow_DiscretePFInputs_MORE
        for (int i = 0; i < NTRACK; ++i) { if (track[i].hwPt == 0) continue;
            l1tpf_impl::PropagatedTrack tk; fw2dpf::convert(track[i], tk); 
            printf("FW  \t track %3d: pt %8d [ %7.2f ]  calo eta %+7d [ %+5.2f ]  calo phi %+7d [ %+5.2f ]  calo ptErr %6d [ %7.2f ] \n", 
                                i, tk.hwPt, tk.floatPt(), tk.hwEta, tk.floatEta(), tk.hwPhi, tk.floatPhi(), tk.hwCaloPtErr, tk.floatCaloPtErr());
        }
        for (int i = 0; i < NEMCALO; ++i) { if (emcalo[i].hwPt == 0) continue;
            l1tpf_impl::CaloCluster em; fw2dpf::convert(emcalo[i], em); 
            printf("FW  \t EM    %3d: pt %8d [ %7.2f ]  calo eta %+7d [ %+5.2f ]  calo phi %+7d [ %+5.2f ]  calo ptErr %6d [ %7.2f ] \n", 
                                i, em.hwPt, em.floatPt(), em.hwEta, em.floatEta(), em.hwPhi, em.floatPhi(), em.hwPtErr, em.floatPtErr());
        } 
        for (int i = 0; i < NCALO; ++i) { if (hadcalo[i].hwPt == 0) continue;
            l1tpf_impl::CaloCluster calo; fw2dpf::convert(hadcalo[i], calo); 
            printf("FW  \t calo  %3d: pt %8d [ %7.2f ]  calo eta %+7d [ %+5.2f ]  calo phi %+7d [ %+5.2f ]  calo emPt %7d [ %7.2f ]   isEM %d \n", 
                                i, calo.hwPt, calo.floatPt(), calo.hwEta, calo.floatEta(), calo.hwPhi, calo.floatPhi(), calo.hwEmPt, calo.floatEmPt(), calo.isEM);
        } 
        for (int i = 0; i < NMU; ++i) { if (mu[i].hwPt == 0) continue;
            l1tpf_impl::Muon muon; fw2dpf::convert(mu[i], muon); 
            printf("FW  \t muon  %3d: pt %8d [ %7.2f ]  muon eta %+7d [ %+5.2f ]  muon phi %+7d [ %+5.2f ]   \n", 
                                i, muon.hwPt, muon.floatPt(), muon.hwEta, muon.floatEta(), muon.hwPhi, muon.floatPhi());
        } 
#endif
    }

    // constants
    const pt_t     TKPT_MAX_LOOSE = PFALGO_TK_MAXINVPT_LOOSE; 
    const pt_t     TKPT_MAX_TIGHT = PFALGO_TK_MAXINVPT_TIGHT; 
    const int      DR2MAX   = PFALGO_DR2MAX_TK_CALO;

    ////////////////////////////////////////////////////
    // TK-MU Linking
    bool isMu[NTRACK];
    pfalgo_mu_ref<NTRACK,NMU>(track, mu, isMu, outmu, g_pfalgo3_debug_ref_);

    ////////////////////////////////////////////////////
    // TK-EM Linking
    bool isEle[NTRACK];
    HadCaloObj hadcalo_subem[NCALO];
    pfalgo3_em_ref(emcalo, hadcalo, track, isMu, isEle, outpho, hadcalo_subem);

    ////////////////////////////////////////////////////
    // TK-HAD Linking

    // initialize sum track pt
    pt_t calo_sumtk[NCALO], calo_subpt[NCALO];
    int  calo_sumtkErr2[NCALO];
    for (int ic = 0; ic < NCALO; ++ic) { calo_sumtk[ic] = 0;  calo_sumtkErr2[ic] = 0;}

    // initialize good track bit
    bool track_good[NTRACK];
    for (int it = 0; it < NTRACK; ++it) { 
        track_good[it] = (track[it].hwPt < (track[it].hwTightQuality ? TKPT_MAX_TIGHT : TKPT_MAX_LOOSE) || isEle[it] || isMu[it]); 
    }

    // initialize output
    for (int ipf = 0; ipf < NTRACK; ++ipf) { outch[ipf].hwPt = 0; outch[ipf].hwEta = 0; outch[ipf].hwPhi = 0; outch[ipf].hwId = 0; outch[ipf].hwZ0 = 0; }
    for (int ipf = 0; ipf < NSELCALO; ++ipf) { outne[ipf].hwPt = 0; outne[ipf].hwEta = 0; outne[ipf].hwPhi = 0; outne[ipf].hwId = 0; }

    // for each track, find the closest calo
    for (int it = 0; it < NTRACK; ++it) {
        if (track[it].hwPt > 0 && !isEle[it] && !isMu[it]) {
            int  ibest = best_match_with_pt_ref<NCALO,DR2MAX,HadCaloObj>(hadcalo_subem, track[it]);
            //int  ibest = best_match_ref<NCALO,DR2MAX,true,HadCaloObj>(hadcalo_subem, track[it]);
            if (ibest != -1) {
                if (g_pfalgo3_debug_ref_) printf("FW  \t track  %3d pt %7d matched to calo %3d pt %7d\n", it, int(track[it].hwPt), ibest, int(hadcalo_subem[ibest].hwPt));
                track_good[it] = 1;
                calo_sumtk[ibest]    += track[it].hwPt;
                calo_sumtkErr2[ibest] += sqr(track[it].hwPtErr);
            }
        }
    }

    for (int ic = 0; ic < NCALO; ++ic) {
        if (calo_sumtk[ic] > 0) {
            pt_t ptdiff = hadcalo_subem[ic].hwPt - calo_sumtk[ic];
            int sigmamult = (calo_sumtkErr2[ic] + (calo_sumtkErr2[ic] >> 1)); // this multiplies by 1.5 = sqrt(1.5)^2 ~ (1.2)^2
            if (g_pfalgo3_debug_ref_ && (hadcalo_subem[ic].hwPt > 0)) {
#ifdef L1Trigger_Phase2L1ParticleFlow_DiscretePFInputs_MORE
                l1tpf_impl::CaloCluster floatcalo; fw2dpf::convert(hadcalo_subem[ic], floatcalo); 
                printf("FW  \t calo  %3d pt %7d [ %7.2f ] eta %+7d [ %+5.2f ] has a sum track pt %7d, difference %7d +- %.2f \n",
                            ic, int(hadcalo_subem[ic].hwPt), floatcalo.floatPt(), int(hadcalo_subem[ic].hwEta), floatcalo.floatEta(), 
                                int(calo_sumtk[ic]), int(ptdiff), std::sqrt(float(int(calo_sumtkErr2[ic]))));
#endif
                        
            }
            if (ptdiff > 0 && ptdiff*ptdiff > sigmamult) {
                calo_subpt[ic] = ptdiff;
            } else {
                calo_subpt[ic] = 0;
            }
        } else {
            calo_subpt[ic] = hadcalo_subem[ic].hwPt;
        }
        if (g_pfalgo3_debug_ref_ && (hadcalo_subem[ic].hwPt > 0)) printf("FW  \t calo  %3d pt %7d ---> %7d \n", ic, int(hadcalo_subem[ic].hwPt), int(calo_subpt[ic]));
    }

    // copy out charged hadrons
    for (int it = 0; it < NTRACK; ++it) {
        if (track_good[it]) {
            outch[it].hwPt = track[it].hwPt;
            outch[it].hwEta = track[it].hwEta;
            outch[it].hwPhi = track[it].hwPhi;
            outch[it].hwZ0 = track[it].hwZ0;
            outch[it].hwId  = isEle[it] ? PID_Electron : (isMu[it] ? PID_Muon : PID_Charged);
        }
    }

    // copy out neutral hadrons
    PFNeutralObj outne_all[NCALO];
    for (int ipf = 0; ipf < NCALO; ++ipf) { outne_all[ipf].hwPt = 0; outne_all[ipf].hwEta = 0; outne_all[ipf].hwPhi = 0; outne_all[ipf].hwId = 0; }
    for (int ic = 0; ic < NCALO; ++ic) {
        if (calo_subpt[ic] > 0) {
            outne_all[ic].hwPt  = calo_subpt[ic];
            outne_all[ic].hwEta = hadcalo_subem[ic].hwEta;
            outne_all[ic].hwPhi = hadcalo_subem[ic].hwPhi;
            outne_all[ic].hwId  = PID_Neutral;
        }
    }

    ptsort_ref<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);
}
