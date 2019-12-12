#include "pfalgo2hgc_ref.h"
#include "pfalgo_common_ref.icc"

#include "DiscretePFInputs.h"
#include "utils/Firmware2DiscretePF.h"
#include <cmath>
#include <algorithm>


int g_pfalgo2hgc_debug_ref_ = 0;

void pfalgo2hgc_ref_set_debug(int debug) { g_pfalgo2hgc_debug_ref_ = debug; }

void pfalgo2hgc_ref(const HadCaloObj calo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) {

    if (g_pfalgo2hgc_debug_ref_) {
#ifdef L1Trigger_Phase2L1ParticleFlow_DiscretePFInputs_MORE
        for (int i = 0; i < NTRACK; ++i) { if (track[i].hwPt == 0) continue;
            l1tpf_impl::PropagatedTrack tk; fw2dpf::convert(track[i], tk); 
            printf("FW  \t track %3d: pt %8d [ %7.2f ]  calo eta %+7d [ %+5.2f ]  calo phi %+7d [ %+5.2f ]  calo ptErr %6d [ %7.2f ] \n", 
                                i, tk.hwPt, tk.floatPt(), tk.hwEta, tk.floatEta(), tk.hwPhi, tk.floatPhi(), tk.hwCaloPtErr, tk.floatCaloPtErr());
        }
        for (int i = 0; i < NCALO; ++i) { if (calo[i].hwPt == 0) continue;
            l1tpf_impl::CaloCluster c; fw2dpf::convert(calo[i], c); 
            printf("FW  \t calo  %3d: pt %8d [ %7.2f ]  calo eta %+7d [ %+5.2f ]  calo phi %+7d [ %+5.2f ]  calo emPt %7d [ %7.2f ]   isEM %d \n", 
                                i, c.hwPt, c.floatPt(), c.hwEta, c.floatEta(), c.hwPhi, c.floatPhi(), c.hwEmPt, c.floatEmPt(), c.isEM);
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
    const int      DR2MAX         = PFALGO_DR2MAX_TK_CALO;

    ////////////////////////////////////////////////////
    // TK-MU Linking

    bool isMu[NTRACK];
    pfalgo_mu_ref<NTRACK,NMU>(track, mu, isMu, outmu, g_pfalgo2hgc_debug_ref_ );


    ////////////////////////////////////////////////////
    // TK-HAD Linking

    // initialize sum track pt
    pt_t calo_sumtk[NCALO], calo_subpt[NCALO];
    int  calo_sumtkErr2[NCALO];
    for (int ic = 0; ic < NCALO; ++ic) { calo_sumtk[ic] = 0;  calo_sumtkErr2[ic] = 0;}

    // initialize good track bit
    bool track_good[NTRACK];
    bool isEle[NTRACK];
    for (int it = 0; it < NTRACK; ++it) { 
        track_good[it] = (track[it].hwPt < (track[it].hwTightQuality ? TKPT_MAX_TIGHT : TKPT_MAX_LOOSE) || isMu[it]); 
        isEle[it] = false;
    }

    // initialize output
    for (int ipf = 0; ipf < NTRACK; ++ipf) clear(outch[ipf]); 
    for (int ipf = 0; ipf < NSELCALO; ++ipf) clear(outne[ipf]);

    // for each track, find the closest calo
    for (int it = 0; it < NTRACK; ++it) {
        if (track[it].hwPt > 0 && !isMu[it]) {
            int  ibest = best_match_with_pt_ref<NCALO,DR2MAX,HadCaloObj>(calo, track[it]);
            if (ibest != -1) {
                if (g_pfalgo2hgc_debug_ref_) printf("FW  \t track  %3d pt %7d matched to calo' %3d pt %7d\n", it, int(track[it].hwPt), ibest, int(calo[ibest].hwPt));
                track_good[it] = 1;
                isEle[it] = calo[ibest].hwIsEM;
                calo_sumtk[ibest]    += track[it].hwPt;
                calo_sumtkErr2[ibest] += sqr(track[it].hwPtErr);
            }
        }
    }

    for (int ic = 0; ic < NCALO; ++ic) {
        if (calo_sumtk[ic] > 0) {
            pt_t ptdiff = calo[ic].hwPt - calo_sumtk[ic];
            int sigmamult = (calo_sumtkErr2[ic] + (calo_sumtkErr2[ic] >> 1)); // this multiplies by 1.5 = sqrt(1.5)^2 ~ (1.2)^2
            if (g_pfalgo2hgc_debug_ref_ && (calo[ic].hwPt > 0)) {
#ifdef L1Trigger_Phase2L1ParticleFlow_DiscretePFInputs_MORE
                l1tpf_impl::CaloCluster floatcalo; fw2dpf::convert(calo[ic], floatcalo); 
                printf("FW  \t calo'  %3d pt %7d [ %7.2f ] eta %+7d [ %+5.2f ] has a sum track pt %7d, difference %7d +- %.2f \n",
                            ic, int(calo[ic].hwPt), floatcalo.floatPt(), int(calo[ic].hwEta), floatcalo.floatEta(), 
                                int(calo_sumtk[ic]), int(ptdiff), std::sqrt(float(int(calo_sumtkErr2[ic]))));
#endif
                        
            }
            if (ptdiff > 0 && ptdiff*ptdiff > sigmamult) {
                calo_subpt[ic] = ptdiff;
            } else {
                calo_subpt[ic] = 0;
            }
        } else {
            calo_subpt[ic] = calo[ic].hwPt;
        }
        if (g_pfalgo2hgc_debug_ref_ && (calo[ic].hwPt > 0)) printf("FW  \t calo'  %3d pt %7d ---> %7d \n", ic, int(calo[ic].hwPt), int(calo_subpt[ic]));
    }

    // copy out charged hadrons
    for (int it = 0; it < NTRACK; ++it) {
        if (track_good[it]) {
            assert(!(isEle[it] && isMu[it]));
            outch[it].hwPt = track[it].hwPt;
            outch[it].hwEta = track[it].hwEta;
            outch[it].hwPhi = track[it].hwPhi;
            outch[it].hwZ0 = track[it].hwZ0;
            outch[it].hwId  = isEle[it] ? PID_Electron : (isMu[it] ? PID_Muon : PID_Charged);
        }
    }

#if NCALO == NSELCALO
    // copy out neutral hadrons directly without sorting
    for (int ic = 0; ic < NCALO; ++ic) {
        clear(outne[ic])
        if (calo_subpt[ic] > 0) {
            outne[ic].hwPt  = calo_subpt[ic];
            outne[ic].hwEta = calo[ic].hwEta;
            outne[ic].hwPhi = calo[ic].hwPhi;
            outne[ic].hwId  = calo[ic].hwIsEM ? PID_Photon : PID_Neutral;
        }
    }
#else
    // copy out neutral hadrons with sorting and cropping
    PFNeutralObj outne_all[NCALO];
    for (int ipf = 0; ipf < NCALO; ++ipf) clear(outne_all[ipf]);
    for (int ic = 0; ic < NCALO; ++ic) {
        if (calo_subpt[ic] > 0) {
            outne_all[ic].hwPt  = calo_subpt[ic];
            outne_all[ic].hwEta = calo[ic].hwEta;
            outne_all[ic].hwPhi = calo[ic].hwPhi;
            outne_all[ic].hwId  = calo[ic].hwIsEM ? PID_Photon : PID_Neutral;
        }
    }

    ptsort_ref<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);
#endif

}
