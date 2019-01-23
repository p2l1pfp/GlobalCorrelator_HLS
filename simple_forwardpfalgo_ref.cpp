#include "firmware/data.h"
#include "firmware/simple_forwardpfalgo.h"
#include "DiscretePFInputs.h"
#include "utils/Firmware2DiscretePF.h"
#include <cmath>
#include <algorithm>

bool g_debug_ = 0;

void pfalgo3_forward_ref_set_debug(bool debug) { g_debug_ = debug; }

template <typename T> int sqr(const T & t) { return t*t; }

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


void pfalgo3_em_ref(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], PFNeutralObj outpho[NPHOTON], HadCaloObj hadcalo_out[NCALO]) {
    // constants
    const int DR2MAX_EH = PFALGO3_DR2MAX_EM_CALO;

    for (int ic = 0; ic < NEMCALO; ++ic) {
        pt_t photonPt = emcalo[ic].hwPt;
        if (g_debug_ && emcalo[ic].hwPt > 0) printf("FW  \t emcalo %3d pt %7d flagged as photon\n", ic, int(emcalo[ic].hwPt));
        outpho[ic].hwPt  = photonPt;
        outpho[ic].hwEta = photonPt ? emcalo[ic].hwEta : etaphi_t(0);
        outpho[ic].hwPhi = photonPt ? emcalo[ic].hwPhi : etaphi_t(0);
        outpho[ic].hwId  = photonPt ? PID_Photon : particleid_t(0);

    }

    int em2calo[NEMCALO];
    for (int ic = 0; ic < NEMCALO; ++ic) {
        em2calo[ic] = best_match_ref<NCALO,DR2MAX_EH>(hadcalo, emcalo[ic]);
        if (g_debug_ && (emcalo[ic].hwPt > 0)) {
             printf("FW  \t emcalo %3d pt %7d matched to hadcalo %7d pt %7d emPt %7d isEM %d\n", 
                                ic, int(emcalo[ic].hwPt), em2calo[ic], (em2calo[ic] >= 0 ? int(hadcalo[em2calo[ic]].hwPt) : -1), 
                                (em2calo[ic] >= 0 ? int(hadcalo[em2calo[ic]].hwEmPt) : -1), (em2calo[ic] >= 0 ? int(hadcalo[em2calo[ic]].hwIsEM) : 0));
        }
    }
    
    for (int ih = 0; ih < NCALO; ++ih) {
        hadcalo_out[ih] = hadcalo[ih];
        pt_t sub = 0; bool keep = false;
        for (int ic = 0; ic < NEMCALO; ++ic) {
            if (em2calo[ic] == ih) {
                sub += emcalo[ic].hwPt;
            }
        }
        pt_t emdiff  = hadcalo[ih].hwEmPt - sub;
        pt_t alldiff = hadcalo[ih].hwPt - sub;
        if (g_debug_ && (hadcalo[ih].hwPt > 0)) {
            printf("FW  \t calo   %3d pt %7d has a subtracted pt of %7d, empt %7d -> %7d   isem %d keep %d \n",
                        ih, int(hadcalo[ih].hwPt), int(alldiff), int(hadcalo[ih].hwEmPt), int(emdiff), int(hadcalo[ih].hwIsEM), keep);
                    
        }
        if (alldiff <= ( hadcalo[ih].hwPt >>  4 ) ) {
            hadcalo_out[ih].hwPt = 0;   // kill
            hadcalo_out[ih].hwEmPt = 0; // kill
            if (g_debug_ && (hadcalo[ih].hwPt > 0)) printf("FW  \t calo   %3d pt %7d --> discarded (zero pt)\n", ih, int(hadcalo[ih].hwPt));
        } else if ((hadcalo[ih].hwIsEM && emdiff <= ( hadcalo[ih].hwEmPt >> 3 )) && !keep) {
            hadcalo_out[ih].hwPt = 0;   // kill
            hadcalo_out[ih].hwEmPt = 0; // kill
            if (g_debug_ && (hadcalo[ih].hwPt > 0)) printf("FW  \t calo   %3d pt %7d --> discarded (zero em)\n", ih, int(hadcalo[ih].hwPt));
        } else {
            hadcalo_out[ih].hwPt   = alldiff;   
            hadcalo_out[ih].hwEmPt = (emdiff > 0 ? emdiff : pt_t(0)); 
        }
    }
}

void pfalgo3_forward_ref(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], MuObj mu[NMU], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) {

    if (g_debug_) {
#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_MORE
        for (int i = 0; i < NEMCALO; ++i) { if (emcalo[i].hwPt == 0) continue;
            l1tpf_int::CaloCluster em; fw2dpf::convert(emcalo[i], em); 
            printf("FW  \t EM    %3d: pt %8d [ %7.2f ]  calo eta %+7d [ %+5.2f ]  calo phi %+7d [ %+5.2f ]  calo ptErr %6d [ %7.2f ] \n", 
                                i, em.hwPt, em.floatPt(), em.hwEta, em.floatEta(), em.hwPhi, em.floatPhi(), em.hwPtErr, em.floatPtErr());
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

    ////////////////////////////////////////////////////
    // TK-MU Linking

    // initialize output
    for (int ipf = 0; ipf < NMU; ++ipf) { outmu[ipf].hwPt = 0; outmu[ipf].hwEta = 0; outmu[ipf].hwPhi = 0; outmu[ipf].hwId  = 0; outmu[ipf].hwZ0  = 0; }

    for (int im = 0; im < NMU; ++im) {
        if (mu[im].hwPt > 0) {
                outmu[im].hwPt = mu[im].hwPt;
                outmu[im].hwEta = mu[im].hwEta;
                outmu[im].hwPhi = mu[im].hwPhi;
                outmu[im].hwId  = PID_Muon;
                if (g_debug_) printf("FW  \t muon %3d passed \n", im);
        }
    }

    ////////////////////////////////////////////////////
    // TK-EM Linking
    HadCaloObj hadcalo_subem[NCALO];
    pfalgo3_em_ref(emcalo, hadcalo, outpho, hadcalo_subem);

    ////////////////////////////////////////////////////
    // TK-HAD Linking

    // initialize output
    for (int ipf = 0; ipf < NSELCALO; ++ipf) { outne[ipf].hwPt = 0; outne[ipf].hwEta = 0; outne[ipf].hwPhi = 0; outne[ipf].hwId = 0; }

    // copy out neutral hadrons
    PFNeutralObj outne_all[NCALO];
    for (int ipf = 0; ipf < NCALO; ++ipf) { outne_all[ipf].hwPt = 0; outne_all[ipf].hwEta = 0; outne_all[ipf].hwPhi = 0; outne_all[ipf].hwId = 0; }
    for (int ic = 0; ic < NCALO; ++ic) {
        if (hadcalo_subem[ic].hwPt > 0) {
            outne_all[ic].hwPt  = hadcalo_subem[ic].hwPt;
            outne_all[ic].hwEta = hadcalo_subem[ic].hwEta;
            outne_all[ic].hwPhi = hadcalo_subem[ic].hwPhi;
            outne_all[ic].hwId  = PID_Neutral;
        }
    }

    ptsort_ref<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);

}
