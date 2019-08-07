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


void pfalgo3_ne_ref(HadCaloObj hadcalo[NCALO], PFNeutralObj outpho[NCALO], PFNeutralObj outne[NCALO]) {

    for (int ih = 0; ih < NCALO; ++ih) {
        if (g_debug_ && (hadcalo[ih].hwPt > 0)) {
            printf("FW  \t calo   %3d has a pt of %7d, empt %7d   isem %d \n",
                        ih, int(hadcalo[ih].hwPt), int(hadcalo[ih].hwEmPt), int(hadcalo[ih].hwIsEM));
                    
        }
        if ((hadcalo[ih].hwPt > 0) && hadcalo[ih].hwIsEM) {
            outpho[ih].hwPt  = hadcalo[ih].hwPt;
            outpho[ih].hwEta = hadcalo[ih].hwEta;
            outpho[ih].hwPhi = hadcalo[ih].hwPhi;
            outpho[ih].hwId  = PID_Photon;
            outne[ih].hwPt  = pt_t(0);
            outne[ih].hwEta = etaphi_t(0);
            outne[ih].hwPhi = etaphi_t(0);
            outne[ih].hwId  = PID_Neutral;
            if (g_debug_ && (hadcalo[ih].hwPt > 0)) printf("FW  \t calo   %3d pt %7d --> PHOTON\n", ih, int(hadcalo[ih].hwPt));
        } else if (hadcalo[ih].hwPt > 0) {
            outpho[ih].hwPt  = pt_t(0);
            outpho[ih].hwEta = etaphi_t(0);
            outpho[ih].hwPhi = etaphi_t(0);
            outpho[ih].hwId  = PID_Photon;
            outne[ih].hwPt  = hadcalo[ih].hwPt;
            outne[ih].hwEta = hadcalo[ih].hwEta;
            outne[ih].hwPhi = hadcalo[ih].hwPhi;
            outne[ih].hwId  = PID_Neutral;
            if (g_debug_ && (hadcalo[ih].hwPt > 0)) printf("FW  \t calo   %3d pt %7d --> NEUTRAL HADRON\n", ih, int(hadcalo[ih].hwPt));
        } else {
            outpho[ih].hwPt  = pt_t(0);
            outpho[ih].hwEta = etaphi_t(0);
            outpho[ih].hwPhi = etaphi_t(0);
            outpho[ih].hwId  = PID_Photon;
            outne[ih].hwPt  = pt_t(0);
            outne[ih].hwEta = etaphi_t(0);
            outne[ih].hwPhi = etaphi_t(0);
            outne[ih].hwId  = PID_Neutral;
        }
    }
}

void pfalgo3_forward_ref(HadCaloObj hadcalo[NCALO], MuObj mu[NMU], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) {

    if (g_debug_) {
#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_MORE
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

    PFNeutralObj outpho_all[NCALO];
    for (int ipf = 0; ipf < NCALO; ++ipf) { outpho_all[ipf].hwPt = 0; outpho_all[ipf].hwEta = 0; outpho_all[ipf].hwPhi = 0; outpho_all[ipf].hwId = 0; }
    PFNeutralObj outne_all[NCALO];
    for (int ipf = 0; ipf < NCALO; ++ipf) { outne_all[ipf].hwPt = 0; outne_all[ipf].hwEta = 0; outne_all[ipf].hwPhi = 0; outne_all[ipf].hwId = 0; }

    ////////////////////////////////////////////////////
    // Converting
    pfalgo3_ne_ref(hadcalo, outpho_all, outne_all);

    // initialize output
    for (int ipf = 0; ipf < NPHOTON; ++ipf) { outpho[ipf].hwPt = 0; outpho[ipf].hwEta = 0; outpho[ipf].hwPhi = 0; outpho[ipf].hwId = 0; }
    for (int ipf = 0; ipf < NSELCALO; ++ipf) { outne[ipf].hwPt = 0; outne[ipf].hwEta = 0; outne[ipf].hwPhi = 0; outne[ipf].hwId = 0; }

    ptsort_ref<PFNeutralObj,NCALO,NPHOTON>(outpho_all, outpho);
    ptsort_ref<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);

}
