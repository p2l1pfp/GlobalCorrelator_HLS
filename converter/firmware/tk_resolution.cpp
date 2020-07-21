#include "tk_resolution.h"

void reso_calo(pt_t pt, etaphi_t eta_calo, pt_t& err){
    // Resolution from:
    // https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-integration-CMSSW_11_1_0_pre6/L1Trigger/Phase2L1ParticleFlow/python/pfTracksFromL1Tracks_cfi.py
    // https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-integration-CMSSW_11_1_0_pre6/L1Trigger/Phase2L1ParticleFlow/python/l1ParticleFlow_cff.py
    // barrel reso at cal
    // etaBins = cms.vdouble( 0.700,  1.200,  1.600),
    // offset  = cms.vdouble( 2.582,  2.191, -0.077),
    // scale   = cms.vdouble( 0.122,  0.143,  0.465),
    // track reso
    // etaBins = cms.vdouble( 0.800,  1.200,  1.500,  2.000,  2.500),
    // offset  = cms.vdouble( 0.007,  0.009,  0.011,  0.015,  0.025),
    // scale   = cms.vdouble( 0.275,  0.404,  0.512,  0.480,  1.132),

    // TODO currently barrel only, and requires 3 MACs
    etaphi_t abs_eta = eta_calo;
    if(eta_calo<0) abs_eta = -eta_calo;
    if(abs_eta < etaphi_t(0.700 * PF_ETAPHI_SCALE)){
        err = bigfix_t(pt) * bigfix_t(0.122) + bigfix_t(2.582 * PF_PT_SCALE);
    } else if(abs_eta < etaphi_t(1.200 * PF_ETAPHI_SCALE)){
        err = bigfix_t(pt) * bigfix_t(0.143) + bigfix_t(2.191 * PF_PT_SCALE);
    } else {
        err = bigfix_t(pt) * bigfix_t(0.465) + bigfix_t(-0.077 * PF_PT_SCALE);
    }
    if (err > pt) err = pt;
    
}
