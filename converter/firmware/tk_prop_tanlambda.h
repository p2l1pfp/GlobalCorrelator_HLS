//

#include "tk_input_converter.h"

//void propagate_tanlam(tkz0_t z0, tanlam_t tanlam, tanlam_t &tanlam_at_det);

// https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-integration-CMSSW_11_1_0_pre6/L1Trigger/Phase2L1ParticleFlow/src/L1TPFUtils.cc

void propagate_tanlam(tkz0_t z0, tanlam_t tanlam, tanlam_t &tanlam_at_det){
    // simplify for barrel
    tanlam_at_det = tanlam + tanlam_help_t(z0) * tanlam_help_t(1./DETR);
    return;

    // using a helper type to handle necessary range
    zdet_t z_at_det = z0 + zdet_t(tanlam) * DETR;
    // std::cout << z0 << "  " << z_at_det << " -- ";
    if(z_at_det < -DETZ){
        // neg endcap
        tanlam_at_det=0; 
        // tanlam_at_det = tanlam / (1 + tanlam_help_t(z0) * tanlam_help_t(1/DETZ));
        // use power expansion for now: 1/(1+z0/z_det) -> 1-z0/z_det
        // tanlam_at_det = tanlam * (1 - tanlam_help_t(z0) * tanlam_help_t(1/DETZ));
    } else if(z_at_det > DETZ){
        // pos endcap
        tanlam_at_det=0; 
        // tanlam_at_det = tanlam / (1 - tanlam_help_t(z0) * tanlam_help_t(1/DETZ));
        // use power expansion for now: 1/(1-z0/z_det) -> 1+z0/z_det
        // tanlam_at_det = tanlam * (1 + tanlam_help_t(z0) * tanlam_help_t(1/DETZ));
    } else {
        tanlam_at_det = tanlam + tanlam_help_t(z0) * tanlam_help_t(1./DETR);
    }
}




    // https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-integration-CMSSW_11_1_0_pre6/L1Trigger/Phase2L1ParticleFlow/plugins/PFTrackProducerFromL1Tracks.cc
    // reco::Candidate::PolarLorentzVector p4p(pt, eta, phi, 0.137); // pion mass
    // reco::Particle::LorentzVector p4(p4p.X(), p4p.Y(), p4p.Z(), p4p.E());
    // reco::Particle::Point vtx(0.,0.,z0);
    // auto caloetaphi = l1tpf::propagateToCalo(p4, math::XYZTLorentzVector(0.,0.,z0,0.), charge, fBz_);

    // https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-integration-CMSSW_11_1_0_pre6/L1Trigger/Phase2L1ParticleFlow/src/L1TPFUtils.cc
    // std::pair<float,float> l1tpf::propagateToCalo(const math::XYZTLorentzVector& iMom, const math::XYZTLorentzVector& iVtx, double iCharge, double iBField) {
    //     BaseParticlePropagator prop = BaseParticlePropagator(RawParticle(iMom,iVtx,iCharge),0.,0.,iBField);
    //     prop.propagateToEcalEntrance(false);
    //     double ecalShowerDepth = reco::PFCluster::getDepthCorrection(prop.particle().momentum().E(),false,false);
    //     math::XYZVector point = math::XYZVector(prop.particle().vertex())+math::XYZTLorentzVector(prop.particle().momentum()).Vect().Unit()*ecalShowerDepth;
    //     return std::make_pair(point.eta(), point.phi());
    // }
    // https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-integration-CMSSW_11_1_0_pre6/CommonTools/BaseParticlePropagator/src/BaseParticlePropagator.cc
    // https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-integration-CMSSW_11_1_0_pre6/CommonTools/BaseParticlePropagator/interface/RawParticle.h
    
    // https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-integration-CMSSW_11_1_0_pre6/CommonTools/BaseParticlePropagator/src/BaseParticlePropagator.cc#L457
    // R Ecal entrance = 129 cm
    // Z Ecal entrance = 320.9 cm

    // Tan lambda propogation
    //
    // Barrel:
    // tam(lambda_det) = z_det / R_det
    // z@det = z0 + R_det * tam(lambda)
    // tan(lambda_det) = z0/R_det + tam(lambda)
    //
    // Endcap: if z@det above is greater than z_det
    // r@zmax = (z_det - z0) / tan(lambda)
    // tan(lambda_det) = tan(lambda) * (z_det / (z_det - z0))
    
    // dPhi for track with curvature radius C and detector radius D
    // detector is x^2 + y^2 = D^2
    // particle trajectory (x-C)^2 + y^2 = C^2
    // x=D^2/2C, y=D*sqrt(1-D^2/4C^2)
    // tan phi = 2C/D*sqrt(1-D^2/4C^2)
