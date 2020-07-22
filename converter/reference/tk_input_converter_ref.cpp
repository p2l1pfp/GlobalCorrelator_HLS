// only needed for the constants
#include "../firmware/tk_input_converter.h"

//inline float tanlam_to_eta_ref(float tanlam){return -log(tan((M_PI/2 - atan(tanlam))/2));}

float propagate_tanlam_ref(float z0, float tanlam){
    float tanlam_calo;
    float z_at_calo = z0 + DETR * tanlam;
    if(z_at_calo < DETZ){
        tanlam_calo = z0/DETR + tanlam;
    } else {
        tanlam_calo = tanlam * DETZ/(DETZ-z0);
    }
    return tanlam_calo;
}

float convert_dphi_ref(float pt){
    // r of curv is pt * (looper radius / 2)/(looper pt)
    float rCurv = pt * (129./2)/(0.735);
    float x = (DETR)/(2*rCurv);
    return atan(x / sqrt(1-x*x));
}

float reso_calo_ref(float pt, float eta_calo){
    if( fabs(eta_calo) < 0.700 ){
        return pt * 0.122 + 2.582;
    } else if( fabs(eta_calo) < 1.200 ){
        return pt * 0.143 + 2.191;
    } else {
        return pt * 0.465 + -0.077;
    }
}

void pf_input_track_conv_ref(float rinv           ,
                             float tkphi          ,
                             float tanlam         ,
                             float tkz0           ,
                             float tkd0           ,
                             float chi2rphi       ,
                             float chi2rz         ,
                             float bendChi2       ,
                             int hit              ,
                             int trackMVA         ,
                             int extraMVA         ,
                             int valid            ,
                             float& pf_pt         , 
                             float& pf_pterr      , 
                             float& pf_eta_at_calo, 
                             float& pf_phi_at_calo, 
                             float& pf_z0         , 
                             bool& pf_TightQuality){
    pf_pt = 1./rinv;
    pf_phi_at_calo = tkphi + (rinv>0?1:-1)*convert_dphi_ref(pf_pt);
    float tl_at_calo = propagate_tanlam_ref(tkz0, tanlam);
    pf_eta_at_calo = tanlam_to_eta_ref( tl_at_calo );
    pf_pterr = reso_calo_ref(pf_pt, pf_eta_at_calo);
    pf_z0 = tkz0;
    pf_TightQuality = 1;
}

