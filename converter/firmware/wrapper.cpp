/*
HLS implementation of input conversion wrapper
*/
#include "wrapper.h"
#ifndef __SYNTHESIS__
#include <cstdio>
using std::cout;
using std::endl;
#endif

void pf_input_track_conv_hw(input_t in, output_t& out, numlink_t nlink){
    #pragma HLS pipeline II=1

    // unpack L1Track format
    rinv_t     rinv    ;
    tkphi_t    tkphi   ;
    tanlam_t   tanlam  ;
    tkz0_t     tkz0    ;
    tkd0_t     tkd0    ;
    chi2rphi_t chi2rphi;
    chi2rz_t   chi2rz  ;
    bendChi2_t bendChi2;
    hit_t      hit     ;
    trackMVA_t trackMVA;
    extraMVA_t extraMVA;
    valid_t    valid   ;

    unpack_L1T_track(in, rinv, tkphi, tanlam, tkz0, tkd0, chi2rphi, chi2rz, bendChi2, hit, trackMVA, extraMVA, valid);

    // selection
    if(!valid){
        out=0;
        return;
    }

    // track propagation to calo surface in eta
    tanlam_t tanlam_at_calo;
    propagate_tanlam(tkz0, tanlam, tanlam_at_calo);

    // track propagation to calo surface in phi
    rinv_t pt_inv = rinv;
    etaphi_t dphi;
    if (rinv<0) pt_inv = -rinv;
    convert_dphi(pt_inv, dphi);

    etaphi_t pf_phi_at_calo;
    if (rinv<0) pf_phi_at_calo = tkphi - dphi;
    else pf_phi_at_calo = tkphi + dphi;


    // converters
    pt_t conv_pt;
    //std::cout << "inverse pt " << rinv << std::endl;
    convert_pt(pt_inv,conv_pt);
    //std::cout << "        pt " << conv_pt << std::endl;

    etaphi_t pf_eta_at_calo;
    convert_eta(tanlam_at_calo, pf_eta_at_calo);



    //     // initialize conversion constants
    // #ifdef __HLS_SYN__
    //     bool initialized = false;
    //     pt_t phi_offsets[TT_NPHI_SECTORS];
    // #else 
    //     static bool initialized = false;
    //     static pt_t phi_offsets[TT_NPHI_SECTORS];
    // #endif
    //     if (!initialized) {
    //         for(int i=0;i<TT_NPHI_SECTORS;i++)
    //             phi_offsets[i] = (PF_ETAPHI_SCALE * M_PI / TT_NPHI_SECTORS)*(2*i-1);
    //         initialized = true;
    //     }
    //  etaphi_t pf_phi = bigfix_t(tkphi)*bigfix_t(PF_ETAPHI_SCALE) + phi_offsets[nlink];

    //etaphi_t pf_phi = tkphi;

    // z0_t pf_z0 = bigfix_t(tkz0)*bigfix_t(PF_Z0_SCALE); // scale=20
    // for now, copy z0 values without explicitly converting... must be replaced
    z0_t pf_z0 = zdet_t(tkz0) * zdet_t(PF_Z0_SCALE);

    // eta phi propogation
    // https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-integration-CMSSW_11_1_0_pre6/L1Trigger/Phase2L1ParticleFlow/src/L1TPFUtils.cc
    // etaphi_t pf_eta_calo;
    // etaphi_t pf_phi_calo;
    // track_propagate(conv_pt, conv_eta, pf_phi, pf_z0, pf_eta_calo, pf_phi_calo);

    // pack in PF format
    pt_t pf_pt = conv_pt;
    pt_t pf_pterr;
    reso_calo(pf_pt, pf_eta_at_calo, pf_pterr); // needs eta calo
    bool pf_TightQuality = 1;
    pack_pf_track(out, pf_pt, pf_pterr, pf_eta_at_calo, pf_phi_at_calo, pf_z0, pf_TightQuality);
    //std::cout << out.to_string(16) << std::endl;
}

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

void pack_L1T_track(ap_uint<kTrackWordSize> &tk,
                    rinv_t     rinv    ,
                    tkphi_t    tkphi   ,
                    tanlam_t   tanlam  ,
                    tkz0_t     tkz0    ,
                    tkd0_t     tkd0    ,
                    chi2rphi_t chi2rphi,
                    chi2rz_t   chi2rz  ,
                    bendChi2_t bendChi2,
                    hit_t      hit     ,
                    trackMVA_t trackMVA,
                    extraMVA_t extraMVA,
                    valid_t    valid   ){

    ap_uint<rinv_t    ::width> word_rinv      = rinv    .range(rinv_t    ::width-1,0);
    ap_uint<tkphi_t   ::width> word_tkphi     = tkphi   .range(tkphi_t   ::width-1,0);
    ap_uint<tanlam_t  ::width> word_tanlam    = tanlam  .range(tanlam_t  ::width-1,0);
    ap_uint<tkz0_t    ::width> word_tkz0      = tkz0    .range(tkz0_t    ::width-1,0);
    ap_uint<tkd0_t    ::width> word_tkd0      = tkd0    .range(tkd0_t    ::width-1,0);
    ap_uint<chi2rphi_t::width> word_chi2rphi  = chi2rphi.range(chi2rphi_t::width-1,0);
    ap_uint<chi2rz_t  ::width> word_chi2rz    = chi2rz  .range(chi2rz_t  ::width-1,0);
    ap_uint<bendChi2_t::width> word_bendChi2  = bendChi2.range(bendChi2_t::width-1,0);
    ap_uint<hit_t     ::width> word_hit       = hit     .range(hit_t     ::width-1,0);
    ap_uint<trackMVA_t::width> word_trackMVA  = trackMVA.range(trackMVA_t::width-1,0);
    ap_uint<extraMVA_t::width> word_extraMVA  = extraMVA.range(extraMVA_t::width-1,0);
    ap_uint<valid_t   ::width> word_valid     = valid   .range(valid_t   ::width-1,0);

    tk = (word_valid, word_extraMVA, word_trackMVA, word_hit, word_bendChi2, word_chi2rz, 
          word_chi2rphi, word_tkd0, word_tkz0, word_tanlam, word_tkphi, word_rinv);    
}

void unpack_L1T_track(ap_uint<kTrackWordSize> in,
                    rinv_t     &rinv    ,
                    tkphi_t    &tkphi   ,
                    tanlam_t   &tanlam  ,
                    tkz0_t     &tkz0    ,
                    tkd0_t     &tkd0    ,
                    chi2rphi_t &chi2rphi,
                    chi2rz_t   &chi2rz  ,
                    bendChi2_t &bendChi2,
                    hit_t      &hit     ,
                    trackMVA_t &trackMVA,
                    extraMVA_t &extraMVA,
                    valid_t    &valid   ){
    unsigned int lo = 0;
    unsigned int len = 0;
    len=rinv_t    ::width; bit_copy(in, rinv    , lo); lo += len;
    len=tkphi_t   ::width; bit_copy(in, tkphi   , lo); lo += len;
    len=tanlam_t  ::width; bit_copy(in, tanlam  , lo); lo += len;
    len=tkz0_t    ::width; bit_copy(in, tkz0    , lo); lo += len;
    len=tkd0_t    ::width; bit_copy(in, tkd0    , lo); lo += len;
    len=chi2rphi_t::width; bit_copy(in, chi2rphi, lo); lo += len;
    len=chi2rz_t  ::width; bit_copy(in, chi2rz  , lo); lo += len;
    len=bendChi2_t::width; bit_copy(in, bendChi2, lo); lo += len;
    len=hit_t     ::width; bit_copy(in, hit     , lo); lo += len;
    len=trackMVA_t::width; bit_copy(in, trackMVA, lo); lo += len;
    len=extraMVA_t::width; bit_copy(in, extraMVA, lo); lo += len;
    len=valid_t   ::width; bit_copy(in, valid   , lo); lo += len;
}

void pack_pf_track(ap_uint<64> &tk,
                   pt_t     pf_pt   ,
                   pt_t     pf_pterr,
                   etaphi_t pf_eta  ,
                   etaphi_t pf_phi  ,
                   z0_t     pf_z0   ,
                   bool     pf_TightQuality){
    bool debug=false;
    if(debug){
        
        std::cout << "pack_pf_track" << std::endl;
        std::cout << "  pf_pt           " << pf_pt          .to_string(16) << std::endl;
        std::cout << "  pf_pterr        " << pf_pterr       .to_string(16) << std::endl;
        std::cout << "  pf_eta          " << pf_eta         .to_string(16) << std::endl;
        std::cout << "  pf_phi          " << pf_phi         .to_string(16) << std::endl;
        std::cout << "  pf_z0           " << pf_z0          .to_string(16) << std::endl;
        std::cout << "  pf_TightQuality " << pf_TightQuality << std::endl;
        
    }
    tk = (pf_TightQuality, pf_z0, pf_phi, pf_eta, pf_pterr, pf_pt);    
}

void unpack_pf_track(ap_uint<64> in,
                   pt_t     &pf_pt   ,
                   pt_t     &pf_pterr,
                   etaphi_t &pf_eta  ,
                   etaphi_t &pf_phi  ,
                   z0_t     &pf_z0   ,
                   bool     &pf_TightQuality){
    unsigned int lo = 0;
    unsigned int len = 0;
    len=pt_t    ::width; bit_copy(in, pf_pt          , lo); lo += len;
    len=pt_t    ::width; bit_copy(in, pf_pterr       , lo); lo += len;
    len=etaphi_t::width; bit_copy(in, pf_eta         , lo); lo += len;
    len=etaphi_t::width; bit_copy(in, pf_phi         , lo); lo += len;
    len=z0_t    ::width; bit_copy(in, pf_z0          , lo); lo += len;
    pf_TightQuality = in[lo];
}

template<class in_t, class out_t> 
void bit_copy(in_t in, out_t &out, int offset){
    for(int i  =out_t::width-1; i>=0; i--){
        out[i] = in[i+offset];
    }    
}

template<class pt_inv_T, class pt_T>
void init_pt_inv_table(pt_T table_out[(1<<PT_INV_TAB_SIZE)]) {
    // index is a uint from 0 to 111...11=2^PT_INV_TAB_SIZE-1 that encodes 0.000 to 0.4999999
    // resulting value pt_T is a uint from 0 to 2^16-1
    table_out[0] = (1<<PT_INV_MAX_BITS);
    for (unsigned int i = 1; i < (1<<PT_INV_TAB_SIZE); i++) {
        float invpt = float(i)/(1<<PT_INV_TAB_SIZE) * 0.5; // in 1/GeV
        table_out[i] = PF_PT_SCALE / invpt;
    }
    return;
}

template<class pt_inv_T, class pt_T>
void convert_pt(pt_inv_T inv, pt_T &pt){
    //pt_inv_T is ap_fixed<> (signed)
    //pt_T is ap_uint<> !

    // Initialize the lookup tables
#ifdef __HLS_SYN__
    bool initialized = false;
    pt_t inv_table[(1<<PT_INV_TAB_SIZE)];
#else 
    static bool initialized = false;
    static pt_t inv_table[(1<<PT_INV_TAB_SIZE)];
#endif
    if (!initialized) {
        init_pt_inv_table<pt_inv_T,pt_T>(inv_table);
        initialized = true;
    }

    // if(inv<0) inv = -inv; assume input is positive (upstream)
    urinv_t uinv = inv;

    // cutoffs at high and low pt
    if(uinv >= urinv_t(0.5)){
        pt = 2. * PF_PT_SCALE;
        return;
    } else if (uinv <= urinv_t(1./(1<<PT_INV_MAX_BITS))){
        pt=(1<<PT_INV_MAX_BITS) * PF_PT_SCALE;
        return;
    }

    ap_uint<PT_INV_TAB_SIZE> index;
    const int offset = 1; // ignore the first bit since 0b0.01111.. = 0.499.. is largest value
    #pragma unroll
    for(int i=0; i<PT_INV_TAB_SIZE; i++){
        index[PT_INV_TAB_SIZE-1-i] = uinv[urinv_t::width-1-i-offset]; //msb down to lowest
    }

    pt = inv_table[index];
}

template<class eta_T>
void init_eta_table(eta_T table_out[(1<<ETA_TAB_SIZE)]) {
    // index is a uint from 0 to 111...11=2^ETA_TAB_SIZE-1 that encodes 0.000 to 8
    // resulting value eta_T is a uint

    for (unsigned int i = 0; i < (1<<ETA_TAB_SIZE); i++) {
        // eta =  -ln(tan((pi/2 - arctan(TANLAM))/2))
        float tanlam = float(i)/(1<<ETA_TAB_SIZE) * 8.;
        float eta = -log(tan((M_PI/2 - atan(tanlam))/2));
        table_out[i] = eta * PF_ETAPHI_SCALE;
        // phi in -511,512 can only hold eta up to 2.23. else saturate for now
        if (eta * PF_ETAPHI_SCALE > (1<<(eta_T::width-1))-1) table_out[i] = (1<<(eta_T::width-1))-1;
    }
    return;
}

template<class tanlam_T, class eta_T>
void convert_eta(tanlam_T tanlam, eta_T &eta){
    // tanlam_T is ap_fixed<16,3>
    // eta_T is ap_int<10>

    // Initialize the lookup tables
#ifdef __HLS_SYN__
    bool initialized = false;
    etaphi_t eta_table[(1<<ETA_TAB_SIZE)];
#else 
    static bool initialized = false;
    static etaphi_t eta_table[(1<<ETA_TAB_SIZE)];
#endif
    if (!initialized) {
        init_eta_table<eta_T>(eta_table);
        initialized = true;
    }
    bool flip = false;
    if(tanlam<0){
        tanlam = -tanlam;
        flip=true;
    }
    utanlam_t utanlam = tanlam;

    ap_uint<ETA_TAB_SIZE> index;
    #pragma unroll
    for(int i=0; i<ETA_TAB_SIZE; i++){
        index[ETA_TAB_SIZE-1-i] = utanlam[utanlam_t::width-1-i]; //msb down to lowest
    }

    eta = eta_table[index];
    if(flip) eta = -eta;
}


template<class phi_T>
void init_dphi_table(phi_T table_out[(1<<DPHI_TAB_SIZE)]) {
    // index is a uint from 0 to 111...11=2^PT_INV_TAB_SIZE-1 that encodes 0.000 to 0.4999999
    // resulting value phi_T is a uint
    table_out[0] = 0;
    for (unsigned int i = 1; i < (1<<DPHI_TAB_SIZE); i++) {
        float invpt = float(i)/(1<<DPHI_TAB_SIZE) * 0.5; // in 1/GeV
        float rCurv = (1/invpt) * (129./2)/(0.735); // curv is pt * (looper radius / 2)/(looper pt)
        float x = (DETR)/(2*rCurv);
        float dPhi = atan(x / sqrt(1-x*x));
        table_out[i] = dPhi * PF_ETAPHI_SCALE;
        // overflow guard. shouldn't happen
        if( dPhi * PF_ETAPHI_SCALE >  (1<<(phi_T::width-1))-1)  table_out[i] = (1<<(phi_T::width-1))-1;
    }
    return;
}

template<class pt_inv_T, class phi_T> 
void convert_dphi(pt_inv_T inv, phi_T &dphi){

    // Initialize the lookup tables
#ifdef __HLS_SYN__
    bool initialized = false;
    phi_T dphi_table[(1<<DPHI_TAB_SIZE)];
#else 
    static bool initialized = false;
    static phi_T dphi_table[(1<<DPHI_TAB_SIZE)];
#endif
    if (!initialized) {
        init_dphi_table<phi_T>(dphi_table);
        initialized = true;
    }

    if(inv<0) inv = -inv;
    urinv_t uinv = inv;

    // cutoffs at high and low pt
    if(uinv >= urinv_t(0.5)){
        dphi = 0.395 * PF_ETAPHI_SCALE; // low-pt: curv @ 2 GeV is 0.395
        return;
    } else if (uinv <= urinv_t(1./(1<<PT_INV_MAX_BITS))){
        dphi=0; // high-pt = straight track
        return;
    }

    ap_uint<DPHI_TAB_SIZE> index;
    const int offset = 1; // ignore the first bit since 0b0.01111.. = 0.499.. is largest value
    #pragma unroll
    for(int i=0; i<PT_INV_TAB_SIZE; i++){
        index[DPHI_TAB_SIZE-1-i] = uinv[urinv_t::width-1-i-offset]; //msb down to lowest
    }

    dphi = dphi_table[index];
}


void propagate_tanlam(tkz0_t z0, tanlam_t tanlam, tanlam_t &tanlam_at_det){
    // using a helper type to handle necessary range
    zdet_t z_at_det = z0 + zdet_t(tanlam) * DETR;
    if(z_at_det < DETZ){
        tanlam_at_det = tanlam + tanlam_help_t(z0) * tanlam_help_t(1/DETR);
    } else {
        // use power expansion for now: 1/(1+z0/z_det) -> 1-z0/z_det
        tanlam_at_det = tanlam * (1 - tanlam_help_t(z0) * tanlam_help_t(1/DETZ));
    }
}

/*
void track_propagate(pt_t pt, etaphi_t eta, etaphi_t phi, z0_t z0,etaphi_t &eta_calo, etaphi_t &phi_calo){
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

    //inputs 
    tanlam_t  tanlam;
    tanlam_t  utanlam=tanlam;
    if(tanlam<0) utanlam = -tanlam;

    zdet_t z_at_det = z0 + z_at_det(tanlam) * DETR;
    tanlam_t  tanlam_at_det;
    if(z_at_det < DETZ){
        tanlam_at_det = tanlam + tanlam_help_t(z0) * tanlam_help_t(1/DETR);
    } else {
        // use power expansion for now: 1/(1+z0/z_det) -> 1-z0/z_det
        tanlam_at_det = tanlam * (1 - tanlam_help_t(z0) * tanlam_help_t(1/DETZ));
    }

    convert_dphi(pt_inv, dphi);
    

    // float z_at_det = z0 + 129 * tanlam;
    // float tanlam_det;
    // if( z_at_det < 320.9 ){
    //     tanlam_det = z0/R_det + tanlam;
    // } else {
    //     tanlam_det = z0/(Z_det-z0) * tanlam;
    // }
    // // compute dphi in bins of curvature
    // float rC = UNIT / ptinv;
    // float phi = atan( 2*rC / R_det * sqrt(1-pow(R_det,2)/(4*pow(rC,2))) );

    // placeholder
    // eta_calo=eta;
    // phi_calo=phi;
}
*/
