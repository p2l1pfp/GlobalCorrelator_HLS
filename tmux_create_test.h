#include <cstdio>
#include <fstream>
#include <iomanip>
#include "firmware/simple_fullpfalgo.h"
#include "vertexing/firmware/simple_vtx.h"
#include "puppi/firmware/simple_puppi.h"
#include "utils/random_inputs.h"
#include "utils/DiscretePFInputs_IO.h"
#include "utils/DiscretePF2Firmware.h"
#include "utils/pattern_serializer.h"
#include "utils/test_utils.h"
#include "firmware/mp7pf_encoding.h"

#include "l1tk_types.h"
using namespace l1tk;

#include "converter/firmware/tk_input_converter.h"


#include <math.h>

#define NTEST 6
#define NLINKS_APX_GEN0 96
#define NFRAMES_APX_GEN0 3
#define NCLK_PER_BX 8
// NFRAMES_APX_GEN0 is the number of 64b words per frame. 
//   Must reserve leading 8b for header, but we simply zero first 32b
// NCLK_PER_BX is the number of frames per bx (320 mhz / 40mhz)
#define TT_NPHI_SECTORS 9

#define NLINKS_PER_TRACK 9
#define NLINKS_PER_CALO 10
#define NLINKS_PER_EMCALO 10
#define NLINKS_PER_MU 2
#define NLINKS_PER_REG (NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO+NLINKS_PER_MU)
#define NLINKS_PER_PFCH 10
#define NLINKS_PER_PFNH 10
#define NLINKS_PER_PFPHO 10
#define NLINKS_PER_PFMU 2
#define NLINKS_PER_PFREG (NLINKS_PER_PFCH+NLINKS_PER_PFNH+NLINKS_PER_PFPHO+NLINKS_PER_PFMU)

#define MAXETA_INT 243
#define MINETA_INT 0
#define MAXPHI_INT 512
#define NPHI_INT 1024

#define ETA_BUFFER 32
#define PHI_BUFFER 32

#define NWORDS_TRACK 3
#define NWORDS_CALO 2
#define NWORDS_EMCALO 2
#define NWORDS_MU 2


// This fn essentially replaces placeholder dpf2fw::convert from DiscretePF2Firmware.h
void track_convert(l1tpf_int::PropagatedTrack track_in, TkObj &track_pf, int linkNo) {
    // Temporarily must get components of L1Tk obj while waiting for the `official` word type
    tanlam_t tan_lambda = (M_PI/2.)-(2.*atan(exp(-1.*track_in.floatEta())));
    float phi_center = M_PI*(2.*float(int((track_in.floatPhi()+M_PI)*float(TT_NPHI_SECTORS)/(2.*M_PI)))+1.)/float(TT_NPHI_SECTORS);
    tkphi_t rel_phi = (track_in.floatPhi()-phi_center)/1.026;
    //rinv_t r_inverse = (track_in.hwCharge ? 1. : -1.)*((0.01/kSynchrotron)*(1./track_in.floatPt())-(1./kRmax))/(((1./kRmin)-(1./kRmax))/(double(1<<14)-1.));
    rinv_t r_inverse = (track_in.hwCharge ? 1. : -1.)*(1./track_in.floatPt());
    tkz0_t z0scaled = double(track_in.hwZ0)/(double(l1tpf_int::InputTrack::Z0_SCALE)*1.066666);
    chi2rphi_t conv_chi2 = (1<<kChi2RZSize)-1;
    for (int ic = (1<<kChi2RZSize)-2; ic >=0; ic++) {
        if (rzphiChi2Bins[ic] <= track_in.hwChi2) {
            conv_chi2 = ic;
            break;
        }
    }
    hit_t hit_pattern = track_in.hwStubs;
    // dummy data for now
    tkd0_t d0dummy = track_in.floatPhi()*0.09; 
    trackMVA_t trackMVA = 2-int(track_in.floatPhi());
    extraMVA_t extraMVA = int(track_in.floatPhi());
    valid_t    valid    = 1;
    chi2rz_t chi2rz = 2;
    bendChi2_t bendChi2 = 3;
    
    // the converter operates on words for now, so need to pack, convert, unpack...
    //   (should be cleaned up)
    ap_uint<96> in_packed;
    ap_uint<64> pf_packed;
    pack_L1T_track(in_packed, r_inverse, rel_phi, tan_lambda, z0scaled, d0dummy, conv_chi2, chi2rz, bendChi2, hit_pattern, trackMVA, extraMVA, valid);
    numlink_t nlink=linkNo;
    pf_input_track_conv_hw(in_packed, pf_packed, nlink);
    pt_t     pf_pt   ;
    pt_t     pf_pterr;
    etaphi_t pf_eta  ;
    etaphi_t pf_phi  ;
    z0_t     pf_z0   ;
    bool     pf_TightQuality;
    unpack_pf_track(pf_packed, pf_pt, pf_pterr, pf_eta, pf_phi, pf_z0, pf_TightQuality);

    // fill components into TkObj
    track_pf.hwPt = pf_pt; 
    track_pf.hwPtErr = pf_pterr;
    track_pf.hwEta = pf_eta;
    track_pf.hwPhi = pf_phi;
    track_pf.hwZ0 = pf_z0;
    track_pf.hwTightQuality = pf_z0;
}

void tp_track_to_words(l1tpf_int::PropagatedTrack track_in, MP7DataWord data[NWORDS_TRACK]) {
    // Temporarily must get components of L1Tk obj while waiting for the `official` word type
    tanlam_t tan_lambda = (M_PI/2.)-(2.*atan(exp(-1.*track_in.floatEta())));
    float phi_center = M_PI*(2.*float(int((track_in.floatPhi()+M_PI)*float(TT_NPHI_SECTORS)/(2.*M_PI)))+1.)/float(TT_NPHI_SECTORS);
    tkphi_t rel_phi = (track_in.floatPhi()-phi_center)/1.026;
    //rinv_t r_inverse = (track_in.hwCharge ? 1. : -1.)*((0.01/kSynchrotron)*(1./track_in.floatPt())-(1./kRmax))/(((1./kRmin)-(1./kRmax))/(double(1<<14)-1.));
    rinv_t r_inverse = (track_in.hwCharge ? 1. : -1.)*(1./track_in.floatPt());
    tkz0_t z0scaled = double(track_in.hwZ0)/(double(l1tpf_int::InputTrack::Z0_SCALE)*1.066666);
    chi2rphi_t conv_chi2 = (1<<kChi2RZSize)-1;
    for (int ic = (1<<kChi2RZSize)-2; ic >=0; ic++) {
        if (rzphiChi2Bins[ic] <= track_in.hwChi2) {
            conv_chi2 = ic;
            break;
        }
    }
    hit_t hit_pattern = track_in.hwStubs;
    // dummy data for now
    tkd0_t d0dummy = track_in.floatPhi()*0.09; 
    trackMVA_t trackMVA = 2-int(track_in.floatPhi());
    extraMVA_t extraMVA = int(track_in.floatPhi());
    valid_t    valid    = 1;
    chi2rz_t chi2rz = 2;
    bendChi2_t bendChi2 = 3;

    ap_uint<96> in_packed;
    pack_L1T_track(in_packed, r_inverse, rel_phi, tan_lambda, z0scaled, d0dummy, conv_chi2, chi2rz, bendChi2, hit_pattern, trackMVA, extraMVA, valid);
    data[2] = in_packed(95,64);
    data[1] = in_packed(63,32);
    data[0] = in_packed(31, 0);
}


void write_track_vector_to_link(std::vector<TkObj> in_vec, std::string datawords[], int offset) {
    int index = 0;
    std::stringstream ss;
    bool held = false;
    MP7DataWord heldword = 0;
    for (auto itr = in_vec.begin(); itr != in_vec.end(); ++itr) {
        MP7DataWord tmpdata[NWORDS_TRACK] = {0, 0, 0};
        //for removing null tracks
        if (int(itr->hwPt) != 0){
            tmpdata[0] = ( itr->hwPtErr, itr->hwPt );
            tmpdata[1] = ( itr->hwZ0, itr->hwPhi, itr->hwEta );
            tmpdata[1][30] = itr->hwTightQuality;
            tmpdata[2] = 0;
        }
        ss.str("");
        ss << "0x";
        if (!held) {
            ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[1]);
            ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[0]);
            datawords[offset+index] = ss.str();
            index++;
            heldword = tmpdata[2];
            held = true;
        } else {
            ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[0]);
            ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (heldword);
            datawords[offset+index] = ss.str();
            index++;
            ss.str("");
            ss << "0x";
            ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[2]);
            ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[1]);
            datawords[offset+index] = ss.str();
            index++;
            held = false;
        }
    }
    if (held) {
        ss.str("");
        ss << "0x00000000" << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (heldword);
        datawords[offset+index] = ss.str();
        index++;
    }
}
void write_track_vector_to_link(std::vector<l1tpf_int::PropagatedTrack> in_vec, std::string datawords[], int offset) {
    int index = 0;
    std::stringstream ss;
    bool held = false;
    MP7DataWord heldword = 0;
    for (auto itr = in_vec.begin(); itr != in_vec.end(); ++itr) {
        MP7DataWord tmpdata[NWORDS_TRACK] = {0, 0, 0};
        //for removing null tracks
        if (int(itr->hwPt) != 0) tp_track_to_words(*itr,tmpdata);
        ss.str("");
        ss << "0x";
        if (!held) {
            ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[1]);
            ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[0]);
            datawords[offset+index] = ss.str();
            index++;
            heldword = tmpdata[2];
            held = true;
        } else {
            ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[0]);
            ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (heldword);
            datawords[offset+index] = ss.str();
            index++;
            ss.str("");
            ss << "0x";
            ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[2]);
            ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[1]);
            datawords[offset+index] = ss.str();
            index++;
            held = false;
        }
    }
    if (held) {
        ss.str("");
        ss << "0x00000000" << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (heldword);
        datawords[offset+index] = ss.str();
        index++;
    }
}

void write_emcalo_vector_to_link(std::vector<l1tpf_int::CaloCluster> in_vec, std::string datawords[], int offset) {
    int index = 0;
    std::stringstream ss;
    for (auto itr = in_vec.begin(); itr != in_vec.end(); ++itr) {
        MP7DataWord tmpdata[NWORDS_EMCALO];
        EmCaloObj emcalo;
        dpf2fw::convert(*itr, emcalo);
        tmpdata[0] = ( emcalo.hwPtErr, emcalo.hwPt );
        tmpdata[1] = ( emcalo.hwPhi,   emcalo.hwEta );
        ss.str("");
        ss << "0x";
        ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[1]);
        ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[0]);
        datawords[offset+index] = ss.str();
        index++;
    }
}

void write_calo_vector_to_link(std::vector<l1tpf_int::CaloCluster> in_vec, std::string datawords[], int offset) {
    int index = 0;
    std::stringstream ss;
    for (auto itr = in_vec.begin(); itr != in_vec.end(); ++itr) {
        MP7DataWord tmpdata[NWORDS_CALO];
        HadCaloObj hadcalo;
        dpf2fw::convert(*itr, hadcalo);
        tmpdata[0] = ( hadcalo.hwEmPt, hadcalo.hwPt );
        tmpdata[1] = ( hadcalo.hwIsEM, hadcalo.hwPhi, hadcalo.hwEta );
        ss.str("");
        ss << "0x";
        ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[1]);
        ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[0]);
        datawords[offset+index] = ss.str();
        index++;
    }
}

void write_mu_vector_to_link(std::vector<l1tpf_int::Muon> in_vec, std::string datawords[], int offset) {
    int index = 0;
    std::stringstream ss;
    for (auto itr = in_vec.begin(); itr != in_vec.end(); ++itr) {
        MP7DataWord tmpdata[NWORDS_MU];
        MuObj mu;
        dpf2fw::convert(*itr, mu);
        tmpdata[0] = ( mu.hwPtErr, mu.hwPt );
        tmpdata[1] = ( mu.hwPhi, mu.hwEta );
        ss.str("");
        ss << "0x";
        ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[1]);
        ss << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (tmpdata[0]);
        datawords[offset+index] = ss.str();
        index++;
    }
}

void mp7wrapped_pack_in_reorder(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], MP7DataWord data[MP7_NCHANN]) {
    // pack inputs
    assert(2*NEMCALO + 2*NTRACK + 2*NCALO + 2*NMU <= MP7_NCHANN);
    constexpr unsigned int EMOFFS = 2*NTRACK;
    constexpr unsigned int HADOFFS = 2*NEMCALO+EMOFFS;
    constexpr unsigned int MUOFFS  = 2*NCALO+HADOFFS;

    mp7_pack<NTRACK,0>(track, data);
    mp7_pack<NEMCALO,EMOFFS>(emcalo, data);
    mp7_pack<NCALO,HADOFFS>(hadcalo, data);
    mp7_pack<NMU,MUOFFS>(mu, data);
}

template<typename T, int NIn, int NOut>
void ptsort_out(T in[NIn], T out[NOut]) {
    T tmp[NOut];

    for (int iout = 0; iout < NOut; ++iout) {
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

void mp7wrapped_pfalgo3_full_sort(MP7DataWord input[MP7_NCHANN], MP7DataWord output[MP7_NCHANN], z0_t Z0) {
    
    PFChargedObj pfch[NTRACK]; PFNeutralObj pfpho[NPHOTON]; PFNeutralObj pfne[NSELCALO]; PFChargedObj pfmu[NMU];
    PFChargedObj pfch_sort[NTRACK]; PFNeutralObj pfpho_sort[NPHOTON]; PFNeutralObj pfne_sort[NSELCALO]; PFChargedObj pfmu_sort[NMU];

    mp7wrapped_pfalgo3_full(input, output, Z0);
    mp7wrapped_unpack_out_necomb(output, pfch, pfpho, pfne, pfmu);
    
    ptsort_out<PFChargedObj,NTRACK,NTRACK>(pfch, pfch_sort);
    ptsort_out<PFNeutralObj,NPHOTON,NPHOTON>(pfpho, pfpho_sort);
    ptsort_out<PFNeutralObj,NSELCALO,NSELCALO>(pfne, pfne_sort);
    ptsort_out<PFChargedObj,NMU,NMU>(pfmu, pfmu_sort);

    mp7wrapped_pack_out(pfch_sort, pfpho_sort, pfne_sort, pfmu_sort, output);

}

bool isInPhiRegion(int test, int loBound, int hiBound, int MAXPHI=MAXPHI_INT, int MINPHI=MAXPHI_INT-NPHI_INT){
    // place all values on the circle 
    while (test <MINPHI) test += (MAXPHI-MINPHI);
    while (test>=MAXPHI) test -= (MAXPHI-MINPHI);
    while (loBound <MINPHI) loBound += (MAXPHI-MINPHI);
    while (loBound>=MAXPHI) loBound -= (MAXPHI-MINPHI);
    while (hiBound <MINPHI) hiBound += (MAXPHI-MINPHI);
    while (hiBound>=MAXPHI) hiBound -= (MAXPHI-MINPHI);
    // consider both orderings
    if (loBound <= hiBound) {
        return (test < hiBound) && (test >= loBound);
    }
    else {
        return (test < hiBound) || (test >= loBound);
    }
}
