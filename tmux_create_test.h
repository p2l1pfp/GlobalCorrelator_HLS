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

#include "converter/firmware/wrapper.h"


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

#define CONVERT_INPUTS 0

template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack_full_em(l1tpf_int::CaloCluster emcalo_in[N], MP7DataWord data[]) {
    EmCaloObj emcalo[N];
    std::vector<l1tpf_int::CaloCluster> emcalo_vec;
    for (unsigned int i = 0; i < N; ++i) {
        if (emcalo_in[i].hwPt != 0) emcalo_vec.push_back(emcalo_in[i]);
        else break;
    }
    dpf2fw::convert<N>(emcalo_vec, emcalo);
    for (unsigned int i = 0; i < N; ++i) {
        data[NWORDS_EMCALO*i+0+OFFS] = ( emcalo[i].hwPtErr, emcalo[i].hwPt );
        data[NWORDS_EMCALO*i+1+OFFS] = ( emcalo[i].hwPhi,   emcalo[i].hwEta );
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack_full_had(l1tpf_int::CaloCluster hadcalo_in[N], MP7DataWord data[]) {
    HadCaloObj hadcalo[N];
    std::vector<l1tpf_int::CaloCluster> hadcalo_vec;
    for (unsigned int i = 0; i < N; ++i) {
        if (hadcalo_in[i].hwPt != 0) hadcalo_vec.push_back(hadcalo_in[i]);
        else break;
    }
    dpf2fw::convert<N>(hadcalo_vec, hadcalo);
    for (unsigned int i = 0; i < N; ++i) {
        data[NWORDS_CALO*i+0+OFFS] = ( hadcalo[i].hwEmPt, hadcalo[i].hwPt );
        data[NWORDS_CALO*i+1+OFFS] = ( hadcalo[i].hwIsEM, hadcalo[i].hwPhi, hadcalo[i].hwEta );
    }
}

/* void track_convert(l1tpf_int::PropagatedTrack track_in, MP7DataWord data[NWORDS_TRACK]) { */
/*     tanlam_t tan_lambda = (M_PI/2.)-(2.*atan(exp(-1.*track_in.floatEta()))); */
/*     float phi_center = M_PI*(2.*float(int((track_in.floatPhi()+M_PI)*float(TT_NPHI_SECTORS)/(2.*M_PI)))+1.)/float(TT_NPHI_SECTORS); */
/*     tkphi_t rel_phi = (track_in.floatPhi()-phi_center)/1.026; */
/*     rinv_t r_inverse = (track_in.hwCharge ? 1. : -1.)*((0.01/kSynchrotron)*(1./track_in.floatPt())-(1./kRmax))/(((1./kRmin)-(1./kRmax))/(double(1<<14)-1.)); */
/*     tkz0_t z0scaled = double(track_in.hwZ0)/(double(l1tpf_int::InputTrack::Z0_SCALE)*1.066666); */
/*     chi2rphi_t conv_chi2 = (1<<kChi2RZSize)-1; */
/*     for (int ic = (1<<kChi2RZSize)-2; ic >=0; ic++) { */
/*         if (rzphiChi2Bins[ic] <= track_in.hwChi2) { */
/*             conv_chi2 = ic; */
/*             break; */
/*         } */
/*     } */
/*     hit_t hit_pattern = track_in.hwStubs; */
/*     ap_uint<5> split0 = tan_lambda.range(4,0); // 32 - 15 - 12 */
/*     ap_uint<11> split01 = tan_lambda.range(15,5); // 16 - 5 */
/*     tkd0_t d0dummy = track_in.floatPhi()*0.09; //dummy data */
/*     ap_uint<9> split12 = d0dummy.range(8,0); // 32 - 11 - 12 */
/*     ap_uint<4> split2 = d0dummy.range(12,8); // 13 - 9 */
/*     data[2] = ( split0, ap_uint<kPhiSize>(rel_phi.range(kPhiSize-1,0)), ap_uint<kPtSize+1>(r_inverse.range(kPtSize,0)) ); */
/*     data[1] = ( split12, ap_uint<kZ0Size>(z0scaled.range(kZ0Size-1,0)), split01 ); */
/*     data[0] = ( valid_t(1), extraMVA_t(int(track_in.floatPhi())), trackMVA_t(2-int(track_in.floatPhi())), hit_t(track_in.hwStubs), bendChi2_t(3+int(track_in.floatPhi())), chi2rz_t(conv_chi2-int(track_in.floatPhi())), conv_chi2, split2 ); //dummy data here too */
/* } */

void track_convert(l1tpf_int::PropagatedTrack track_in, MP7DataWord data[NWORDS_TRACK]) {

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

    // dummy data
    tkd0_t d0dummy = track_in.floatPhi()*0.09; 
    trackMVA_t trackMVA = 2-int(track_in.floatPhi());
    extraMVA_t extraMVA = int(track_in.floatPhi());
    valid_t    valid    = 1;
    chi2rz_t chi2rz = 2;
    bendChi2_t bendChi2 = 3;

    ap_uint<96> packed;
    pack_L1T_track(packed, r_inverse, rel_phi, tan_lambda, z0scaled, d0dummy, conv_chi2, chi2rz, bendChi2, hit_pattern, trackMVA, extraMVA, valid);
    data[2] = packed(95,64);
    data[1] = packed(63,32);
    data[0] = packed(31, 0);
}

template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack_full(l1tpf_int::PropagatedTrack track_in[N], MP7DataWord data[]) {
    /*TkObj track[N];
    std::vector<l1tpf_int::PropagatedTrack> track_vec;
    for (unsigned int i = 0; i < N; ++i) {
        if (track_in[i].hwPt != 0) track_vec.push_back(track_in[i]);
        else break;
    }
    dpf2fw::convert<N>(track_vec, track);
    for (unsigned int i = 0; i < N; ++i) {
        data[2*i+0+OFFS] = ( track[i].hwPtErr, track[i].hwPt );
        data[2*i+1+OFFS] = ( track[i].hwZ0, track[i].hwPhi, track[i].hwEta );
        data[2*i+1+OFFS][30] = track[i].hwTightQuality;
    }*/
    for (unsigned int i = 0; i < N; ++i) {
        if (track_in[i].hwPt != 0) {
            MP7DataWord tmpdata[NWORDS_TRACK];
            track_convert(track_in[i], tmpdata);
            data[NWORDS_TRACK*i+2+OFFS] = tmpdata[2];
            data[NWORDS_TRACK*i+1+OFFS] = tmpdata[1];
            data[NWORDS_TRACK*i+0+OFFS] = tmpdata[0]; //dummy data here too
        } else {
            data[NWORDS_TRACK*i+2+OFFS] = 0;
            data[NWORDS_TRACK*i+1+OFFS] = 0;
            data[NWORDS_TRACK*i+0+OFFS] = 0;
        }
        /* if(CONVERT_INPUTS){ */
        /*     ap_uint<96> in = (data[NWORDS_TRACK*i+2+OFFS], */
        /*                       data[NWORDS_TRACK*i+1+OFFS], */
        /*                       data[NWORDS_TRACK*i+0+OFFS]); */
        /*     ap_uint<64> out; */
        /*     pf_input_track_conv_hw(in, out); */
        /*     data[NWORDS_TRACK*i+2+OFFS] = 0; */
        /*     data[NWORDS_TRACK*i+1+OFFS] = out(63,32); */
        /*     data[NWORDS_TRACK*i+0+OFFS] = out(31, 0); */
        /* } */
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack_full(l1tpf_int::Muon mu_in[N], MP7DataWord data[]) {
    MuObj mu[N];
    std::vector<l1tpf_int::Muon> mu_vec;
    for (unsigned int i = 0; i < N; ++i) {
        if (mu_in[i].hwPt != 0) mu_vec.push_back(mu_in[i]);
        else break;
    }
    dpf2fw::convert<N>(mu_vec, mu);
    for (unsigned int i = 0; i < N; ++i) {
        data[NWORDS_MU*i+0+OFFS] = ( mu[i].hwPtErr, mu[i].hwPt );
        data[NWORDS_MU*i+1+OFFS] = ( mu[i].hwPhi, mu[i].hwEta );
    }
}

void write_track_vector_to_link(std::vector<l1tpf_int::PropagatedTrack> in_vec, std::string datawords[], int offset, unsigned int link_no) {
    int index = 0;
    std::stringstream ss;
    bool held = false;
    MP7DataWord heldword = 0;
    for (auto itr = in_vec.begin(); itr != in_vec.end(); ++itr) {
        MP7DataWord tmpdata[NWORDS_TRACK] = {0, 0, 0};
        //for removing null tracks
        if (int(itr->hwPt) != 0) track_convert(*itr,tmpdata);
        if(CONVERT_INPUTS){
            ap_uint<96> in = (tmpdata[2],
                              tmpdata[1],
                              tmpdata[0]);
            ap_uint<64> out;
            pf_input_track_conv_hw(in, out, link_no); // reqs linkNo to determine phi offset
            //std::cout << "converted " << in.to_string(2) << " to " << out.to_string(2) << std::endl;
            if(in!=0) std::cout << "converted " << in.to_string(16) << " to " << out.to_string(16) << std::endl;
            tmpdata[2] = 0;
            tmpdata[1] = out(63,32);
            tmpdata[0] = out(31, 0);
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


void mp7wrapped_pack_in_full(l1tpf_int::CaloCluster emcalo[NEMCALO], l1tpf_int::CaloCluster hadcalo[NCALO], l1tpf_int::PropagatedTrack track[NTRACK], l1tpf_int::Muon mu[NMU], MP7DataWord data[MP7_NCHANN]) {
    // pack inputs
    assert(NWORDS_EMCALO*NEMCALO + NWORDS_TRACK*NTRACK + NWORDS_CALO*NCALO + NWORDS_MU*NMU <= MP7_NCHANN);
    constexpr unsigned int EMOFFS = NWORDS_TRACK*NTRACK;
    constexpr unsigned int HADOFFS =  NWORDS_EMCALO*NEMCALO+EMOFFS;
    constexpr unsigned int MUOFFS  = NWORDS_CALO*NCALO+HADOFFS;

    mp7_pack_full<NTRACK,0>(track, data);
    mp7_pack_full_em<NEMCALO,EMOFFS>(emcalo, data);
    mp7_pack_full_had<NCALO,HADOFFS>(hadcalo, data);
    mp7_pack_full<NMU,MUOFFS>(mu, data);
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
