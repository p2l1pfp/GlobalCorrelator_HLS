#ifndef FASTPUPPI_NTUPLERPRODUCER_FIRMWARE2DISCRETEPF_H
#define FASTPUPPI_NTUPLERPRODUCER_FIRMWARE2DISCRETEPF_H

/// NOTE: this include is not standalone, since the path to DiscretePFInputs is different in CMSSW & Vivado_HLS

#include "../firmware/data.h"
#include <vector>
#include <cassert>

namespace fw2dpf {

    // convert inputs from discrete to firmware
    inline void convert(const PFChargedObj & src, const l1tpf_int::PropagatedTrack & track, std::vector<l1tpf_int::PFParticle> &out) {
        l1tpf_int::PFParticle pf;
        pf.hwPt = src.hwPt;
        pf.hwEta = src.hwEta;
        pf.hwPhi = src.hwPhi;
        pf.hwVtxEta = src.hwEta; // FIXME: get from the track
        pf.hwVtxPhi = src.hwPhi; // before propagation
        pf.track = track; // FIXME: ok only as long as there is a 1-1 mapping
        pf.cluster.hwPt = 0;
        switch(src.hwId) {
            case PID_Electron: pf.hwId =  1; break;
            case PID_Muon: pf.hwId =  4; break;
            default: pf.hwId = 0; break;
        };
        pf.hwStatus = 0;
        out.push_back(pf);
    }
    inline void convert(const PFNeutralObj & src, std::vector<l1tpf_int::PFParticle> &out) {
        l1tpf_int::PFParticle pf;
        pf.hwPt = src.hwPt;
        pf.hwEta = src.hwEta;
        pf.hwPhi = src.hwPhi;
        pf.hwVtxEta = src.hwEta;
        pf.hwVtxPhi = src.hwPhi;
        pf.track.hwPt = 0;
        pf.cluster.hwPt = src.hwPt;
        switch(src.hwId) {
            case PID_Photon: pf.hwId = 3; break;
            default: pf.hwId = 2; break;
        }
        pf.hwStatus = 0;
        out.push_back(pf);
    }

    // convert inputs from discrete to firmware
    inline void convert(const TkObj & in, l1tpf_int::PropagatedTrack & out) {
        out.hwPt = in.hwPt;
        out.hwCaloPtErr = in.hwPtErr;
        out.hwEta = in.hwEta; // @calo
        out.hwPhi = in.hwPhi; // @calo
        out.hwZ0 = in.hwZ0;
    }
    inline void convert(const HadCaloObj & in, l1tpf_int::CaloCluster & out) {
        out.hwPt = in.hwPt;
        out.hwEmPt = in.hwEmPt;
        out.hwEta = in.hwEta;
        out.hwPhi = in.hwPhi;
        out.isEM = in.hwIsEM;
    }
    inline void convert(const EmCaloObj & in, l1tpf_int::CaloCluster & out) {
        out.hwPt = in.hwPt;
        out.hwPtErr = in.hwPtErr;
        out.hwEta = in.hwEta;
        out.hwPhi = in.hwPhi;
    }
    inline void convert(const MuObj & in, l1tpf_int::Muon & out) {
        out.hwPt = in.hwPt;
        out.hwEta = in.hwEta; // @calo
        out.hwPhi = in.hwPhi; // @calo
    }


    template<unsigned int NMAX, typename In>
    void convert(const In in[NMAX], std::vector<l1tpf_int::PFParticle> &out) {
        for (unsigned int i = 0; i < NMAX; ++i) {
            if (in[i].hwPt > 0) convert(in[i], out);
        }
    } 
    template<unsigned int NMAX>
    void convert(const PFChargedObj in[NMAX], std::vector<l1tpf_int::PropagatedTrack> srctracks, std::vector<l1tpf_int::PFParticle> &out) {
        for (unsigned int i = 0; i < NMAX; ++i) {
            if (in[i].hwPt > 0) {
                assert(i < srctracks.size());
                convert(in[i], srctracks[i], out);
            }
        }
    }

} // namespace

#endif
