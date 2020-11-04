#ifndef multififo_regionizer_h
#define multififo_regionizer_h

#include "../../firmware/data.h"
#include "../../firmware/l1pf_encoding.h"

struct GlbMuObj {
	pt_t hwPt, hwPtErr;
	glbeta_t hwEta; //global coordinates
	glbphi_t hwPhi; // 
};
inline void clear(GlbMuObj & c) {
    c.hwPt = 0; c.hwPtErr = 0; c.hwEta = 0; c.hwPhi = 0; 
}
inline ap_uint<64> l1pf_pattern_pack_one(const GlbMuObj & mu) {
    #pragma HLS inline
    ap_uint<64> data  = 0;
    data(31+32,16+32) = mu.hwPtErr;
    data(15+32, 0+32) = mu.hwPt;
    data(11, 0) = mu.hwEta;
    data(22,12) = mu.hwPhi;
    return data;
}
inline void l1pf_pattern_unpack_one(const ap_uint<64> & data, GlbMuObj & mu) {
    #pragma HLS inline
    mu.hwPt    = data(15+32, 0+32);
    mu.hwPtErr = data(31+32,16+32);
    mu.hwEta   = data(11, 0);
    mu.hwPhi   = data(22,12);
}



typedef ap_uint<64> PackedTkObj;
typedef ap_uint<64> PackedCaloObj;
typedef ap_uint<64> PackedMuObj;
inline void clear(ap_uint<64> & o) { o = 0; }

template<typename P>
inline P phiShifted(const P & t, int phi_shift) {
    #pragma HLS inline
    P ret = t;
    ret.hwPhi += phi_t(phi_shift);
    return ret;
}

#define REGIONIZERNCLOCKS 54
#define NPFREGIONS 9
#define PFREGION_PHI_SIZE 160  // size of a phi sector (in L1PF units, LSB = 0.25 degrees)
#define PFREGION_PHI_BORDER 58 // size of the phi border of a PF region (0.25 rad = 58, 0.30 rad = 69)
#define PFREGION_ETA_SIZE 230  // size of an eta sector: 1.0 rad => 229, round up to 230 be even  
#define PFREGION_ETA_BORDER 58 // size of the eta border of a PF region (0.25 rad = 58, 0.30 rad = 69)
#define PFLOWII  4

#define NTKSECTORS 9
#define NTKFIBERS  2
#define NTKFIFOS   6
#define NTKSORTED  NTRACK
#define NTKSTREAMS ((NTKSORTED+PFLOWII-1)/PFLOWII)

#define NCALOSECTORS 3
#define NCALOFIBERS  4 
#define NCALOFIFOS0 NCALOFIBERS
#define NCALOFIFOS12 (2*NCALOFIBERS)
#define NCALOFIFOS (NCALOFIFOS0+2*NCALOFIFOS12)
#define NCALOSORTED  NCALO
#define NCALOSTREAMS ((NCALOSORTED+PFLOWII-1)/PFLOWII)

#define NMUFIBERS  2
#define NMUSORTED  NMU
#define NMUSTREAMS ((NMUSORTED+PFLOWII-1)/PFLOWII)

#if defined(ROUTER_NOMERGE)
    #define NTKOUT (NTKSECTORS*NTKFIFOS)
    #define NCALOOUT (NCALOSECTORS*NCALOFIFOS)
    #define NMUOUT NPFREGIONS
#elif defined(ROUTER_NOMUX)
    #define NTKOUT NPFREGIONS
    #define NCALOOUT NPFREGIONS
    #define NMUOUT NPFREGIONS
#elif defined(ROUTER_NOSTREAM)
    #define NTKOUT NTKSORTED
    #define NCALOOUT NCALOSORTED
    #define NMUOUT NMUSORTED
#else
    #define NTKOUT NTKSTREAMS
    #define NCALOOUT NCALOSTREAMS
    #define NMUOUT NMUSTREAMS
#endif

bool tk_router(bool newevent, const PackedTkObj tracks_in[NTKSECTORS][NTKFIBERS], PackedTkObj tracks_out[NTKOUT], bool & newevent_out);
bool calo_router(bool newevent, const PackedCaloObj calo_in[NCALOSECTORS][NCALOFIBERS], PackedCaloObj calo_out[NCALOOUT], bool & newevent_out);
bool mu_router(bool newevent, const glbeta_t etaCenter, const PackedMuObj mu_in[NMUFIBERS], PackedMuObj mu_out[NMUOUT], bool & newevent_out);

#ifndef __SYNTHESIS__
#include <cstdio>

inline void printTrack(FILE *f, const TkObj & t) { 
    fprintf(f,"%3d % 4d % 4d  ", t.hwPt.to_int(), t.hwEta.to_int(), t.hwPhi.to_int()); // note no leading +'s or 0's, they confuse VHDL text parser
}
inline void printCalo(FILE *f, const HadCaloObj & t) { 
    fprintf(f,"%3d % 4d % 4d  ", t.hwPt.to_int(), t.hwEta.to_int(), t.hwPhi.to_int()); // note no leading +'s or 0's, they confuse VHDL text parser
}
inline void printMu(FILE *f, const MuObj & t) { 
    fprintf(f,"%3d % 4d % 4d  ", t.hwPt.to_int(), t.hwEta.to_int(), t.hwPhi.to_int()); // note no leading +'s or 0's, they confuse VHDL text parser
}
inline void printMu(FILE *f, const GlbMuObj & t) { 
    fprintf(f,"%3d % 4d % 4d  ", t.hwPt.to_int(), t.hwEta.to_int(), t.hwPhi.to_int()); // note no leading +'s or 0's, they confuse VHDL text parser
}

inline void printTrackShort(FILE *f, const TkObj & t) { 
   fprintf(f,"%3d%+04d ", t.hwPt.to_int(), t.hwPhi.to_int());
}
inline void printCaloShort(FILE *f, const HadCaloObj & t) { 
    fprintf(f,"%3d%+04d ", t.hwPt.to_int(), t.hwPhi.to_int());
}
#endif



#endif

