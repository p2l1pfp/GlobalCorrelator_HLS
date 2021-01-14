#ifndef FIRMWARE_DATA_H
#define FIRMWARE_DATA_H

#include <ap_int.h>
#include <cassert>

typedef ap_int<16> pt_t;
typedef ap_int<10>  eta_t;
typedef ap_int<10>  phi_t;
typedef ap_int<12>  glbeta_t;
typedef ap_int<11>  glbphi_t;
typedef ap_int<5>  vtx_t;
typedef ap_uint<3>  particleid_t;
typedef ap_int<10> z0_t;  // 40cm / 0.1
typedef ap_uint<9> puppiWgt_t; // 256 = 1.0 
typedef ap_uint<14> tk2em_dr_t;
typedef ap_uint<14> tk2calo_dr_t;
typedef ap_uint<10> em2calo_dr_t;
typedef ap_uint<13> tk2calo_dq_t;

enum PID { PID_Charged=0, PID_Neutral=1, PID_Photon=2, PID_Electron=3, PID_Muon=4 };

// DEFINE MULTIPLICITIES
#if defined(REG_HGCal)
    #define NTRACK 30
    #define NCALO 20
    #define NMU 4
    #define NSELCALO 20
    #define NALLNEUTRALS NSELCALO
    // dummy
    #define NEMCALO 1
    #define NPHOTON NEMCALO
    // not used but must be there because used in header files
    #define NNEUTRALS 1
//--------------------------------
#elif defined(REG_HGCalNoTK)
    #define NCALO 12
    #define NNEUTRALS 8
    #define NALLNEUTRALS NCALO
    // dummy
    #define NMU 1
    #define NTRACK 1
    #define NEMCALO 1
    #define NPHOTON NEMCALO
    #define NSELCALO 1
//--------------------------------
#elif defined(REG_HF)
    #define NCALO 18
    #define NNEUTRALS 10
    #define NALLNEUTRALS NCALO
    // dummy
    #define NMU 1
    #define NTRACK 1
    #define NEMCALO 1
    #define NPHOTON NEMCALO
    #define NSELCALO 1
//--------------------------------
#else // BARREL
   #ifndef REG_Barrel
     #ifndef CMSSW_GIT_HASH
       #warning "No region defined, assuming it's barrel (#define REG_Barrel to suppress this)"
     #endif
   #endif
   #if defined(BOARD_MP7)
       #warning "MP7 NOT SUPPORTED ANYMORE"
       #define NTRACK 14
       #define NCALO 10
       #define NMU 2
       #define NEMCALO 10
       #define NPHOTON NEMCALO
       #define NSELCALO 10
       #define NALLNEUTRALS (NPHOTON+NSELCALO)
       #define NNEUTRALS 15
   #elif defined(BOARD_CTP7)
       #error "NOT SUPPORTED ANYMORE"
   #elif defined(BOARD_KU15P)
       #define NTRACK 14
       #define NCALO 10
       #define NMU 2
       #define NEMCALO 10
       #define NPHOTON NEMCALO
       #define NSELCALO 10
       #define NALLNEUTRALS (NPHOTON+NSELCALO)
       #define NNEUTRALS 15
   #elif defined(BOARD_VCU118)
       #define NTRACK 22
       #define NCALO 15
       #define NEMCALO 13
       #define NMU 2
       #define NPHOTON NEMCALO
       #define NSELCALO 10
       #define NALLNEUTRALS (NPHOTON+NSELCALO)
       #define NNEUTRALS 25
   #else
       #define NTRACK 22
       #define NCALO 15
       #define NEMCALO 13
       #define NMU 2
       #define NPHOTON NEMCALO
       #define NSELCALO 10
       #define NALLNEUTRALS (NPHOTON+NSELCALO)
       #define NNEUTRALS 25
   #endif

#endif // region

#if defined(BOARD_MP7)
    #define PACKING_DATA_SIZE 32
    #define PACKING_NCHANN    72
#elif defined(BOARD_KU15P)
    #define PACKING_DATA_SIZE 64
    #define PACKING_NCHANN    42
#elif defined(BOARD_VCU118)
    #define PACKING_DATA_SIZE 64
    #define PACKING_NCHANN    120
#elif defined(BOARD_APD1)
    #define PACKING_DATA_SIZE 64
    #define PACKING_NCHANN    96
#endif


template<int N> struct ct_log2_ceil { enum { value = ct_log2_ceil<(N/2)+(N%2)>::value + 1 }; };
template<> struct ct_log2_ceil<2> { enum { value = 1 }; };
template<> struct ct_log2_ceil<1> { enum { value = 0 }; };


struct CaloObj {
	pt_t hwPt;
	eta_t hwEta; // relative to the region center, at calo
	phi_t hwPhi; // relative to the region center, at calo
};
struct HadCaloObj : public CaloObj {
	pt_t hwEmPt;
   	bool hwIsEM;
};
inline void clear(HadCaloObj & c) {
    c.hwPt = 0; c.hwEta = 0; c.hwPhi = 0; c.hwEmPt = 0; c.hwIsEM = 0; 
}

struct EmCaloObj {
	pt_t hwPt, hwPtErr;
	eta_t hwEta; // relative to the region center, at calo
	phi_t hwPhi; // relative to the region center, at calo
};
inline void clear(EmCaloObj & c) {
    c.hwPt = 0; c.hwPtErr = 0; c.hwEta = 0; c.hwPhi = 0; 
}



struct TkObj {
	pt_t hwPt, hwPtErr;
	eta_t hwEta; // relative to the region center, at calo
	phi_t hwPhi; // relative to the region center, at calo
	bool hwCharge; // 1 = positive, 0 = negative
	z0_t hwZ0;
	bool hwTightQuality;
};
inline void clear(TkObj & c) {
    c.hwPt = 0; c.hwPtErr = 0; c.hwEta = 0; c.hwPhi = 0; c.hwZ0 = 0; c.hwCharge = 0; c.hwTightQuality = 0;
}

struct MuObj {
	pt_t hwPt, hwPtErr;
	eta_t hwEta; // relative to the region center, at calo
	phi_t hwPhi; // relative to the region center, at calo
};
inline void clear(MuObj & c) {
    c.hwPt = 0; c.hwPtErr = 0; c.hwEta = 0; c.hwPhi = 0; 
}


struct PFChargedObj {
	pt_t hwPt;
	eta_t hwEta; // relative to the region center, at calo
	phi_t hwPhi; // relative to the region center, at calo
	particleid_t hwId;
	z0_t hwZ0;
};
inline void clear(PFChargedObj & c) {
    c.hwPt = 0; c.hwEta = 0; c.hwPhi = 0; c.hwId = 0; c.hwZ0 = 0;
}

struct PFNeutralObj {
	pt_t hwPt;
	eta_t hwEta; // relative to the region center, at calo
	phi_t hwPhi; // relative to the region center, at calo
	particleid_t hwId;
};

inline void clear(PFNeutralObj & c) {
    c.hwPt = 0; c.hwEta = 0; c.hwPhi = 0; c.hwId = 0; 
}

struct PuppiObj {
	pt_t hwPt;
	eta_t hwEta; // relative to the region center, at calo
	phi_t hwPhi; // relative to the region center, at calo
	particleid_t hwId;
	ap_uint<12> hwData;

        inline z0_t hwZ0() const { 
            #ifndef __SYNTHESIS__
            assert(hwId == PID_Charged || hwId == PID_Electron || hwId == PID_Muon);
            #endif
            return z0_t(hwData(9,0)); 
        }
        inline void setHwZ0(z0_t z0) { 
            #ifndef __SYNTHESIS__
            assert(hwId == PID_Charged || hwId == PID_Electron || hwId == PID_Muon);
            #endif
            hwData(9,0) = z0(9,0); 
        }
        inline puppiWgt_t hwPuppiW() const { 
            #ifndef __SYNTHESIS__
            assert(hwId == PID_Neutral || hwId == PID_Photon);
            #endif
            return puppiWgt_t(hwData(8,0)); 
        }
        inline void setHwPuppiW(puppiWgt_t w) { 
            #ifndef __SYNTHESIS__
            assert(hwId == PID_Neutral || hwId == PID_Photon);
            #endif
            hwData(8,0) = w(8,0); 
        }
};
inline void clear(PuppiObj & c) {
    c.hwPt = 0; c.hwEta = 0; c.hwPhi = 0; c.hwId = 0; c.hwData = 0;
}
inline void fill(PuppiObj &out, const PFChargedObj &src) {
    out.hwEta = src.hwEta;
    out.hwPhi = src.hwPhi;
    out.hwId  = src.hwId;
    out.hwPt  = src.hwPt;
    out.hwData = 0;
    out.setHwZ0(src.hwZ0);
}
inline void fill(PuppiObj &out, const PFNeutralObj &src, pt_t puppiPt, puppiWgt_t puppiWgt) {
    out.hwEta = src.hwEta;
    out.hwPhi = src.hwPhi;
    out.hwId  = src.hwId;
    out.hwPt  = puppiPt;
    out.hwData = 0;
    out.setHwPuppiW(puppiWgt);
}
inline void fill(PuppiObj &out, const HadCaloObj &src, pt_t puppiPt, puppiWgt_t puppiWgt) {
    out.hwEta = src.hwEta;
    out.hwPhi = src.hwPhi;
    out.hwId  = src.hwIsEM ? PID_Photon : PID_Neutral;
    out.hwPt  = puppiPt;
    out.hwData = 0;
    out.setHwPuppiW(puppiWgt);
}





//TMUX
#define NETA_TMUX 2
#define NPHI_TMUX 1
/* #define TMUX_IN 36 */
/* #define TMUX_OUT 18 */
#define TMUX_IN 18
#define TMUX_OUT 6
#define NTRACK_TMUX (NTRACK*TMUX_OUT*NETA_TMUX*NPHI_TMUX)
#define NCALO_TMUX (NCALO*TMUX_OUT*NETA_TMUX*NPHI_TMUX)
#define NEMCALO_TMUX (NEMCALO*TMUX_OUT*NETA_TMUX*NPHI_TMUX)
#define NMU_TMUX (NMU*TMUX_OUT*NETA_TMUX*NPHI_TMUX)

#endif
