#ifndef SIMPLE_PFLOW_DATA_H
#define SIMPLE_PFLOW_DATA_H

#include "ap_int.h"

typedef ap_int<16> pt_t;
typedef ap_int<9>  etaphi_t;
typedef ap_int<5>  vtx_t;
typedef ap_uint<4>  particleid_t;
typedef ap_int<10> z0_t;  // 40cm / 0.1
		
enum PID { PID_Charged=0, PID_Neutral=1, PID_Photon=2, PID_Electron=3, PID_Muon=4 };

// VERTEXING
#define NVTXBINS  15
#define NPOW 6
#define NALLTRACK 1 << NPOW
#define NSECTOR 1
#define VTXPTMAX  200

// PF
#ifdef TESTMP7  // reduced input size to fit in a board
   #define NTRACK 7
   #define NCALO 5
   #define NMU 2
   #define NEMCALO 5
   #define NPHOTON NEMCALO
   #define NSELCALO 4
#else
   #define NTRACK 15
   #define NCALO 15
   #define NMU 4
   #define NEMCALO 15
   #define NPHOTON NEMCALO
   #define NSELCALO 10
#endif

// PUPPI & CHS
#define NPVTRACK 7


struct CaloObj {
	pt_t hwPt;
	etaphi_t hwEta, hwPhi; // relative to the region center, at calo
};
struct HadCaloObj : public CaloObj {
	pt_t hwEmPt;
   	bool hwIsEM;
};
struct EmCaloObj {
	pt_t hwPt, hwPtErr;
	etaphi_t hwEta, hwPhi; // relative to the region center, at calo
};
struct TkObj {
	pt_t hwPt, hwPtErr;
	etaphi_t hwEta, hwPhi; // relative to the region center, at calo
	z0_t hwZ0;
};
struct MuObj {
	pt_t hwPt, hwPtErr;
	etaphi_t hwEta, hwPhi; // relative to the region center, at vtx(?)
};
struct PFChargedObj {
	pt_t hwPt;
	etaphi_t hwEta, hwPhi; // relative to the region center, at calo
	particleid_t hwId;
	z0_t hwZ0;
};
struct PFNeutralObj {
	pt_t hwPt;
	etaphi_t hwEta, hwPhi; // relative to the region center, at calo
	particleid_t hwId;
};
struct VtxObj {
	pt_t  hwSumPt;
	z0_t  hwZ0;
	vtx_t mult;
	particleid_t hwId;
};

#define MP7_NCHANN 72
typedef ap_uint<32> MP7DataWord;

#endif
