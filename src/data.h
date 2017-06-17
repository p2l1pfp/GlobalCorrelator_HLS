#ifndef SIMPLE_PFLOW_DATA_H
#define SIMPLE_PFLOW_DATA_H

#include "ap_int.h"

typedef ap_int<16> pt_t;
typedef ap_int<9>  etaphi_t;
typedef ap_int<5>  vtx_t;
typedef ap_int<2>  particleid_t;
typedef ap_int<10> z0_t;  // 40cm / 0.1
		
enum PID { PID_Charged=0, PID_Neutral=1 };

#define NVTXPOW   4 //Granularity of vtx binning + 4 for Z0_Scale
#define NVTXBINS  (16 >> (NVTXPOW-4))  //Currently just 8 Vtx bins
#define NPOW 6      // Number of tracksin powers of 2, trying 64 with 15 PV
#define NALLTRACK (1 << NPOW)
#define NSECTOR 1
#define VTXPTMAX  200
#define NTRACK 8
#define NCALO 12
#define NSELCALO 15

#define NPVTRACK 7
#define PT_SCALE 4.0
#define ETAPHI_SCALE (4*180/M_PI)
#define Z0_SCALE 16

struct CaloObj {
	pt_t hwPt;
	etaphi_t hwEta, hwPhi; // relative to the region center, at calo
};
struct TkObj {
	pt_t hwPt, hwPtErr;
	etaphi_t hwEta, hwPhi; // relative to the region center, at calo
	z0_t hwZ0;
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


#endif
