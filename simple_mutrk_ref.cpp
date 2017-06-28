#include "src/data.h"
#include "src/simple_mutrk.h"
//#include <hls_half.h>
#include <cmath>
#include <algorithm>
#include <cstdio>

template <typename T> int sqr(const T & t) { return t*t; }

void simple_mutrk_parallel_ref(MuObj mu[NMU], TkObj track[NTRACK], PFMuonObj pfmuout[NMU]) {
		
	// const etaphi_t BOX_SIZE = 81; // ETAPHI_SCALE * 0.2 * std::sqrt(M_PI/4);
	const pt_t     TKPT_MAX = 80; // 20 * PT_SCALE;
	const int      DR2MAX   = 2101;

	// initialize good track bit
	bool mu_good[NMU];
	for (int im = 0; im < NMU; ++im) { mu_good[im] = (mu[im].hwPt < TKPT_MAX); }

	// initialize output
	for (int ipf = 0; ipf < NMU; ++ipf) { pfmuout[ipf].hwPt = 0; }

	// for each muon, find the closest track
	for (int im = 0; im < NMU; ++im) {
		if (mu[im].hwPt > 0) {
			pt_t tkPtMin = mu[im].hwPt - 2*(mu[im].hwPtErr);
			int  drmin = DR2MAX, ibest = -1;
			for (int it = 0; it < NTRACK; ++it) {
				if (track[it].hwPt <= tkPtMin) continue;
				int dr = dr2_int(mu[im].hwEta, mu[im].hwPhi, track[it].hwEta, track[it].hwPhi);
				if (dr < drmin) { drmin = dr; ibest = it; }
			}
			if (ibest != -1) {
				pfmuout[im].hwPt = track[ibest].hwPt;
				pfmuout[im].hwEta = track[ibest].hwEta;
				pfmuout[im].hwPhi = track[ibest].hwPhi;
				pfmuout[im].hwId  = PID_Muon;
				pfmuout[im].hwZ0 = track[ibest].hwZ0;				
			}
			else{
				pfmuout[im].hwPt  = 0;
				pfmuout[im].hwEta = 0;
				pfmuout[im].hwPhi = 0;
				pfmuout[im].hwId  = 0;
				pfmuout[im].hwZ0  = 0;				
			}
		}
	}

}

