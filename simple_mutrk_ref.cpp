#include "src/data.h"
#include "src/simple_mutrk.h"
//#include <hls_half.h>
#include <cmath>
#include <algorithm>
#include <cstdio>

template <typename T> int sqr(const T & t) { return t*t; }

void simple_mutrk_parallel_ref(MuObj mu[NMU], TkObj track[NTRACK], PFMuonObj outmu[NMU], TkObj outtrk[NTRACK]) {
	
	std::cout << "mu trk ref" << std::endl;
	
	// // constants
	// const etaphi_t BOX_SIZE = 81; // ETAPHI_SCALE * 0.2 * std::sqrt(M_PI/4);
	// const pt_t     TKPT_MAX = 80; // 20 * PT_SCALE;
	// const int      DR2MAX   = 2101;

	// // initialize sum track pt
	// pt_t calo_sumtk[NCALO], calo_subpt[NCALO];
	// int  calo_sumtkErr2[NCALO];
	// for (int ic = 0; ic < NCALO; ++ic) { calo_sumtk[ic] = 0;  calo_sumtkErr2[ic] = 0;}

	// // initialize good track bit
	// bool track_good[NTRACK];
	// for (int it = 0; it < NTRACK; ++it) { track_good[it] = (track[it].hwPt < TKPT_MAX); }

	// // initialize output
	// for (int ipf = 0; ipf < NTRACK; ++ipf) { outch[ipf].hwPt = 0; }
	// for (int ipf = 0; ipf < NCALO; ++ipf) { outne[ipf].hwPt = 0; }

	// // for each track, find the closest calo
	// for (int it = 0; it < NTRACK; ++it) {
	// 	if (track[it].hwPt > 0) {
	// 		pt_t caloPtMin = track[it].hwPt - 2*(track[it].hwPtErr);
	// 		int  drmin = DR2MAX, ibest = -1;
	// 		for (int ic = 0; ic < NCALO; ++ic) {
	// 			if (calo[ic].hwPt <= caloPtMin) continue;
	// 			int dr = dr2_int(track[it].hwEta, track[it].hwPhi, calo[ic].hwEta, calo[ic].hwPhi);
	// 			if (dr < drmin) { drmin = dr; ibest = ic; }
	// 		}
	// 		if (ibest != -1) {
	// 			track_good[it] = 1;
	// 			calo_sumtk[ibest]    += track[it].hwPt;
	// 			calo_sumtkErr2[ibest] += sqr(track[it].hwPtErr);
	// 		}
	// 	}
	// }


	// for (int ic = 0; ic < NCALO; ++ic) {
	// 	if (calo_sumtk[ic] > 0) {
	// 		pt_t ptdiff = calo[ic].hwPt - calo_sumtk[ic];
	// 		if (ptdiff > 0 && ptdiff*ptdiff > 4*calo_sumtkErr2[ic]) {
	// 			calo_subpt[ic] = ptdiff;
	// 		} else {
	// 			calo_subpt[ic] = 0;
	// 		}
	// 	} else {
	// 		calo_subpt[ic] = calo[ic].hwPt;
	// 	}
	// }

	// // copy out charged hadrons
	// for (int it = 0; it < NTRACK; ++it) {
	// 	if (track_good[it]) {
	// 		outch[it].hwPt = track[it].hwPt;
	// 		outch[it].hwEta = track[it].hwEta;
	// 		outch[it].hwPhi = track[it].hwPhi;
	// 		outch[it].hwZ0 = track[it].hwZ0;
	// 		outch[it].hwId  = PID_Charged;
	// 	}
	// }

	// // copy out neutral hadrons
	// for (int ic = 0; ic < NCALO; ++ic) {
	// 	if (calo_subpt[ic] > 0) {
	// 		outne[ic].hwPt  = calo_subpt[ic];
	// 		outne[ic].hwEta = calo[ic].hwEta;
	// 		outne[ic].hwPhi = calo[ic].hwPhi;
	// 		outne[ic].hwId  = PID_Neutral;
	// 	}
	// }
}

