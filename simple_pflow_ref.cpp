#include "simple_pflow.h"

void simple_pflow_parallel_ref(CaloObj calo[NCALO], TkObj track[NTRACK], PFObj out[NPF]) {
	// constants
	const etaphi_t BOX_SIZE = 81; // ETAPHI_SCALE * 0.2 * std::sqrt(M_PI/4);
	const pt_t     TKPT_MAX = 80; // 20 * PT_SCALE;
	const int      DR2MAX   = 2101;

	// initialize sum track pt
	pt_t calo_sumtk[NCALO], calo_subpt[NCALO];
	int  calo_sumtkErr2[NCALO];
	for (int ic = 0; ic < NCALO; ++ic) { calo_sumtk[ic] = 0;  calo_sumtkErr2[ic] = 0;}

	// initialize good track bit
	bool track_good[NTRACK];
	for (int it = 0; it < NTRACK; ++it) { track_good[it] = (track[it].hwPt < TKPT_MAX); }

	// initialize output
	for (int ipf = 0; ipf < NPF; ++ipf) { out[ipf].hwPt = 0; }

	// for each track, find the closest calo
	for (int it = 0; it < NTRACK; ++it) {
		if (track[it].hwPt > 0) {
			pt_t caloPtMin = track[it].hwPt - 2*(track[it].hwPtErr);
			int  drmin = DR2MAX, ibest = -1;
			for (int ic = 0; ic < NCALO; ++ic) {
				if (calo[ic].hwPt <= caloPtMin) continue;
				int dr = dr2_int(track[it].hwEta, track[it].hwPhi, calo[ic].hwEta, calo[ic].hwPhi);
				if (dr < drmin) { drmin = dr; ibest = ic; }
			}
			if (ibest != -1) {
//#ifndef __SYNTHESIS__
//				printf("ref: track %2d pt % 7d linked to calo %2d pt % 7d\n", it, int(track[it].hwPt), ibest, int(calo[ibest].hwPt));
//#endif
				track_good[it] = 1;
				calo_sumtk[ibest]    += track[it].hwPt;
				calo_sumtkErr2[ibest] += track[it].hwPtErr*track[it].hwPtErr;
			}
		}
	}


	for (int ic = 0; ic < NCALO; ++ic) {
		if (calo_sumtk[ic] > 0) {
			pt_t ptdiff = calo[ic].hwPt - calo_sumtk[ic];
			if (ptdiff > 0 && ptdiff*ptdiff > 4*calo_sumtkErr2[ic]) {
				calo_subpt[ic] = ptdiff;
			} else {
				calo_subpt[ic] = 0;
			}
		} else {
			calo_subpt[ic] = calo[ic].hwPt;
		}
	}

	// copy out charged hadrons
	for (int it = 0, ipf = 0; it < NTRACK; ++it, ++ipf) {
		if (track_good[it]) {
			out[ipf].hwPt = track[it].hwPt;
			out[ipf].hwEta = track[it].hwEta;
			out[ipf].hwPhi = track[it].hwPhi;
			out[ipf].hwId  = PID_Charged;
		}
	}

	// copy out neutral hadrons
	for (int ic = 0, ipf = NTRACK; ic < NCALO; ++ic, ++ipf) {
		if (calo_subpt[ic] > 0) {
			out[ipf].hwPt  = calo_subpt[ic];
			out[ipf].hwEta = calo[ic].hwEta;
			out[ipf].hwPhi = calo[ic].hwPhi;
			out[ipf].hwId  = PID_Neutral;
		}
	}
}