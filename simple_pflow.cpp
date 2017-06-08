#include "simple_pflow.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

bool match_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, etaphi_t boxSize) {
	etaphi_t deta = (eta1-eta2);
	etaphi_t dphi = (phi1-phi2);
	return (deta <= boxSize && deta >= -boxSize && dphi <= boxSize && dphi >= -boxSize);
}
etaphi_t dr_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
	etaphi_t deta = eta1 > eta2 ? (eta1-eta2) : (eta2-eta1);
	etaphi_t dphi = phi1 > phi2 ? (phi1-phi2) : (phi2-phi1);
	return (deta > dphi ? deta : dphi);
}
int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
	etaphi_t deta = (eta1-eta2);
	etaphi_t dphi = (phi1-phi2);
	return deta*deta + dphi*dphi;
}
ap_uint<12> dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<12> max) {
	etaphi_t deta = (eta1-eta2);
	etaphi_t dphi = (phi1-phi2);
	int dr2 = deta*deta + dphi*dphi;
	return (dr2 < int(max) ? ap_uint<12>(dr2) : max);
}

void simple_pflow_parallel_ref(CaloObj calo[NCALO], TkObj track[NTRACK], PFObj out[NPF]) {
	// constants
	const etaphi_t BOX_SIZE = 81; // ETAPHI_SCALE * 0.2 * std::sqrt(M_PI/4);
	const pt_t     TKPT_MAX = 80; // 20 * PT_SCALE;
	const int      DR2MAX   = 2101;

	// initialize sum track pt
	pt_t calo_sumtk[NCALO], calo_sumtkErr[NCALO], calo_subpt[NCALO];
	for (int ic = 0; ic < NCALO; ++ic) { calo_sumtk[ic] = 0;  calo_sumtkErr[ic] = 0;}

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
				calo_sumtkErr[ibest] += track[it].hwPtErr;
			}
		}
	}


	for (int ic = 0; ic < NCALO; ++ic) {
		if (calo_sumtk[ic] > 0) {
			pt_t ptdiff = calo[ic].hwPt - calo_sumtk[ic];
			if (ptdiff > 2*calo_sumtkErr[ic]) {
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

void _spfph_tk2calo_link(CaloObj calo[NCALO], TkObj track[NTRACK], bool calo_track_link_bit[NTRACK][NCALO]) {
	//const ap_uint<12> DR2MAX = 2101;
	const int DR2MAX = 2101;

	for (int it = 0; it < NTRACK; ++it) {
		pt_t caloPtMin = track[it].hwPt - 2*(track[it].hwPtErr);

//		ap_uint<12> drvals[NCALO];
		int drvals[NCALO];

		R1:
		for (int icalo = 0; icalo < NCALO; ++icalo) {
			if (calo[icalo].hwPt > caloPtMin) {
				drvals[icalo] = dr2_int(track[it].hwEta, track[it].hwPhi, calo[icalo].hwEta, calo[icalo].hwPhi); //, DR2MAX);
			} else {
				drvals[icalo] = DR2MAX;
			}
			if (drvals[icalo] < DR2MAX) calo_track_link_bit[it][icalo] = 1;
			else 					    calo_track_link_bit[it][icalo] = 0;
		}
/*#ifndef __SYNTHESIS__
		printf("hwopt: track %2d pt % 7d drvals: ",it, int(track[it].hwPt));
		for (int i = 0; i < NCALO; ++i) {
			if (calo[i].hwPt > caloPtMin) {
				printf("#%02d = %5d,   ", i, drvals[i]);
			} else {
				printf("#%02d = lowpt,   ", i);
			}
		}
		printf("\n");

		printf("hwopt: track %2d pt % 7d pre-linked to: ",it, int(track[it].hwPt));
		for (int i = 0; i < NCALO; ++i) { if (calo_track_link_bit[it][i]) printf("calo %2d pt % 7d     ", i, int(calo[i].hwPt)); }
		printf("\n");
#endif */

		R2: // N^2 loop, but in theory fully parallel
		for (int i = 0; i < NCALO; ++i) {
			for (int j = 0; j < NCALO; ++j) {
				if (i < j) {
					if (drvals[i] >  drvals[j]) calo_track_link_bit[it][i] = 0;
				} else if (i > j) {
					if (drvals[i] >= drvals[j]) calo_track_link_bit[it][i] = 0;
				}
			}
		}

/*		#ifndef __SYNTHESIS__
			printf("hwopt: track %2d pt % 7d linked to: ",it, int(track[it].hwPt));
			for (int i = 0; i < NCALO; ++i) { if (calo_track_link_bit[it][i]) printf("calo %2d pt % 7d     ", i, int(calo[i].hwPt)); }
			printf("\n");
		#endif*/

	}
}
void _spfph_sumtk(TkObj track[NTRACK], bool calo_track_link_bit[NTRACK][NCALO], pt_t sumtk[NCALO], pt_t sumtkerr[NCALO]) {
	for (int icalo = 0; icalo < NCALO; ++icalo) {
		pt_t sum = 0;
		pt_t sumerr = 0;
		for (int it = 0; it < NTRACK; ++it) {
			if (calo_track_link_bit[it][icalo]) { sum += track[it].hwPt; sumerr += track[it].hwPtErr; }
		}
		sumtk[icalo] = sum;
		sumtkerr[icalo] = sumerr;
	}
}

void _spfph_tkalgo(TkObj track[NTRACK], bool calo_track_link_bit[NTRACK][NCALO], PFObj pfout[NPF]) {
	const pt_t TKPT_MAX = 80; // 20 * PT_SCALE;
	for (int it = 0; it < NTRACK; ++it) {
		bool good = (track[it].hwPt < TKPT_MAX);
		for (int icalo = 0; icalo < NCALO; ++icalo) {
			if (calo_track_link_bit[it][icalo]) good = true;
		}
		if (good) {
			pfout[it].hwPt = track[it].hwPt;
			pfout[it].hwEta = track[it].hwEta;
			pfout[it].hwPhi = track[it].hwPhi;
			pfout[it].hwId  = PID_Charged;
		} else {
			pfout[it].hwPt  = 0;
			pfout[it].hwEta = 0;
			pfout[it].hwPhi = 0;
			pfout[it].hwId  = 0;
		}
	}
}

void _spfph_caloalgo(CaloObj calo[NCALO], pt_t sumtk[NCALO], pt_t sumtkerr[NCALO], PFObj pfout[NPF]) {
	for (int icalo = 0, ipf = NTRACK; icalo < NCALO; ++icalo, ++ipf) {
		pt_t calopt;
		if (sumtk[icalo] == 0) {
			calopt = calo[icalo].hwPt;
		} else {
			pt_t ptdiff = calo[icalo].hwPt - sumtk[icalo];
			calopt = ptdiff > 2*sumtkerr[icalo] ? ptdiff : pt_t(0);
		}
		pfout[ipf].hwPt  = calopt;
		pfout[ipf].hwEta = calopt ? calo[icalo].hwEta : etaphi_t(0);
		pfout[ipf].hwPhi = calopt ? calo[icalo].hwPhi : etaphi_t(0);
		pfout[ipf].hwId  = calopt ? PID_Neutral : 0;
	}
}

void simple_pflow_parallel_hwopt(CaloObj calo[NCALO], TkObj track[NTRACK], PFObj out[NPF]) {
	#pragma HLS ARRAY_PARTITION variable=calo complete
	#pragma HLS ARRAY_PARTITION variable=track complete
	#pragma HLS ARRAY_PARTITION variable=out complete

	#pragma HLS pipeline II=5 rewind

	bool calo_track_link_bit[NTRACK][NCALO];
	#pragma HLS ARRAY_PARTITION variable=calo_track_link_bit dim=0 complete

	_spfph_tk2calo_link(calo, track, calo_track_link_bit);

	pt_t sumtk[NCALO]; pt_t sumtkerr[NCALO];
	#pragma HLS ARRAY_PARTITION variable=sumtk complete
    #pragma HLS ARRAY_PARTITION variable=sumtkerr complete

	_spfph_tkalgo(track, calo_track_link_bit, out);
	_spfph_sumtk(track, calo_track_link_bit, sumtk, sumtkerr);
	_spfph_caloalgo(calo, sumtk, sumtkerr, out);
}
