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
void _spfph_tkerr2(TkObj track[NTRACK], int tkerr2[NTRACK]) {
	for (int it = 0; it < NTRACK; ++it) {
		tkerr2[it] = track[it].hwPtErr * track[it].hwPtErr;
	}
}
void _spfph_sumtk(TkObj track[NTRACK], int tkerr2[NTRACK], bool calo_track_link_bit[NTRACK][NCALO], pt_t sumtk[NCALO], int sumtkerr2[NCALO]) {
	for (int icalo = 0; icalo < NCALO; ++icalo) {
		pt_t sum = 0;
		int sumerr = 0;
		for (int it = 0; it < NTRACK; ++it) {
			if (calo_track_link_bit[it][icalo]) { sum += track[it].hwPt; sumerr += tkerr2[it]; }
		}
		sumtk[icalo] = sum;
		sumtkerr2[icalo] = sumerr;
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

void _spfph_caloalgo(CaloObj calo[NCALO], pt_t sumtk[NCALO], int sumtkerr2[NCALO], PFObj pfout[NPF]) {
	for (int icalo = 0, ipf = NTRACK; icalo < NCALO; ++icalo, ++ipf) {
		pt_t calopt;
		if (sumtk[icalo] == 0) {
			calopt = calo[icalo].hwPt;
		} else {
			pt_t ptdiff = calo[icalo].hwPt - sumtk[icalo];
			if (ptdiff > 0 && ptdiff*ptdiff > 4*sumtkerr2[icalo]) {
				calopt = ptdiff;
			} else {
				calopt = 0;
			}
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

	int  tkerr2[NTRACK];
	#pragma HLS ARRAY_PARTITION variable=tkerr2 complete
	_spfph_tkerr2(track, tkerr2);

	pt_t sumtk[NCALO]; int sumtkerr2[NCALO];
	#pragma HLS ARRAY_PARTITION variable=sumtk complete
    #pragma HLS ARRAY_PARTITION variable=sumtkerr2 complete

	_spfph_tkalgo(track, calo_track_link_bit, out);
	_spfph_sumtk(track, tkerr2, calo_track_link_bit, sumtk, sumtkerr2);
	_spfph_caloalgo(calo, sumtk, sumtkerr2, out);
}
