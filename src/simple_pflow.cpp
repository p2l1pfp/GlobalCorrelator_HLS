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
		tkerr2[it] = (track[it].hwPtErr * track[it].hwPtErr) << 2;
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

void _spfph_tkalgo(TkObj track[NTRACK], bool calo_track_link_bit[NTRACK][NCALO], PFChargedObj pfout[NPF]) {
	const pt_t TKPT_MAX = 80; // 20 * PT_SCALE;
	for (int it = 0; it < NTRACK; ++it) {
		bool good = (track[it].hwPt < TKPT_MAX);
		for (int icalo = 0; icalo < NCALO; ++icalo) {
			if (calo_track_link_bit[it][icalo]) good = true;
		}
		if (good) {
			pfout[it].hwPt  = track[it].hwPt;
			pfout[it].hwEta = track[it].hwEta;
			pfout[it].hwPhi = track[it].hwPhi;
			pfout[it].hwId  = PID_Charged;
			pfout[it].hwZ0  = track[it].hwZ0;
		} else {
			pfout[it].hwPt  = 0;
			pfout[it].hwEta = 0;
			pfout[it].hwPhi = 0;
			pfout[it].hwId  = 0;
			pfout[it].hwZ0  = 0;
		}
	}
}

void _spfph_caloalgo(CaloObj calo[NCALO], pt_t sumtk[NCALO], int sumtkerr2[NCALO], PFNeutralObj pfout[NPF]) {
	for (int icalo = 0; icalo < NCALO; ++icalo) {
		pt_t calopt;
		if (sumtk[icalo] == 0) {
			calopt = calo[icalo].hwPt;
		} else {
			pt_t ptdiff = calo[icalo].hwPt - sumtk[icalo];
			if (ptdiff > 0 && (ptdiff*ptdiff) > sumtkerr2[icalo]) {
				calopt = ptdiff;
			} else {
				calopt = 0;
			}
		}
		pfout[icalo].hwPt  = calopt;
		pfout[icalo].hwEta = calopt ? calo[icalo].hwEta : etaphi_t(0);
		pfout[icalo].hwPhi = calopt ? calo[icalo].hwPhi : etaphi_t(0);
		pfout[icalo].hwId  = calopt ? PID_Neutral : 0;
	}
}

void simple_pflow_parallel_hwopt(CaloObj calo[NCALO], TkObj track[NTRACK], PFChargedObj outch[NTRACK], PFNeutralObj outne[NCALO]) {
	#pragma HLS ARRAY_PARTITION variable=calo complete
	#pragma HLS ARRAY_PARTITION variable=track complete
	#pragma HLS ARRAY_PARTITION variable=outch complete
	#pragma HLS ARRAY_PARTITION variable=outne complete

	#pragma HLS pipeline II=5

	bool calo_track_link_bit[NTRACK][NCALO];
	#pragma HLS ARRAY_PARTITION variable=calo_track_link_bit dim=0 complete

	_spfph_tk2calo_link(calo, track, calo_track_link_bit);

	int  tkerr2[NTRACK];
	#pragma HLS ARRAY_PARTITION variable=tkerr2 complete
	_spfph_tkerr2(track, tkerr2);

	pt_t sumtk[NCALO]; int sumtkerr2[NCALO];
	#pragma HLS ARRAY_PARTITION variable=sumtk complete
        #pragma HLS ARRAY_PARTITION variable=sumtkerr2 complete

	_spfph_tkalgo(track, calo_track_link_bit, outch);
	_spfph_sumtk(track, tkerr2, calo_track_link_bit, sumtk, sumtkerr2);
	_spfph_caloalgo(calo, sumtk, sumtkerr2, outne);
}

void medium_pflow_parallel_hwopt(CaloObj calo[NCALO], TkObj track[NTRACK], PFChargedObj outch[NTRACK], PFNeutralObj outne[NSELCALO]) {
	#pragma HLS ARRAY_PARTITION variable=calo complete
	#pragma HLS ARRAY_PARTITION variable=track complete
	#pragma HLS ARRAY_PARTITION variable=outch complete
	#pragma HLS ARRAY_PARTITION variable=outne complete

	#pragma HLS pipeline II=5
	PFNeutralObj outne_all[NCALO];
	#pragma HLS ARRAY_PARTITION variable=outne_all complete
	simple_pflow_parallel_hwopt(calo, track, outch, outne_all);
	ptsort_hwopt<PFNeutralObj,NCALO,NSELCALO>(outne_all, outne);
}

void simple_chs_hwopt(PFChargedObj pfch[NTRACK], z0_t pvZ, z0_t pvZCut, bool isPV[NTRACK]) {
	#pragma HLS ARRAY_PARTITION variable=pfch complete
	#pragma HLS ARRAY_PARTITION variable=isPV complete
	#pragma HLS pipeline II=5
        for (int it = 0; it < NTRACK; ++it) {
		if (pfch[it].hwPt == 0) {
			isPV[it] = false;
		} else {
			z0_t dz = pfch[it].hwZ0 - pvZ;
			isPV[it] = (-pvZCut <= dz && dz <= pvZCut);
		}
	}
}

void _lut_invert_init(ap_uint<16> _table[512]) {
	_table[0] = 32768;
	for (int i = 1; i <= 511; ++i) {
		_table[i] = (32768 / i);
	}
}
int _lut_divide(ap_uint<17> num, ap_uint<9> den) {
	ap_uint<16> _table[512];
	_lut_invert_init(_table);
	return (num * _table[den]);
}
void simple_puppi_hwopt(PFChargedObj pfch[NTRACK], bool isPV[NTRACK], PFNeutralObj pfne[NCALO], pt_t puppiPt[NCALO]) {
	#pragma HLS ARRAY_PARTITION variable=pfch complete
	#pragma HLS ARRAY_PARTITION variable=isPV complete
	#pragma HLS ARRAY_PARTITION variable=pfne complete
	#pragma HLS ARRAY_PARTITION variable=puppiPt complete
	#pragma HLS pipeline II=5

	const int DR2MAX = 8404; // 0.4 cone
	ap_uint<17> pt2[NTRACK];
	for (int it = 0; it < NTRACK; ++it) {
		int mypt2 = (pfch[it].hwPt*pfch[it].hwPt) >> 5;
		pt2[it] = (mypt2 < 131071 ? mypt2 : 131071);
	}
	for (int ic = 0; ic < NCALO; ++ic) {
		int sum = 0; pt_t ret = 0;
		for (int it = 0; it < NTRACK; ++it) {
			int dr2 = dr2_int(pfch[it].hwEta, pfch[it].hwPhi, pfne[ic].hwEta, pfne[ic].hwPhi);
			if (isPV[it] && dr2 <= DR2MAX) {
				ap_uint<9> dr2short = dr2 >> 5;
				//if (dr2short == 0) dr2short = 1;
				//sum += half(pfne[ic].hwPt*pfne[ic].hwPt)/half(dr2);
				//sum += ((pfch[it].hwPt*pfch[it].hwPt))/(dr2short);
				sum += _lut_divide(pt2[it], dr2short);
			}
		}
		if (sum != 0) {
			//sum = hls::log(sum);
			if (sum > 120) ret = pfne[ic].hwPt;
		}
		puppiPt[ic] = ret;
	}
}

void apply_chs_hwopt(PFChargedObj allch[NTRACK], bool isPV[NTRACK], PFChargedObj pvch[NPVTRACK]) {
#pragma HLS ARRAY_PARTITION variable=allch complete
#pragma HLS ARRAY_PARTITION variable=pvch complete
#pragma HLS ARRAY_PARTITION variable=isPV complete
#pragma HLS pipeline II=5

	PFChargedObj pvch_tmp[NPVTRACK];
	#pragma HLS ARRAY_PARTITION variable=pvch_tmp complete

	for (int iout = 0; iout < NPVTRACK; ++iout) {
		#pragma HLS unroll
		pvch_tmp[iout].hwPt = 0;
		pvch_tmp[iout].hwEta = 0;
		pvch_tmp[iout].hwPhi = 0;
		pvch_tmp[iout].hwId  = 0;
		pvch_tmp[iout].hwZ0  = 0;
	}

	int nout = 0;
	for (int it = 0; it < NTRACK; ++it) {
		if (isPV[it]) nout++;
		for (int iout = NPVTRACK-1; iout >= 0; --iout) {
			if (isPV[it] && nout <= NPVTRACK) {
				if (iout == 0) {
					pvch_tmp[iout] = allch[it];
				} else {
					pvch_tmp[iout] = pvch_tmp[iout-1];
				}
			//} else {
			//	pvch_tmp[iout] = pvch_tmp[iout];
			}
		}

	}
	for (int iout = 0; iout < NPVTRACK; ++iout) {
		pvch[iout] = pvch_tmp[iout];
	}
}

void apply_chs_sort_hwopt(PFChargedObj allch[NTRACK], bool isPV[NTRACK], PFChargedObj pvch[NPVTRACK]) {
#pragma HLS ARRAY_PARTITION variable=allch complete
#pragma HLS ARRAY_PARTITION variable=pvch complete
#pragma HLS ARRAY_PARTITION variable=isPV complete
#pragma HLS pipeline II=5

	PFChargedObj pvch_tmp[NPVTRACK];
	#pragma HLS ARRAY_PARTITION variable=pvch_tmp complete

	for (int iout = 0; iout < NPVTRACK; ++iout) {
		#pragma HLS unroll
		pvch_tmp[iout].hwPt = 0;
		pvch_tmp[iout].hwEta = 0;
		pvch_tmp[iout].hwPhi = 0;
		pvch_tmp[iout].hwId  = 0;
		pvch_tmp[iout].hwZ0  = 0;
	}

	for (int it = 0; it < NTRACK; ++it) {
		for (int iout = NPVTRACK-1; iout >= 0; --iout) {
			if (isPV[it] && pvch_tmp[iout].hwPt <= allch[it].hwPt) {
				if (iout == 0 || pvch_tmp[iout-1].hwPt > allch[it].hwPt) {
					pvch_tmp[iout] = allch[it];
				} else {
					pvch_tmp[iout] = pvch_tmp[iout-1];
				}
			}
		}

	}
	for (int iout = 0; iout < NPVTRACK; ++iout) {
		pvch[iout] = pvch_tmp[iout];
	}
}

void ptsort_pfneutral_hwopt(PFNeutralObj in[NCALO], PFNeutralObj out[NSELCALO]) {
#pragma HLS ARRAY_PARTITION variable=in complete
#pragma HLS ARRAY_PARTITION variable=out complete
#pragma HLS pipeline II=5 rewind
	ptsort_hwopt<PFNeutralObj,NCALO,NSELCALO>(in, out);
}
