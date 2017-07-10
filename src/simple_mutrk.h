#ifndef SIMPLE_PFLOW_SIMPLE_PFLOW_H
#define SIMPLE_PFLOW_SIMPLE_PFLOW_H

#include "data.h"

// bool match_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, etaphi_t boxSize) ;
// etaphi_t dr_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) ;
int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) ;
ap_uint<12> dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<12> max) ;

void simple_mutrk_parallel_ref(MuObj mu[NMU], TkObj track[NTRACK], PFChargedObj outmu[NMU]) ;
void simple_mutrk_parallel_hwopt(MuObj mu[NMU], TkObj track[NTRACK], PFChargedObj outmu[NMU]) ;

template<typename T, int NIn, int NOut>
void ptsort_ref(T in[NIn], T out[NOut]) {
	for (int iout = 0; iout < NOut; ++iout) {
		out[iout].hwPt = 0;
	}
	int nout = 0;
	for (int it = 0; it < NIn; ++it) {
		for (int iout = 0; iout < NOut; ++iout) {
			if (in[it].hwPt >= out[iout].hwPt) {
				for (int i2 = NOut-1; i2 > iout; --i2) {
					out[i2] = out[i2-1];
				}
				out[iout] = in[it];
				break;
			}
		}
	}
}

template<typename T, int NIn, int NOut>
void ptsort_hwopt(T in[NIn], T out[NOut]) {
#pragma HLS ARRAY_PARTITION variable=in complete
#pragma HLS ARRAY_PARTITION variable=out complete
#pragma HLS pipeline II=5 rewind

	T tmp[NOut];
	#pragma HLS ARRAY_PARTITION variable=tmp complete

	for (int iout = 0; iout < NOut; ++iout) {
		#pragma HLS unroll
		tmp[iout].hwPt = 0;
	}

	for (int it = 0; it < NIn; ++it) {
		for (int iout = NOut-1; iout >= 0; --iout) {
			if (tmp[iout].hwPt <= in[it].hwPt) {
				if (iout == 0 || tmp[iout-1].hwPt > in[it].hwPt) {
					tmp[iout] = in[it];
				} else {
					tmp[iout] = tmp[iout-1];
				}
			}
		}

	}
	for (int iout = 0; iout < NOut; ++iout) {
		out[iout] = tmp[iout];
	}

}
// void ptsort_pfneutral_hwopt(PFNeutralObj inne[NCALO], PFNeutralObj outne[NSELCALO]);

#endif
