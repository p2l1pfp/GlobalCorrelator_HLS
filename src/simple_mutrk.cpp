#include "simple_mutrk.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

// bool match_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, etaphi_t boxSize) {
// 	etaphi_t deta = (eta1-eta2);
// 	etaphi_t dphi = (phi1-phi2);
// 	return (deta <= boxSize && deta >= -boxSize && dphi <= boxSize && dphi >= -boxSize);
// }
// etaphi_t dr_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
// 	etaphi_t deta = eta1 > eta2 ? (eta1-eta2) : (eta2-eta1);
// 	etaphi_t dphi = phi1 > phi2 ? (phi1-phi2) : (phi2-phi1);
// 	return (deta > dphi ? deta : dphi);
// }
// int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
// 	etaphi_t deta = (eta1-eta2);
// 	etaphi_t dphi = (phi1-phi2);
// 	return deta*deta + dphi*dphi;
// }

template<unsigned NB>
ap_uint<NB> dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) {
	etaphi_t deta = (eta1-eta2);
	etaphi_t dphi = (phi1-phi2);
	int dr2 = deta*deta + dphi*dphi;
	return (dr2 < int(max) ? ap_uint<NB>(dr2) : max);
}

#define mu2trk_dr_t ap_uint<4>
void spfph_mu2trk_drvals(MuObj mu[NCALO], TkObj track[NTRACK], mu2trk_dr_t mu_track_drval[NMU][NTRACK]) {

	const mu2trk_dr_t DR2MAX = 2101;
	for (int im = 0; im < NMU; ++im) {
		pt_t muPtMin = mu[im].hwPt - 2*(mu[im].hwPtErr);
		for (int it = 0; it < NTRACK; ++it) {
			if (track[it].hwPt > muPtMin) {
				mu_track_drval[im][it] = dr2_int_cap<4>(mu[im].hwEta, mu[im].hwPhi, track[it].hwEta, track[it].hwPhi, DR2MAX);
			} else {
				mu_track_drval[im][it] = DR2MAX;
			}
		}
	}
}

void spfph_mu2trk_linkstep(mu2trk_dr_t mu_track_drval[NMU][NTRACK], ap_uint<NMU> mu_track_link_bit[NTRACK]) {
	
	const mu2trk_dr_t DR2MAX = 2101;
	for (int im = 0; im < NMU; ++im) {
		for (int it = 0; it < NTRACK; ++it) {
			mu2trk_dr_t mydr = mu_track_drval[im][it];
			bool link = (mydr != DR2MAX);
			for (int j = 0; j < NTRACK; ++j) {
				if (it <= j) link = link && (mu_track_drval[im][j] >= mydr);
				else         link = link && (mu_track_drval[im][j] >  mydr);
			}
			mu_track_link_bit[it][im] = link;
		}
	}
}

void spfph_mutrk_link(MuObj mu[NMU], TkObj track[NTRACK], ap_uint<NMU> mu_track_link_bit[NTRACK]) {
	
	#pragma HLS ARRAY_PARTITION variable=mu complete
	#pragma HLS ARRAY_PARTITION variable=track complete
	#pragma HLS ARRAY_PARTITION variable=mu_track_link_bit complete dim=0

	#pragma HLS pipeline II=5

	mu2trk_dr_t drvals[NMU][NTRACK];
	#pragma HLS ARRAY_PARTITION variable=drvals complete dim=0

	spfph_mu2trk_drvals(mu, track, drvals);
	spfph_mu2trk_linkstep(drvals, mu_track_link_bit);
}

void spfph_mualgo(MuObj mu[NMU], TkObj track[NTRACK], ap_uint<NMU> mu_track_link_bit[NTRACK], PFMuonObj pfmuout[NMU]) {
	
	const pt_t TKPT_MAX = 80; // 20 * PT_SCALE;
	for (int im = 0; im < NMU; ++im) {
		bool good = false;
		for (int it = 0; it < NTRACK; ++it) {
			if (mu_track_link_bit[it][im]) good = true;
			track[it].hwIsMu  = true;
		}
		if (good) {
			pfmuout[im].hwPt  = track[im].hwPt;
			pfmuout[im].hwEta = track[im].hwEta;
			pfmuout[im].hwPhi = track[im].hwPhi;
			pfmuout[im].hwId  = PID_Muon;
			pfmuout[im].hwZ0  = track[im].hwZ0;
		} else {
			pfmuout[im].hwPt  = 0;
			pfmuout[im].hwEta = 0;
			pfmuout[im].hwPhi = 0;
			pfmuout[im].hwId  = 0;
			pfmuout[im].hwZ0  = 0;
		}
	}
}

void simple_mutrk_parallel_hwopt(MuObj mu[NMU], TkObj track[NTRACK], PFMuonObj outmu[NMU]) {
	
	#pragma HLS ARRAY_PARTITION variable=mu complete
	#pragma HLS ARRAY_PARTITION variable=track complete
	#pragma HLS ARRAY_PARTITION variable=outmu complete

	#pragma HLS pipeline II=5

	ap_uint<NMU> mu_track_link_bit[NTRACK];
	#pragma HLS ARRAY_PARTITION variable=mu_track_link_bit complete

	spfph_mutrk_link(mu, track, mu_track_link_bit);

	// int  tkerr2[NTRACK];
	// #pragma HLS ARRAY_PARTITION variable=tkerr2 complete
	// _spfph_tkerr2(track, tkerr2);

	// pt_t sumtk[NCALO]; int sumtkerr2[NCALO];
	// #pragma HLS ARRAY_PARTITION variable=sumtk complete
 //    #pragma HLS ARRAY_PARTITION variable=sumtkerr2 complete

	spfph_mualgo(mu, track, mu_track_link_bit, outmu);
	// _spfph_sumtk(track, tkerr2, calo_track_link_bit, sumtk, sumtkerr2);
	// _spfph_caloalgo(calo, sumtk, sumtkerr2, outne);
}

void ptsort_pfneutral_hwopt(PFNeutralObj in[NCALO], PFNeutralObj out[NSELCALO]) {
#pragma HLS ARRAY_PARTITION variable=in complete
#pragma HLS ARRAY_PARTITION variable=out complete
#pragma HLS pipeline II=5 rewind
	ptsort_hwopt<PFNeutralObj,NCALO,NSELCALO>(in, out);
}
