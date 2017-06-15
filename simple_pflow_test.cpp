#include <cstdio>
#include "src/simple_pflow.h"
#include "random_inputs.h"
#include "DiscretePFInputs_IO.h"

#define NTEST 500

bool pf_equals(const PFChargedObj &out_ref, const PFChargedObj &out, const char *what, int idx) {
	bool ret;
	if (out_ref.hwPt == 0) {
		ret = (out.hwPt == 0);
	} else {
		ret = (out_ref.hwPt == out.hwPt && out_ref.hwEta == out.hwEta && out_ref.hwPhi == out.hwPhi && out_ref.hwId  == out.hwId && out_ref.hwZ0  == out.hwZ0);
	}
	if  (!ret) {
		printf("Mismatch at %s[%3d], hwPt % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwId %1d %1d      hwZ0 %+7d %+7d   \n", what, idx,
				int(out_ref.hwPt), int(out.hwPt),
				int(out_ref.hwEta), int(out.hwEta),
				int(out_ref.hwPhi), int(out.hwPhi),
				int(out_ref.hwId), int(out.hwId),
				int(out_ref.hwZ0), int(out.hwZ0));
	}
	return ret;
}
bool pf_equals(const PFNeutralObj &out_ref, const PFNeutralObj &out, const char *what, int idx) {
	bool ret;
	if (out_ref.hwPt == 0) {
		ret = (out.hwPt == 0);
	} else {
		ret = (out_ref.hwPt == out.hwPt && out_ref.hwEta == out.hwEta && out_ref.hwPhi == out.hwPhi && out_ref.hwId  == out.hwId);
	}
	if  (!ret) {
		printf("Mismatch at %s[%3d], hwPt % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwId %1d %1d \n", what, idx,
				int(out_ref.hwPt), int(out.hwPt),
				int(out_ref.hwEta), int(out.hwEta),
				int(out_ref.hwPhi), int(out.hwPhi),
				int(out_ref.hwId), int(out.hwId));
	}
	return ret;
}
int main() {

    //RandomPFInputs inputs(37); // 37 is a good random number
	DiscretePFInputs inputs("regions_TTbar_PU140.dump");
	
	CaloObj calo[NCALO]; TkObj track[NTRACK]; z0_t hwZPV, hwZ0Cut = 7;
    PFChargedObj outch[NTRACK], outch_ref[NTRACK];
    PFNeutralObj outne[NCALO], outne_ref[NCALO];

	for (int test = 1; test <= NTEST; ++test) {
		for (int i = 0; i < NTRACK; ++i) {
			track[i].hwPt = 0; track[i].hwPtErr = 0; track[i].hwEta = 0; track[i].hwPhi = 0; track[i].hwZ0 = 0;
		}
		for (int i = 0; i < NCALO; ++i) {
			calo[i].hwPt = 0; calo[i].hwEta = 0; calo[i].hwPhi = 0;
		}

		if (!inputs.nextRegion(calo, track, hwZPV)) break;

		//bool calo_track_link_bit[NTRACK][NCALO];
		//link_pflow_parallel_hwopt(calo, track, calo_track_link_bit);
		ap_uint<NCALO> calo_track_link_bit[NTRACK];
		//link_pflow_parallel_hwopt2(calo, track, calo_track_link_bit);
		spfph_tk2calo_link_v2(calo, track, calo_track_link_bit);
		//spfph_tk2calo_link_v3(calo, track, calo_track_link_bit);
		ap_uint<NCALO> calo_track_link_bit_ref[NTRACK];
		link_pflow_parallel_ref(calo, track, calo_track_link_bit_ref);
		//simple_pflow_iterative_ref(calo, track, out_ref);
		//simple_pflow_iterative_hwopt(calo, track, out);
		//simple_pflow_parallel_ref(calo, track, outch_ref, outne_ref);
		simple_pflow_parallel_hwopt(calo, track, outch, outne);
		//medium_pflow_parallel_ref(calo, track, outch_ref, outne_ref);
		//medium_pflow_parallel_hwopt(calo, track, outch, outne);

		//PFNeutralObj outne_sorted_ref[NSELCALO], outne_sorted[NSELCALO];
		//ptsort_ref<PFNeutralObj,NCALO,NSELCALO>(outne_ref, outne_sorted_ref);
		//ptsort_pfneutral_hwopt(outne, outne_sorted);

// ---------------- COMPARE WITH EXPECTED ----------------
/*
		int errors = 0; int ntot = 0, nch = 0, nneu = 0;
		for (int i = 0; i < NTRACK; ++i) { for (int j = 0; j < NCALO; ++j) {
			if (calo_track_link_bit[i][j] != calo_track_link_bit_ref[i][j]) {
				printf("mismatch track-calo link[%3d][%3d] = %1d (hwopt), %1d (ref)\n",
						i,j,int(calo_track_link_bit[i][j]),int(calo_track_link_bit_ref[i][j]));
				errors++;
			}
		} }
		if (errors) {
			printf("          ");
			for (int j = 0; j < NCALO; ++j) printf("C%02d   ", j);
			printf("\n");
			for (int i = 0; i < NTRACK; ++i) {
				printf("TRACK %2d: ",i);
				for (int j = 0; j < NCALO; ++j) printf("%1d %1d   ", int(calo_track_link_bit[i][j]),int(calo_track_link_bit_ref[i][j]));
				printf("\n");
			}
		}
		for (int i = 0; i < NTRACK; ++i) {
			if (!pf_equals(outch_ref[i], outch[i], "PF Charged", i)) errors++;
			if (outch_ref[i].hwPt > 0) { ntot++; nch++; }
		}
		for (int i = 0; i < NCALO; ++i) {
			if (!pf_equals(outne_ref[i], outne[i], "PF Neutral", i)) errors++;
			if (outne_ref[i].hwPt > 0) { ntot++; nneu++; }
		}
		//for (int i = 0; i < NSELCALO; ++i) {
		//	if (!pf_equals(outne_sorted_ref[i], outne_sorted[i], "PF sorted neutral", i)) errors++;
		//}

		// ------- run CHS and PUPPI ------
		bool isPV[NTRACK], isPV_ref[NTRACK];
		pt_t puppiPt[NCALO], puppiPt_ref[NCALO];

		simple_chs_ref(outch_ref, hwZPV, hwZ0Cut, isPV_ref) ;
		simple_chs_hwopt(outch, hwZPV, hwZ0Cut, isPV) ;
		simple_puppi_ref(outch_ref, isPV_ref, outne_ref, puppiPt_ref) ;
		simple_puppi_hwopt(outch, isPV, outne, puppiPt) ;

		PFChargedObj pvch[NTRACK], pvch_ref[NTRACK];
		apply_chs_sort_ref(outch_ref, isPV_ref, pvch_ref);
		apply_chs_sort_hwopt(outch, isPV, pvch);


// ---------------- COMPARE WITH EXPECTED ----------------

		int nchpv = 0, nnepv = 0;
		for (int i = 0; i < NTRACK; ++i) {
			if (outch_ref[i].hwPt > 0 && isPV[i] != isPV_ref[i]) errors++;
			if (isPV_ref[i]) nchpv++;
		}
		for (int i = 0; i < NCALO; ++i) {
			if (outne_ref[i].hwPt > 0 && puppiPt[i] != puppiPt_ref[i]) errors++;
			if (puppiPt[i]) nnepv++;
		}
		for (int i = 0; i < NPVTRACK; ++i) {
			if (!pf_equals(pvch_ref[i], pvch[i], "PF CHS", i)) errors++;
		}



		if (errors != 0) {
			printf("Error in computing test %d (%d)\n", test, errors);
			for (int i = 0; i < NCALO; ++i) {
				printf("calo  %3d, hwPt % 7d   hwPtErr % 7d    hwEta %+7d   hwPhi %+7d\n", i, int(calo[i].hwPt), 0, int(calo[i].hwEta), int(calo[i].hwPhi));
			}
			for (int i = 0; i < NTRACK; ++i) {
				printf("track %3d, hwPt % 7d   hwPtErr % 7d    hwEta %+7d   hwPhi %+7d     hwZ0 %+7d\n", i, int(track[i].hwPt), int(track[i].hwPtErr), int(track[i].hwEta), int(track[i].hwPhi), int(track[i].hwZ0));
			}
			for (int i = 0; i < NTRACK; ++i) {
				printf("charged pf %3d, hwPt % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwId %1d %1d      hwZ0 %+7d %+7d    isPV %1d %1d\n", i,
					int(outch_ref[i].hwPt), int(outch[i].hwPt), int(outch_ref[i].hwEta), int(outch[i].hwEta),
					int(outch_ref[i].hwPhi), int(outch[i].hwPhi), int(outch_ref[i].hwId), int(outch[i].hwId),
					int(outch_ref[i].hwZ0), int(outch[i].hwZ0), int(isPV_ref[i]), int(isPV[i]));
			}
			for (int i = 0; i < NCALO; ++i) {
				printf("neutral pf %3d, hwPt % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwId %1d %1d      puppiPt % 7d % 7d \n", i,
					int(outne_ref[i].hwPt), int(outne[i].hwPt), int(outne_ref[i].hwEta), int(outne[i].hwEta),
					int(outne_ref[i].hwPhi), int(outne[i].hwPhi), int(outne_ref[i].hwId), int(outne[i].hwId),
					int(puppiPt_ref[i]), int(puppiPt[i]));
			}
			for (int i = 0; i < NPVTRACK; ++i) {
				printf("CHS charged pf %3d, hwPt % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwId %1d %1d      hwZ0 %+7d %+7d\n", i,
					int(pvch_ref[i].hwPt), int(pvch[i].hwPt), int(pvch_ref[i].hwEta), int(pvch[i].hwEta),
					int(pvch_ref[i].hwPhi), int(pvch[i].hwPhi), int(pvch_ref[i].hwId), int(pvch[i].hwId),
					int(pvch_ref[i].hwZ0), int(pvch[i].hwZ0));
			}
			return 1;
		} else {
			printf("Passed test %d (%d, %d, %d)\n", test, ntot, nch, nneu);
		}
	*/

	}
	return 0;
}
