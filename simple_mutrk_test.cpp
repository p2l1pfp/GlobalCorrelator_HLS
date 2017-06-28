#include <cstdio>
#include "src/simple_mutrk.h"
#include "random_inputs.h"
#include "DiscretePFInputs_IO.h"

#define NTEST 10

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
bool pf_equals(const PFMuonObj &out_ref, const PFMuonObj &out, const char *what, int idx) {
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

// ----------------------------------------------------------------------------------------
// ---------------- MAIN ------------------------------------------------------------------

int main() {

    RandomPFInputs inputs(37); // 37 is a good random number
	//DiscretePFInputs inputs("regions_TTbar_PU140.dump");
	
	CaloObj calo[NCALO]; TkObj track[NTRACK]; z0_t hwZPV, hwZ0Cut = 7;
	MuObj mu[NMU];
    PFChargedObj outch[NTRACK], outch_ref[NTRACK];
    PFNeutralObj outne_unsorted_ref[NSELCALO], outne_unsorted[NSELCALO];
    PFNeutralObj outne_ref[NSELCALO], outne[NSELCALO];
    PFMuonObj outmupf[NMU], outmupf_ref[NMU];
    TkObj outtrk[NTRACK], outtrk_ref[NTRACK];
	
	for (int test = 1; test <= NTEST; ++test) {
		
		for (int i = 0; i < NTRACK; ++i) {
			track[i].hwPt = 0; track[i].hwPtErr = 0; track[i].hwEta = 0; track[i].hwPhi = 0; track[i].hwZ0 = 0; track[i].hwIsMu = false;
		}
		for (int i = 0; i < NCALO; ++i) {
			calo[i].hwPt = 0; calo[i].hwEta = 0; calo[i].hwPhi = 0;
		}
		for (int i = 0; i < NMU; ++i) {
			mu[i].hwPt = 0; calo[i].hwEta = 0; calo[i].hwPhi = 0;
		}

		if (!inputs.nextRegion(calo, track, mu, hwZPV)) break;

		// ---------------- TOP FUNCTIONS and REF ----------------
		simple_mutrk_parallel_ref(mu,track,outmupf_ref);
		simple_mutrk_parallel_hwopt(mu,track,outmupf);

		
		// ---------------- COMPARE WITH EXPECTED ----------------
		int errors = 0; int ntot = 0, nmu = 0, nneu = 0;

		// == MU-TRK COMPARISON ===
		std::cout << ">> Test number = " << test << std::endl;
		
		// for (int i = 0; i < NMU; ++i) {
		// 	printf("muon  %3d, hwPt % 7d   hwPtErr % 7d    hwEta %+7d   hwPhi %+7d\n", i, int(mu[i].hwPt), int(mu[i].hwPtErr), int(mu[i].hwEta), int(mu[i].hwPhi));
		// }
		// for (int i = 0; i < NTRACK; ++i) {
		// 	printf("track %3d, hwPt % 7d   hwPtErr % 7d    hwEta %+7d   hwPhi %+7d     hwZ0 %+7d\n", i, int(track[i].hwPt), int(track[i].hwPtErr), int(track[i].hwEta), int(track[i].hwPhi), int(track[i].hwZ0));
		// }		
		for (int i = 0; i < NMU; ++i) {
			if (!pf_equals(outmupf_ref[i], outmupf[i], "PF Muon", i)) errors++;
			if (outmupf_ref[i].hwPt > 0) { ntot++; nmu++; }
		}

		if (errors != 0) {
			printf("Error in computing test %d (%d)\n", test, errors);
			for (int i = 0; i < NMU; ++i) {
				printf("muon  %3d, hwPt % 7d   hwPtErr % 7d    hwEta %+7d   hwPhi %+7d\n", i, int(mu[i].hwPt), int(mu[i].hwPtErr), int(mu[i].hwEta), int(mu[i].hwPhi));
			}
			for (int i = 0; i < NTRACK; ++i) {
				printf("track %3d, hwPt % 7d   hwPtErr % 7d    hwEta %+7d   hwPhi %+7d     hwZ0 %+7d\n", i, int(track[i].hwPt), int(track[i].hwPtErr), int(track[i].hwEta), int(track[i].hwPhi), int(track[i].hwZ0));
			}

			return 1;
		} else {
			printf("<< Passed test %d (%d, %d)\n", test, ntot, nmu);
		}

	}
	return 0;
}
