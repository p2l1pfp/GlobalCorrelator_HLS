#include <cstdio>
#include "simple_pflow.h"

#define NTEST 50

int main() {

	srand(37); // 37 is a good random number
	
	CaloObj calo[NCALO]; TkObj track[NTRACK];
    PFObj out[NPF], out_ref[NPF];

	for (int test = 1; test <= NTEST; ++test) {
		for (int i = 0; i < NTRACK; ++i) {
			track[i].hwPt = 0; track[i].hwPtErr = 0; track[i].hwEta = 0; track[i].hwPhi = 0;
		}
		for (int i = 0; i < NCALO; ++i) {
			calo[i].hwPt = 0; calo[i].hwEta = 0; calo[i].hwPhi = 0;
		}
		int ncharged = (rand() % NTRACK/2) + NTRACK/2;
		int nneutral = (rand() % ((3*NCALO)/4));
		for (int i = 1; i < nneutral && i < NCALO; i += 2) {
			float pt = (rand()/float(RAND_MAX))*80+1, eta = (rand()/float(RAND_MAX))*2.0-1.0, phi = (rand()/float(RAND_MAX))*2.0-1.0;
			calo[i].hwPt  = pt * PT_SCALE;
			calo[i].hwEta = eta * ETAPHI_SCALE;
			calo[i].hwPhi = phi * ETAPHI_SCALE;
		}
		for (int i = 0; i < ncharged && i < NTRACK; ++i) {
			float pt = (rand()/float(RAND_MAX))*50+2, eta = (rand()/float(RAND_MAX))*2.0-1.0, phi = (rand()/float(RAND_MAX))*2.0-1.0;
			track[i].hwPt    = pt * PT_SCALE;
			track[i].hwPtErr = (0.2*pt+4) * PT_SCALE; 
			track[i].hwEta = eta * ETAPHI_SCALE;
			track[i].hwPhi = phi * ETAPHI_SCALE;
			int icalo = rand() % NCALO;
			if (i % 3 == 1 || icalo >= NCALO) continue;
			float dpt_calo = ((rand()/float(RAND_MAX))*3-1.5) * (0.2*pt+4);
			float deta_calo = ((rand()/float(RAND_MAX))*0.3-0.15), dphi_calo = ((rand()/float(RAND_MAX))*0.3-0.15);
			if (pt + dpt_calo > 0) {
				calo[icalo].hwPt  += (pt + dpt_calo) * PT_SCALE;
				calo[icalo].hwEta = (eta + deta_calo) * ETAPHI_SCALE;
				calo[icalo].hwPhi = (phi + dphi_calo) * ETAPHI_SCALE;
			}
		}

		//simple_pflow_iterative_ref(calo, track, out_ref);
		//simple_pflow_iterative_hwopt(calo, track, out);
		simple_pflow_parallel_ref(calo, track, out_ref);
		simple_pflow_parallel_hwopt(calo, track, out);

// ---------------- COMPARE WITH EXPECTED ----------------

		int errors = 0; int ntot = 0, nch = 0, nneu = 0;
		for (int i = 0; i < NPF; ++i) {
			if (out_ref[i].hwPt == 0) {
				if (out[i].hwPt != 0) {
					printf("Mismatch at %3d, hwPt % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwId %1d %1d \n", i,
							int(out_ref[i].hwPt), int(out[i].hwPt),
							int(out_ref[i].hwEta), int(out[i].hwEta),
							int(out_ref[i].hwId), int(out[i].hwPhi),
							int(out_ref[i].hwId), int(out[i].hwId));
					errors++;
				}
				continue;
			}
			ntot++; 
			nch  += out_ref[i].hwId == PID_Charged;
			nneu += out_ref[i].hwId == PID_Neutral;
			if (out_ref[i].hwPt  != out[i].hwPt || out_ref[i].hwEta != out[i].hwEta || out_ref[i].hwPhi != out[i].hwPhi || out_ref[i].hwId  != out[i].hwId) {
				printf("Mismatch at %3d, hwPt % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwId %1d %1d \n", i,
						int(out_ref[i].hwPt), int(out[i].hwPt),
						int(out_ref[i].hwEta), int(out[i].hwEta),
						int(out_ref[i].hwId), int(out[i].hwPhi),
						int(out_ref[i].hwId), int(out[i].hwId));
				errors++;
			}
		}
		if (errors != 0) {
			printf("Error in computing test %d (%d)\n", test, errors);
			for (int i = 0; i < NCALO; ++i) {
				printf("calo  %3d, hwPt % 7d   hwPtErr % 7d    hwEta %+7d   hwPhi %+7d\n", i, int(calo[i].hwPt), 0, int(calo[i].hwEta), int(calo[i].hwPhi));
			}
			for (int i = 0; i < NTRACK; ++i) {
				printf("track %3d, hwPt % 7d   hwPtErr % 7d    hwEta %+7d   hwPhi %+7d\n", i, int(track[i].hwPt), int(track[i].hwPtErr), int(track[i].hwEta), int(track[i].hwPhi));
			}
			for (int i = 0; i < NPF; ++i) {
				printf("pf %3d, hwPt % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwId %1d %1d \n", i,
					int(out_ref[i].hwPt), int(out[i].hwPt), int(out_ref[i].hwEta), int(out[i].hwEta),
					int(out_ref[i].hwId), int(out[i].hwPhi), int(out_ref[i].hwId), int(out[i].hwId));
			}
			return 1;
		} else {
			printf("Passed test %d (%d, %d, %d)\n", test, ntot, nch, nneu);
		}
	}	
	return 0;
}
