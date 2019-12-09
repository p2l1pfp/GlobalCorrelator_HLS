#include <cstdio>
#include "firmware/simple_fullpfalgo.h"
#include "utils/DiscretePFInputsReader.h"
#include "utils/pattern_serializer.h"
#include "utils/test_utils.h"

#define NTEST 1000


int main() {

    DiscretePFInputsReader inputs("TTbar_PU200_Barrel.dump");
    
    // input TP objects
    HadCaloObj calo[NCALO]; EmCaloObj emcalo[NEMCALO]; TkObj track[NTRACK]; z0_t hwZPV;
    MuObj mu[NMU];

    // output PF objects
    PFChargedObj outch[NTRACK], outch_ref[NTRACK];
    PFNeutralObj outpho[NPHOTON], outpho_ref[NPHOTON];
    PFNeutralObj outne[NSELCALO], outne_ref[NSELCALO];
    PFChargedObj outmupf[NMU], outmupf_ref[NMU];
    HumanReadablePatternSerializer serHR("human_readable_patterns.txt");
    HumanReadablePatternSerializer debugHR("-", /*zerosuppress=*/true); // this will print on stdout, we'll use it for errors;


    // -----------------------------------------
    // run multiple tests
    for (int test = 1; test <= NTEST; ++test) {

        // get the inputs from the input object
        if (!inputs.nextRegion(calo, emcalo, track, mu, hwZPV)) break;

        pfalgo3_full_ref_set_debug(test == 209);
        pfalgo3_full_ref(emcalo, calo, track, mu, outch_ref, outpho_ref, outne_ref, outmupf_ref);
        pfalgo3_full_set_debug(test == 209);
        pfalgo3_full(emcalo, calo, track, mu, outch, outpho, outne, outmupf);

        // write out human-readable patterns
        serHR(emcalo, calo, track, mu, outch, outpho, outne, outmupf);

        // -----------------------------------------
        // validation against the reference algorithm
        int errors = 0; int ntot = 0, npho = 0, nch = 0, nneu = 0, nmu = 0;

        // check charged hadrons
        for (int i = 0; i < NTRACK; ++i) {
            if (!pf_equals(outch_ref[i], outch[i], "PF Charged", i)) errors++;
            if (outch_ref[i].hwPt > 0) { ntot++; nch++; }
        }
        // check photon 
        for (int i = 0; i < NPHOTON; ++i) {
            if (!pf_equals(outpho_ref[i], outpho[i], "Photon", i)) errors++;
            if (outpho_ref[i].hwPt > 0) { ntot++; npho++; }
        }
        for (int i = 0; i < NSELCALO; ++i) {
            if (!pf_equals(outne_ref[i], outne[i], "PF Neutral", i)) errors++;
            if (outne_ref[i].hwPt > 0) { ntot++; nneu++; }
        }
        for (int i = 0; i < NMU; ++i) {
            if (!pf_equals(outmupf_ref[i], outmupf[i], "PF Muon", i)) errors++;
            if (outmupf_ref[i].hwPt > 0) { ntot++; nmu++; }
        }        

        if (errors != 0) {
            printf("Error in computing test %d (%d)\n", test, errors);
            printf("Inputs: \n"); debugHR.dump_inputs(emcalo, calo, track, mu);
            printf("Reference output: \n"); debugHR.dump_outputs(outch_ref, outpho_ref, outne_ref, outmupf_ref);
            printf("Current output: \n"); debugHR.dump_outputs(outch, outpho, outne, outmupf);
            return 1;
        } else {
            printf("Passed test %d (%d, %d, %d, %d)\n", test, ntot, nch, npho, nneu);
        }

    }
    return 0;
}
