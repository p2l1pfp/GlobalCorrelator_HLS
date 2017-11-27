#include <cstdio>
#include "firmware/simple_fullpfalgo.h"
#include "firmware/simple_puppi.h"
#include "random_inputs.h"
#include "DiscretePFInputs_IO.h"
#include "pattern_serializer.h"
#include "test_utils.h"

#define NTEST 500


int main() {

    // input format: could be random or coming from simulation
    //RandomPFInputs inputs(37); // 37 is a good random number
    DiscretePFInputs inputs("regions_TTbar_PU140.dump");
    
    // input TP objects
    HadCaloObj calo[NCALO]; EmCaloObj emcalo[NEMCALO]; TkObj track[NTRACK]; z0_t hwZPV;
    HadCaloObj calo_subem[NCALO], calo_subem_ref[NCALO]; 
    MuObj mu[NMU];

    // output PF objects
    PFChargedObj outch[NTRACK], outch_ref[NTRACK];
    PFNeutralObj outpho[NPHOTON], outpho_ref[NPHOTON];
    PFNeutralObj outne[NSELCALO], outne_ref[NSELCALO];
    PFChargedObj outmupf[NMU], outmupf_ref[NMU];

    HumanReadablePatternSerializer serHR("human_readable_patterns.txt");
    HumanReadablePatternSerializer debugHR("-"); // this will print on stdout, we'll use it for errors

    // -----------------------------------------
    // run multiple tests
    for (int test = 1; test <= NTEST; ++test) {

        // initialize TP objects
        for (int i = 0; i < NTRACK; ++i) {
            track[i].hwPt = 0; track[i].hwPtErr = 0; track[i].hwEta = 0; track[i].hwPhi = 0; track[i].hwZ0 = 0; 
        }
        for (int i = 0; i < NCALO; ++i) {
            calo[i].hwPt = 0; calo[i].hwEmPt = 0; calo[i].hwEta = 0; calo[i].hwPhi = 0; calo[i].hwIsEM = 0; 
        }
        for (int i = 0; i < NEMCALO; ++i) {
            emcalo[i].hwPt = 0; emcalo[i].hwPtErr = 0;  emcalo[i].hwEta = 0; emcalo[i].hwPhi = 0;
        }
        for (int i = 0; i < NMU; ++i) {
            mu[i].hwPt = 0; mu[i].hwPtErr = 0; mu[i].hwEta = 0; mu[i].hwPhi = 0;
        }

        // get the inputs from the input object
        if (!inputs.nextRegion(calo, emcalo, track, mu, hwZPV)) break;

        pfalgo3_full_ref(emcalo, calo, track, mu, outch_ref, outpho_ref, outne_ref, outmupf_ref);
        // pfalgo3_full(emcalo, calo, track, mu, outch, outpho, outne, outmupf);

        // simple_puppi_ref(outch_ref, outpho_ref, outne_ref, pz, out)
        simple_puppi_hw(outch_ref);

        // write out human-readable patterns
        // serHR(emcalo, calo, track, mu, outch, outpho, outne, outmupf);


        // // -----------------------------------------
        // // validation against the reference algorithm
        // int errors = 0; int ntot = 0, npho = 0, nch = 0, nneu = 0, nmu = 0;

        // // check charged hadrons
        // for (int i = 0; i < NTRACK; ++i) {
        //     if (!pf_equals(outch_ref[i], outch[i], "PF Charged", i)) errors++;
        //     if (outch_ref[i].hwPt > 0) { ntot++; nch++; }
        // }
        // // check photon 
        // for (int i = 0; i < NPHOTON; ++i) {
        //     if (!pf_equals(outpho_ref[i], outpho[i], "Photon", i)) errors++;
        //     if (outpho_ref[i].hwPt > 0) { ntot++; npho++; }
        // }
        // for (int i = 0; i < NSELCALO; ++i) {
        //     if (!pf_equals(outne_ref[i], outne[i], "PF Neutral", i)) errors++;
        //     if (outne_ref[i].hwPt > 0) { ntot++; nneu++; }
        // }
        // for (int i = 0; i < NMU; ++i) {
        //     if (!pf_equals(outmupf_ref[i], outmupf[i], "PF Muon", i)) errors++;
        //     if (outmupf_ref[i].hwPt > 0) { ntot++; nmu++; }
        // }        

        // if (errors != 0) {
        //     printf("Error in computing test %d (%d)\n", test, errors);
        //     printf("Inputs: \n"); debugHR.dump_inputs(emcalo, calo, track, mu);
        //     printf("Reference output: \n"); debugHR.dump_outputs(outch_ref, outpho_ref, outne_ref, outmupf_ref);
        //     printf("Current output: \n"); debugHR.dump_outputs(outch, outpho, outne, outmupf);
        //     return 1;
        // } else {
        //     printf("Passed test %d (%d, %d, %d, %d)\n", test, ntot, nch, npho, nneu);
        // }

    }
    return 0;
}
