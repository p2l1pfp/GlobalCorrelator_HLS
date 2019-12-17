#include <cstdio>
#include "firmware/pfalgo3.h"
#include "ref/pfalgo3_ref.h"
#include "utils/DiscretePFInputsReader.h"
#include "utils/pattern_serializer.h"
#include "utils/test_utils.h"

#define NTEST 1000

int main() {
    HumanReadablePatternSerializer debugHR("-", /*zerosuppress=*/true); // this will print on stdout, we'll use it for errors

    DiscretePFInputsReader inputs("TTbar_PU200_Barrel.dump");
    
    // input TP objects
    HadCaloObj hadcalo[NCALO]; EmCaloObj emcalo[NEMCALO]; TkObj track[NTRACK]; z0_t hwZPV;
    MuObj mu[NMU];

    // output PF objects
    PFChargedObj outch[NTRACK], outch_ref[NTRACK];
    PFNeutralObj outne[NSELCALO], outne_ref[NSELCALO];
    PFNeutralObj outpho[NPHOTON], outpho_ref[NPHOTON];
    PFChargedObj outmupf[NMU], outmupf_ref[NMU];

    pfalgo3_config cfg(NTRACK,NEMCALO,NCALO,NMU, 
                       NPHOTON,NSELCALO,NALLNEUTRALS,
                       PFALGO_DR2MAX_TK_MU, PFALGO_DR2MAX_TK_EM, PFALGO_DR2MAX_EM_CALO, PFALGO_DR2MAX_TK_CALO,
                       PFALGO_TK_MAXINVPT_LOOSE, PFALGO_TK_MAXINVPT_TIGHT);

#if defined(PACKING_DATA_SIZE) && defined(PACKING_NCHANN)
    PatternSerializer serPatternsIn("pfalgo3_input_patterns.txt"), serPatternsOut("pfalgo3_output_patterns.txt");
    ap_uint<PACKING_DATA_SIZE> packed_input[PACKING_NCHANN], packed_output[PACKING_NCHANN];
    for (unsigned int i = 0; i < PACKING_NCHANN; ++i) { packed_input[i] = 0; packed_output[i] = 0; }
#endif

    // -----------------------------------------
    // run multiple tests
    for (int test = 1; test <= NTEST; ++test) {
        // get the inputs from the input object
        if (!inputs.nextRegion(hadcalo, emcalo, track, mu, hwZPV)) break;

        bool verbose = false; // can set this on to get detailed printout of some test
        pfalgo3_set_debug(verbose);
        pfalgo3_ref_set_debug(verbose);

#if defined(PACKING_DATA_SIZE) && defined(PACKING_NCHANN)
        pfalgo3_pack_in(emcalo, hadcalo, track, mu, packed_input); 
        serPatternsIn(packed_input);
        packed_pfalgo3(packed_input, packed_output);
        serPatternsOut(packed_output);
        pfalgo3_unpack_out(packed_output, outch, outpho, outne, outmupf);
#else
        pfalgo3(emcalo, hadcalo, track, mu, outch, outpho, outne, outmupf);
#endif

        pfalgo3_ref(cfg, emcalo, hadcalo, track, mu, outch_ref, outpho_ref, outne_ref, outmupf_ref);

        // -----------------------------------------
        // validation against the reference algorithm
        int errors = 0; int ntot = 0, nch = 0, npho = 0, nneu = 0, nmu = 0;

        for (int i = 0; i < NTRACK; ++i) {
            if (!pf_equals(outch_ref[i], outch[i], "PF Charged", i)) errors++;
            if (outch_ref[i].hwPt > 0) { ntot++; nch++; }
        }
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
            printf("Inputs: \n"); debugHR.dump_inputs(emcalo, hadcalo, track, mu);
            printf("Reference output: \n"); debugHR.dump_outputs(outch_ref, outpho_ref, outne_ref, outmupf_ref);
            printf("Current output: \n"); debugHR.dump_outputs(outch, outpho, outne, outmupf);
            return 1;
        } else {
            printf("Passed test %d (%d, %d, %d, %d)\n", test, ntot, nch, npho, nneu);
        }

    }
    return 0;
}
