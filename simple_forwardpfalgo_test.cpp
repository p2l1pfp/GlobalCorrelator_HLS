#include <cstdio>
#include "firmware/simple_forwardpfalgo.h"
#include "puppi/firmware/simple_puppi_forward.h"
#include "utils/random_inputs.h"
#include "utils/DiscretePFInputs_IO.h"
#include "utils/pattern_serializer.h"
#include "utils/test_utils.h"

#define NTEST 10


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
#if defined(TESTMP7)
    MP7PatternSerializer serInPatterns( "mp7_input_patterns.txt", HLS_pipeline_II,HLS_pipeline_II-1); // mux each event into HLS_pipeline_II frames
    MP7PatternSerializer serOutPatterns("mp7_output_patterns.txt",HLS_pipeline_II,HLS_pipeline_II-1); // assume only one PF core running per chip,
    MP7PatternSerializer serInPatterns2( "mp7_input_patterns_magic.txt", HLS_pipeline_II,-HLS_pipeline_II+1); // mux each event into HLS_pipeline_II frames
    MP7PatternSerializer serOutPatterns2("mp7_output_patterns_magic.txt",HLS_pipeline_II,-HLS_pipeline_II+1); // assume only one PF core running per chip,
    MP7PatternSerializer serInPatterns3( "mp7_input_patterns_nomux.txt");  // 
    MP7PatternSerializer serOutPatterns3("mp7_output_patterns_nomux.txt"); // ,
#endif
#if defined(TESTCTP7)
    CTP7PatternSerializer serInPatterns4( "ctp7_input_patterns_nomux.txt",CTP7_NCHANN_IN, true);  // 
    CTP7PatternSerializer serOutPatterns4("ctp7_output_patterns_nomux.txt",CTP7_NCHANN_OUT, false); // fill the rest of the lines with empty events for now
#endif
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


#if defined(TESTMP7) // Full PF, with MP7 wrapping 
        MP7DataWord data_in[MP7_NCHANN], data_out[MP7_NCHANN];
        // initialize
        for (int i = 0; i < MP7_NCHANN; ++i) {
            data_in[i] = 0;
            data_out[i] = 0;
        }
        mp7wrapped_pack_in(emcalo, calo, mu, data_in);
        MP7_TOP_FUNC(data_in, data_out);
        //mp7wrapped_unpack_out(data_out, outpho, outne, outmupf);
        mp7wrapped_unpack_out_necomb(data_out, outpho, outne, outmupf);
		// for (int ii = 0; ii < 72; ++ii){ std::cout << ii << ", " << data_in[ii] << std::endl; }
		
        MP7_REF_FUNC(emcalo, calo, mu, outpho_ref, outne_ref, outmupf_ref);

        // write out patterns for MP7 board hardware or simulator test
        serInPatterns(data_in); serOutPatterns(data_out);
        serInPatterns2(data_in); serOutPatterns2(data_out);
        serInPatterns3(data_in); serOutPatterns3(data_out);

#elif defined(TESTCTP7) // Full PF, with CTP7 wrapping
        MP7DataWord data_in[CTP7_NCHANN_IN], data_out[CTP7_NCHANN_OUT];
        // initialize
        for (int i = 0; i < CTP7_NCHANN_IN; ++i) { data_in[i] = 0; }
        for (int i = 0; i < CTP7_NCHANN_OUT; ++i) { data_out[i] = 0; }
        mp7wrapped_pack_in(emcalo, calo, track, mu, data_in);
        MP7_TOP_FUNC(data_in, data_out);
        //mp7wrapped_unpack_out(data_out, outch, outpho, outne, outmupf);
        mp7wrapped_unpack_out_necomb(data_out, outch, outpho, outne, outmupf);
    
        MP7_REF_FUNC(emcalo, calo, track, mu, outch_ref, outpho_ref, outne_ref, outmupf_ref);
        // write out patterns for CTP7 board hardware or simulator test
        serInPatterns4(data_in,CTP7_NCHANN_IN); serOutPatterns4(data_out,CTP7_NCHANN_OUT);       

#else // standard PFAlgo test without MP7 packing
        pfalgo3_full_ref(emcalo, calo, track, mu, outch_ref, outpho_ref, outne_ref, outmupf_ref);
        pfalgo3_full(emcalo, calo, track, mu, outch, outpho, outne, outmupf);
#endif
        // write out human-readable patterns
        serHR(emcalo, calo, track, mu, outch, outpho, outne, outmupf);

        PFNeutralObj outallne_ref[NNEUTRALS];
        // sort/merge neutrals (at the end of the PF algo?) - do it in test bench for now
        for (int ipfne = 0; ipfne < NPHOTON; ++ipfne){
            outpho_ref[ipfne].hwPtPuppi = 0;
            outallne_ref[ipfne] = outpho_ref[ipfne];
        }
        for (int ipfne = NPHOTON; ipfne < NNEUTRALS; ++ipfne){
            outne_ref[ipfne-NPHOTON].hwPtPuppi = 0;
            outallne_ref[ipfne] = outne_ref[ipfne-NPHOTON];
        }

        simple_puppi_forward_ref( outallne_ref);

#ifdef TESTMP7
        if (!MP7_VALIDATE) continue;
#endif
#ifdef TESTCTP7
        if (!CTP7_VALIDATE) continue;
#endif

        // -----------------------------------------
        // validation against the reference algorithm
        int errors = 0; int ntot = 0, npho = 0, nch = 0, nneu = 0, nmu = 0;

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

        int puperrors = 0;
        for (int i = 0; i < NPHOTON; ++i){
            printf("hwpt = %i, hwptpuppi = %i, refptpuppi = %i, hw-ref_ptpuppi = %i \n", (int) outpho[i].hwPt, (int) outpho[i].hwPtPuppi, (int) outallne_ref[i].hwPtPuppi, int(outpho[i].hwPtPuppi-outallne_ref[i].hwPtPuppi));
            if (outpho[i].hwPtPuppi-outallne_ref[i].hwPtPuppi != 0 && outpho[i].hwPt>0) puperrors++;
        }
        for (int i = 0; i < NSELCALO; ++i){
            printf("hwpt = %i, hwptpuppi = %i, refptpuppi = %i, hw-ref_ptpuppi = %i \n", (int) outne[i].hwPt, (int) outne[i].hwPtPuppi, (int) outallne_ref[i+NPHOTON].hwPtPuppi, int(outne[i].hwPtPuppi-outallne_ref[i+NPHOTON].hwPtPuppi));
            if (outne[i].hwPtPuppi-outallne_ref[i+NPHOTON].hwPtPuppi != 0 && outne[i].hwPt>0) puperrors++;
        }
        std::cout << "end of test ---- " << test << std::endl;

        if (puperrors>0) {
            printf("Found %i errors in puppi test!\n", puperrors);
            //return errors;
        }

    }
    return 0;
}
