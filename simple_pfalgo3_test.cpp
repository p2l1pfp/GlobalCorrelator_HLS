#include <cstdio>
#include "src/simple_pfalgo3.h"
#include "random_inputs.h"
#include "DiscretePFInputs_IO.h"
#include "pattern_serializer.h"

#define NTEST 500

bool em_equals(const EmCaloObj &out_ref, const EmCaloObj &out, const char *what, int idx) {
    bool ret;
    if (out_ref.hwPt == 0) {
        ret = (out.hwPt == 0);
    } else {
        ret = (out_ref.hwPt == out.hwPt && out_ref.hwPtErr == out.hwPtErr && out_ref.hwEta == out.hwEta && out_ref.hwPhi == out.hwPhi);
    }
    if  (!ret) {
        printf("Mismatch at %s[%3d], hwPt % 7d % 7d   hwPtErr % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d\n", what, idx,
                int(out_ref.hwPt), int(out.hwPt), int(out_ref.hwPtErr), int(out.hwPtErr), int(out_ref.hwEta), int(out.hwEta), int(out_ref.hwPhi), int(out.hwPhi));
    }
    return ret;
}
bool trk_equals(const TkObj &out_ref, const TkObj &out, const char *what, int idx) {
    bool ret;
    if (out_ref.hwPt == 0) {
        ret = (out.hwPt == 0);
    } else {
        ret = (out_ref.hwPt == out.hwPt && out_ref.hwPtErr == out.hwPtErr && out_ref.hwEta == out.hwEta && out_ref.hwPhi == out.hwPhi && out_ref.hwZ0  == out.hwZ0);
    }
    if  (!ret) {
        printf("Mismatch at %s[%3d], hwPt % 7d % 7d   hwPtErr % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwZ0 % 7d % 7d\n", what, idx,
                int(out_ref.hwPt), int(out.hwPt), int(out_ref.hwPtErr), int(out.hwPtErr), int(out_ref.hwEta), int(out.hwEta), int(out_ref.hwPhi), int(out.hwPhi), int(out_ref.hwZ0), int(out.hwZ0));
    }
    return ret;
}

bool had_equals(const HadCaloObj &out_ref, const HadCaloObj &out, const char *what, int idx) {
    bool ret;
    if (out_ref.hwPt == 0) {
        ret = (out.hwPt == 0);
    } else {
        ret = (out_ref.hwPt == out.hwPt && out_ref.hwEmPt == out.hwEmPt && out_ref.hwEta == out.hwEta && out_ref.hwPhi == out.hwPhi && out_ref.hwIsEM  == out.hwIsEM);
    }
    if  (!ret) {
        printf("Mismatch at %s[%3d], hwPt % 7d % 7d   hwEmPt % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwId %1d %1d\n", what, idx,
                int(out_ref.hwPt), int(out.hwPt), int(out_ref.hwEmPt), int(out.hwEmPt), int(out_ref.hwEta), int(out.hwEta), int(out_ref.hwPhi), int(out.hwPhi), int(out_ref.hwIsEM), int(out.hwIsEM));
    }
    return ret;
}
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
    
    HadCaloObj calo[NCALO]; EmCaloObj emcalo[NEMCALO]; TkObj track[NTRACK]; z0_t hwZPV;
    HadCaloObj calo_subem[NCALO], calo_subem_ref[NCALO]; 
    PFChargedObj outch[NTRACK], outch_ref[NTRACK];
    PFNeutralObj outpho[NPHOTON], outpho_ref[NPHOTON];
    PFNeutralObj outne[NSELCALO], outne_ref[NSELCALO];
#if defined(TESTMP7)
    MP7PatternSerializer serInPatterns("mp7_input_patterns.txt"), serOutPatterns("mp7_output_patterns.txt");
#endif
    HumanReadablePatternSerializer serHR("human_readable_patterns.txt");
    HumanReadablePatternSerializer debugHR("-"); // this will print on stdout, we'll use it for errors

    for (int test = 1; test <= NTEST; ++test) {
        for (int i = 0; i < NTRACK; ++i) {
            track[i].hwPt = 0; track[i].hwPtErr = 0; track[i].hwEta = 0; track[i].hwPhi = 0; track[i].hwZ0 = 0;
        }
        for (int i = 0; i < NCALO; ++i) {
            calo[i].hwPt = 0; calo[i].hwEmPt = 0; calo[i].hwEta = 0; calo[i].hwPhi = 0; calo[i].hwIsEM = 0; 
        }
        for (int i = 0; i < NEMCALO; ++i) {
            emcalo[i].hwPt = 0; emcalo[i].hwPtErr = 0;  emcalo[i].hwEta = 0; emcalo[i].hwPhi = 0;
        }


        if (!inputs.nextRegion(calo, emcalo, track, hwZPV)) break;


#if defined(TESTMP7) // Full PF, with MP7 wrapping 
        MP7DataWord data_in[MP7_NCHANN], data_out[MP7_NCHANN];
        mp7wrapped_pack_in(emcalo, calo, track, data_in);

#ifndef TESTMP7FAST // fast fake PF
        mp7wrapped_pfalgo3_full(data_in, data_out);
        pfalgo3_full_ref(emcalo, calo, track, outch_ref, outpho_ref, outne_ref);
#else  // standard PF
        mp7wrapped_pfalgo3_fast(data_in, data_out);
        pfalgo3_fast_ref(emcalo, calo, track, outch_ref, outpho_ref, outne_ref);
#endif
        mp7wrapped_unpack_out(data_out, outch, outpho, outne);
        // write out patterns for MP7 board hardware or simulator test
        serInPatterns(data_in); serOutPatterns(data_out);

#else // standard PFAlgo test without MP7 packing
        pfalgo3_full_ref(emcalo, calo, track, outch_ref, outpho_ref, outne_ref);
        pfalgo3_full(emcalo, calo, track, outch, outpho, outne);
#endif

        // write out human-readable patterns
        serHR(emcalo, calo, track, outch, outpho, outne);

        int errors = 0; int ntot = 0, npho = 0, nch = 0, nneu = 0;

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

        if (errors != 0) {
            printf("Error in computing test %d (%d)\n", test, errors);
            printf("Inputs: \n"); debugHR.dump_inputs(emcalo, calo, track);
            printf("Reference output: \n"); debugHR.dump_outputs(outch_ref, outpho_ref, outne_ref);
            printf("Current output: \n"); debugHR.dump_outputs(outch, outpho, outne);
            return 1;
        } else {
            printf("Passed test %d (%d, %d, %d, %d)\n", test, ntot, nch, npho, nneu);
        }

    }
    return 0;
}
