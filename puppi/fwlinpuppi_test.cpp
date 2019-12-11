#include <cstdio>
#include "firmware/linpuppi.h"
#include "../utils/DiscretePFInputsReader.h"
#include "../utils/pattern_serializer.h"
#include "../utils/test_utils.h"

#define NTEST 100


int main() {

    //DiscretePFInputsReader inputs("TTbar_PU200_HGCalNoTK.dump");
    DiscretePFInputsReader inputs("VBFHToBB_PU200_HGCalNoTK.dump");
    
    // input TP objects (used)
    HadCaloObj calo[NCALO];

    // input TP objects (unused, but needed to read the inputs)
    EmCaloObj emcalo[NEMCALO]; TkObj track[NTRACK]; MuObj mu[NMU]; z0_t hwZPV;

    // input/output PFPUPPI objects
    PFNeutralObj outselne[NNEUTRALS];
    PFNeutralObj outselne_ref[NNEUTRALS];
    PFNeutralObj outselne_flt[NNEUTRALS];
    PFNeutralObj outallne[NCALO];
    PFNeutralObj outallne_ref[NCALO];
    PFNeutralObj outallne_flt[NCALO];


    int flt_nok = 0, flt_1bit = 0, flt_almostok = 0, flt_nz = 0, flt_nmiss = 0;
    float flt_sumDiff, flt_sumAbsDiff = 0;
    for (int test = 1; test <= NTEST; ++test) {
        // get the inputs from the input object
        if (!inputs.nextRegion(calo, emcalo, track, mu, hwZPV)) break;

#ifdef TEST_PT_CUT
        float minpt = 0;
        for (unsigned int i = 0; i < NCALO; ++i) minpt += calo[i].hwPt*LINPUPPI_ptLSB;
        if (minpt < TEST_PT_CUT) { 
            //std::cout << "Skipping region with total calo pt " << minpt << " below threshold." << std::endl; 
            --test; continue; 
        }
#endif

        bool verbose = true;
        if (verbose) printf("test case %d\n", test);
        fwdlinpuppi_set_debug(verbose);

#if defined(TEST_PUPPI_NOCROP)
        fwdlinpuppiNoCrop(calo, outallne);
#else
        fwdlinpuppi(calo, outselne);
#endif
        fwdlinpuppi_ref(calo, outallne_ref, outselne_ref, verbose);
        fwdlinpuppi_flt(calo, outallne_flt, outselne_flt, verbose);

        for (int i = 0; i < NCALO; ++i){
            if (calo[i].hwPt > 0) {
                if (outallne_flt[i].hwPtPuppi > 2/LINPUPPI_ptLSB) flt_nz++;
                if (outallne_flt[i].hwPtPuppi == outallne_ref[i].hwPtPuppi) {
                    flt_nok++; 
                } else {
                    if (std::abs<int>(outallne_flt[i].hwPtPuppi - outallne_ref[i].hwPtPuppi)*LINPUPPI_ptLSB < 1 + 0.01 * outallne_flt[i].hwPtPuppi) {
                        if (std::abs<int>(outallne_flt[i].hwPtPuppi - outallne_ref[i].hwPtPuppi) == 1) {
                            flt_1bit++; 
                        } else {
                            flt_almostok++;
                        }
                    }
                    flt_nmiss++;
                }
                float ptDiff = (outallne_flt[i].hwPtPuppi - outallne_ref[i].hwPtPuppi) * LINPUPPI_ptLSB;
                flt_sumDiff += ptDiff;
                flt_sumAbsDiff += std::abs(ptDiff);
                if (verbose)  printf("particle %02d pT %7.2f  em %1d :  puppiPt_ref %7.2f   puppiPt_flt %7.2f    diff %+7.2f\n",
                                        i, calo[i].hwPt * LINPUPPI_ptLSB, int(calo[i].hwIsEM), 
                                        outallne_ref[i].hwPtPuppi * LINPUPPI_ptLSB, outallne_flt[i].hwPtPuppi * LINPUPPI_ptLSB, ptDiff);
            }
#if defined(TEST_PUPPI_NOCROP)
            // apply pT cut to the reference
            if (outallne_ref[i].hwPtPuppi < LINPUPPI_ptCut/LINPUPPI_ptLSB) clear(outallne_ref[i]);
            if (!puppi_equals(outallne_ref[i], outallne[i], "Puppi", i)) {
                printf("Mismatch in test %d: particle %02d pT %7.2f  em %1d :  puppiPt_hw %7.2f    puppiPt_ref %7.2f   puppiPt_flt %7.2f\n", test,
                    i, calo[i].hwPt * LINPUPPI_ptLSB, int(calo[i].hwIsEM), outallne[i].hwPtPuppi * LINPUPPI_ptLSB, outallne_ref[i].hwPtPuppi * LINPUPPI_ptLSB, outallne_flt[i].hwPtPuppi * LINPUPPI_ptLSB);
                return 1;
            }
#else
            if (i <= NSELCALO && !puppi_equals(outselne_ref[i], outselne[i], "Puppi", i)) {
                printf("Mismatch in test %d: particle %02d pT %7.2f  em %1d :  puppiPt_hw %7.2f    puppiPt_ref %7.2f   puppiPt_flt %7.2f\n", test,
                    i, calo[i].hwPt * LINPUPPI_ptLSB, int(calo[i].hwIsEM), outselne[i].hwPtPuppi * LINPUPPI_ptLSB, outselne_ref[i].hwPtPuppi * LINPUPPI_ptLSB, outselne_flt[i].hwPtPuppi * LINPUPPI_ptLSB);
                return 1;
            }
#endif
        }

        if (verbose) printf("\n");
        else         printf("passed test %d\n", test);

    }

    printf("Integer vs floating point accuracy (%d non-empty candidates out of %d regions cropped at NCALO=%d, %d candidates after PUPPI pT > 2):\n", flt_nok+flt_nmiss, NTEST, NCALO, flt_nz);
    printf("  - exact match: %6d  (%6.2f%% ):\n", flt_nok,   flt_nok * 100.0 / (flt_nok+flt_nmiss));
    printf("  - mismatch   : %6d  (%6.2f%% ):\n", flt_nmiss, flt_nmiss * 100.0 / (flt_nok+flt_nmiss));
    printf("  -   by 1*LSB : %6d  (%6.2f%% )   [ 1 unit, %.2f GeV ]\n", flt_1bit, flt_1bit * 100.0 / (flt_nok+flt_nmiss), LINPUPPI_ptLSB);
    printf("  -      small : %6d  (%6.2f%% )   [ %.2f < delta(pt) <= 1 GeV + 1% ]\n", flt_almostok, flt_almostok * 100.0 / (flt_nok+flt_nmiss), LINPUPPI_ptLSB);
    printf("  -      big   : %6d  (%6.2f%% )   [ delta(pt) > 1 GeV + 1% ]\n", (flt_nmiss-flt_almostok-flt_1bit), (flt_nmiss-flt_almostok-flt_1bit) * 100.0 / (flt_nok+flt_nmiss));
    printf("  - average pT  diff   %+8.4f  (on all)    %+8.4f  (on mismatch)\n", flt_sumDiff/(flt_nok+flt_nmiss), flt_sumDiff/std::max(flt_nmiss,1));
    printf("  - average pT |diff|  % 8.4f  (on all)    % 8.4f  (on mismatch)\n", flt_sumAbsDiff/(flt_nok+flt_nmiss), flt_sumAbsDiff/std::max(flt_nmiss,1));
    return 0;
}
