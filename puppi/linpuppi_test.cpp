#include <cstdio>
#include "firmware/linpuppi.h"
#include "linpuppi_ref.h"
#include "../utils/DiscretePFInputsReader.h"
#include "../utils/pattern_serializer.h"
#include "../utils/test_utils.h"
#include "puppi_checker.h"

#if defined(REG_Barrel)
    #include "../pfalgo3_ref.h"
#elif defined(REG_HGCal)
    #include "../pfalgo2hgc_ref.h"
#endif

#define NTEST 100


int main() {
#if defined(REG_Barrel)
    DiscretePFInputsReader inputs("TTbar_PU200_Barrel.dump");
    pfalgo3_config pfcfg(NTRACK,NEMCALO,NCALO,NMU, 
                         NPHOTON,NSELCALO,NALLNEUTRALS,
                         PFALGO_DR2MAX_TK_MU, PFALGO_DR2MAX_TK_EM, PFALGO_DR2MAX_EM_CALO, PFALGO_DR2MAX_TK_CALO,
                         PFALGO_TK_MAXINVPT_LOOSE, PFALGO_TK_MAXINVPT_TIGHT);
    linpuppi_config pucfg(NTRACK, NALLNEUTRALS, NNEUTRALS,
                          LINPUPPI_DR2MIN, LINPUPPI_DR2MAX, LINPUPPI_ptMax, LINPUPPI_dzCut,
                          LINPUPPI_ptSlopeNe, LINPUPPI_ptSlopePh, LINPUPPI_ptZeroNe, LINPUPPI_ptZeroPh, 
                          LINPUPPI_alphaSlope, LINPUPPI_alphaZero, LINPUPPI_alphaCrop, 
                          LINPUPPI_priorNe, LINPUPPI_priorPh,
                          LINPUPPI_ptCut);
#elif defined(REG_HGCal)
    DiscretePFInputsReader inputs("TTbar_PU200_HGCal.dump");
    pfalgo_config pfcfg(NTRACK,NCALO,NMU, NSELCALO,
                        PFALGO_DR2MAX_TK_MU, PFALGO_DR2MAX_TK_CALO,
                        PFALGO_TK_MAXINVPT_LOOSE, PFALGO_TK_MAXINVPT_TIGHT);
#endif
    
    // input TP objects and PV
    HadCaloObj hadcalo[NCALO]; EmCaloObj emcalo[NEMCALO]; TkObj track[NTRACK]; MuObj mu[NMU]; 
    z0_t hwZPV;

    // PF objects
    PFChargedObj pfch[NTRACK], pfmu[NMU];
    PFNeutralObj pfpho[NPHOTON], pfne[NSELCALO], pfallne[NALLNEUTRALS];

    // Puppi objects
    PFChargedObj outallch[NTRACK], outallch_ref[NTRACK];
    PFNeutralObj outallne[NALLNEUTRALS], outallne_ref_nocut[NALLNEUTRALS], outallne_ref[NALLNEUTRALS], outallne_flt_nocut[NALLNEUTRALS], outallne_flt[NALLNEUTRALS];
    PFNeutralObj outselne[NNEUTRALS], outselne_ref[NNEUTRALS], outselne_flt[NNEUTRALS];

    PuppiChecker checker;

    for (int test = 1; test <= NTEST; ++test) {
        // get the inputs from the input object
        if (!inputs.nextRegion(hadcalo, emcalo, track, mu, hwZPV)) break;

#ifdef TEST_PT_CUT
        float minpt = 0;
        for (unsigned int i = 0; i < NTRACK; ++i) minpt += track[i].hwPt*LINPUPPI_ptLSB;
        if (minpt < TEST_PT_CUT) { 
            //std::cout << "Skipping region with total calo pt " << minpt << " below threshold." << std::endl; 
            --test; continue; 
        }
#endif

#if defined(REG_Barrel)
        pfalgo3_ref(pfcfg, emcalo, hadcalo, track, mu, pfch, pfpho, pfne, pfmu);
        pfalgo3_merge_neutrals_ref(pfcfg, pfpho, pfne, pfallne);
#elif defined(REG_HGCal)
        pfalgo2hgc(pfcfg, hadcalo, track, mu, pfch, pfallne, pfmu); 
#endif

        bool verbose = true;
        if (verbose) printf("test case %d\n", test);
        linpuppi_set_debug(verbose);

        linpuppi_chs(hwZPV, pfch, outallch);
#if defined(TEST_PUPPI_NOCROP)
        linpuppiNoCrop(track, hwZPV, pfallne, outallne);
#else
        linpuppi(track, hwZPV, pfallne, outselne);
#endif

        linpuppi_chs_ref(pucfg, hwZPV, pfch, outallch_ref);
        linpuppi_ref(pucfg, track, hwZPV, pfallne, outallne_ref_nocut, outallne_ref, outselne_ref, verbose);
        linpuppi_flt(pucfg, track, hwZPV, pfallne, outallne_flt_nocut, outallne_flt, outselne_flt, verbose);

        // validate numerical accuracy 
        checker.checkIntVsFloat<PFNeutralObj,NALLNEUTRALS>(pfallne, outallne_ref_nocut, outallne_flt_nocut, verbose);

        bool ok = checker.checkChs<NTRACK>(hwZPV, outallch, outallch_ref) && 
#if defined(TEST_PUPPI_NOCROP)
                  checker.check<NALLNEUTRALS>(outallne, outallne_ref, outallne_flt);
#else
                  checker.check<NNEUTRALS>(outselne, outselne_ref, outselne_flt);
#endif
        if (!ok) {
            printf("FAILED test %d\n", test);
            HumanReadablePatternSerializer dumper("-", true);
#if defined(TEST_PUPPI_NOCROP)
            dumper.dump_puppi(NALLNEUTRALS, "all    ", outallne);
            dumper.dump_puppi(NALLNEUTRALS, "all ref", outallne_ref);
#else
            dumper.dump_puppi(NNEUTRALS,    "sel    ", outselne);
            dumper.dump_puppi(NNEUTRALS,    "sel ref", outselne_ref);
#endif
            dumper.dump_puppi(NALLNEUTRALS, "all rnc", outallne_ref_nocut);
            dumper.dump_puppi(NALLNEUTRALS, "all flt", outallne_flt_nocut);
            return 1;
        }

        if (verbose) printf("\n");
        else         printf("passed test %d\n", test);

    }

    printf("Report for %d regions (cropped at N=%d):\n", NTEST, NALLNEUTRALS);
    checker.printIntVsFloatReport();
    return 0;
}
