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
    linpuppi_config pucfg(NTRACK, NCALO, NNEUTRALS,
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
        for (unsigned int i = 0; i < NCALO; ++i) minpt += hadcalo[i].hwPt*LINPUPPI_ptLSB;
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

        linpuppi_chs_ref(pucfg, hwZPV, pfch, outallch_ref);
        linpuppi_ref(pucfg, track, hwZPV, pfallne, outallne_ref_nocut, outallne_ref, outselne_ref, verbose);
        linpuppi_flt(pucfg, track, hwZPV, pfallne, outallne_flt_nocut, outallne_flt, outselne_flt, verbose);

        // validate numerical accuracy 
        checker.checkIntVsFloat<PFNeutralObj,NALLNEUTRALS>(pfallne, outallne_ref_nocut, outallne_flt_nocut, verbose);

        bool ok = true;
                  //checker.check<NTRACK>(outallch, outallch_ref, outallch_flt) && 
                  //checker.check<NALLNEUTRALS>(outallne, outallne_ref, outallne_flt);
        if (!ok) {
            printf("FAILED test %d\n", test);
            HumanReadablePatternSerializer dumper("-", true);
            dumper.dump_puppi(NALLNEUTRALS, "    ", outallne);
            dumper.dump_puppi(NALLNEUTRALS, "ref ", outallne_ref);
            dumper.dump_puppi(NALLNEUTRALS, "rnc ", outallne_ref_nocut);
            dumper.dump_puppi(NALLNEUTRALS, "flt ", outallne_flt);
            return 1;
        }

        if (verbose) printf("\n");
        else         printf("passed test %d\n", test);

    }

    printf("Report for %d regions (cropped at N=%d):\n", NTEST, NALLNEUTRALS);
    checker.printIntVsFloatReport();
    return 0;
}
