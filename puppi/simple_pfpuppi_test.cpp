#include <cstdio>
#include "../firmware/simple_fullpfalgo.h"
#include "../vertexing/firmware/simple_vtx.h"
#include "firmware/simple_puppi.h"
#include "../utils/random_inputs.h"
#include "../utils/DiscretePFInputs_IO.h"
#include "../utils/pattern_serializer.h"
#include "../utils/test_utils.h"

#define NTEST 50


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
    ap_int<16>   trk_z0;

    // input/output PFPUPPI objects
    PFNeutralObj outallne[NNEUTRALS];
    PFNeutralObj outallne_ref[NNEUTRALS];
    // PFChargedObj pupch_pv[NTRACK],    pupch_pv_ref[NTRACK];
    // PFNeutralObj pupne_pv[NNEUTRALS], pupne_pv_ref[NNEUTRALS];

    // HumanReadablePatternSerializer serHR("human_readable_patterns.txt");
    // HumanReadablePatternSerializer debugHR("-"); // this will print on stdout, we'll use it for errors

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

        // sort/merge neutrals (at the end of the PF algo?) - do it in test bench for now
        for (int ipfne = 0; ipfne < NPHOTON; ++ipfne){
            outpho_ref[ipfne].hwPtPuppi = 0;
            outallne_ref[ipfne] = outpho_ref[ipfne];
            outallne[ipfne] = outpho_ref[ipfne];
        }
        for (int ipfne = NPHOTON; ipfne < NNEUTRALS; ++ipfne){
            outne_ref[ipfne-NPHOTON].hwPtPuppi = 0;
            outallne_ref[ipfne] = outne_ref[ipfne-NPHOTON];
            outallne[ipfne] = outne_ref[ipfne-NPHOTON];
        }

        std::cout << "test " << test << std::endl;

        ap_uint<8> weights[NNEUTRALS];
        VtxObj curvtx;    
        simple_vtx_ref(track,&curvtx);
        simple_puppi_ref( outch_ref, outallne_ref, curvtx.hwZ0);
        simple_puppi_hw(  outch_ref, outallne,     curvtx.hwZ0);
        // weight_t curweight;
        // compute_puppi_weight_hw( 100, curweight );
        // std::cout << "curweight = " << curweight << std::endl;

        for (int i = 0; i < NNEUTRALS; ++i){
            printf("hwpt = %i, hwptpuppi = %i, hwptpuppi-ref = %i \n", (int) outallne[i].hwPt, (int) outallne[i].hwPtPuppi, (int) outallne_ref[i].hwPtPuppi);
        }
        std::cout << "end of test ---- " << test << std::endl;

    }


    return 0;
}
