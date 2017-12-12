#include "simple_puppi.h"
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

typedef ap_uint<7> tk2em_dr_t;
typedef ap_uint<10> tk2calo_dr_t;
typedef ap_uint<10> em2calo_dr_t;
typedef ap_uint<12> tk2calo_dq_t;
typedef ap_uint<12> mu2trk_dr_t;

int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}

void simple_puppi_hw(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne[NNEUTRALS], z0_t Z0) {


    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete

    #pragma HLS pipeline II=1

    z0_t DZMAX = 256; // placeholder, big window (z0 is 10-bits)

    // input to the table is a 12-bit number, 2^12 = 4096
    ap_uint<8> puppiweight_table[4096];
    lut_puppiweight_init< ap_uint<8>, 4096 >( puppiweight_table );

    // compute alpha
    const int DR2MAX = 8404; // 0.4 cone
    AlphaCalc: for (int in = 0; in < NNEUTRALS; ++in) {
        
        if (pfallne[in].hwPt == 0) continue;
        
        int sum = 0;
        for (int it = 0; it < NTRACK; ++it) {
            if ((Z0 - pfch[it].hwZ0 > DZMAX) || (Z0 - pfch[it].hwZ0 < -DZMAX)) continue; // if track is PV
            int dr2 = dr2_int(pfch[it].hwEta, pfch[it].hwPhi, pfallne[in].hwEta, pfallne[in].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
                ap_uint<9> dr2short = dr2 >> 5; // why?
                int pt2 = pfch[it].hwPt*pfch[it].hwPt;
                int term = (std::min(pt2 >> 5, 131071) << 15)/(dr2short > 0 ? dr2short : ap_uint<9>(1));
                sum += term;
            }
        }
      
		int weight = 0;
        if (sum > 0) { // get the weight if e^alpha is not 0

            // also need to check if alpha > alphaMED
            if (sum < 22026) weight = 0; // < e^10 where that is the median
            else if (sum > 1202604) weight = 256; // e^14 where that is the median
            else{
                // bitshift down the sum by 10 bits to get it in the right range
                int sum_short = sum >> 10;
                // -- compute weight --
                // sum is the current e^alpha
                // med of alpha = 10, rms of alpha = 2
                // pass in sum, alphamed, alpharms into a LUT
                int index = sum_short >= 4096 ? 4096 : sum_short; // 2^16
                // ap_uint<8> weight_short = puppiweight_table[index];            
                weight = (int) puppiweight_table[index];
                // std::cout << "sum = " << sum << ", weight = " << weight << std::endl;
            }
        }

        int pT_new = ( pfallne[in].hwPt * weight ) >> 8;
        pfallne[in].hwPtPuppi = (pt_t) pT_new;
    }

}


