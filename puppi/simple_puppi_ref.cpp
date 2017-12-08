#include "../firmware/data.h"
#include "../firmware/simple_fullpfalgo.h"
#include "firmware/simple_puppi.h"
#include <cmath>
#include <algorithm>

void simple_puppi_ref(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne[NNEUTRALS], z0_t Z0, 
                      PFChargedObj out_pfch_pv[NTRACK], PFNeutralObj out_pfne_pv[NNEUTRALS]) {

    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS ARRAY_PARTITION variable=out_pfch_pv complete
    #pragma HLS ARRAY_PARTITION variable=out_pfne_pv complete
    
    #pragma HLS pipeline II=1


    z0_t DZMAX = 128; // placeholder

    // e^10 = 22000, e^16 = 8,000,000, e^14 = 1,200,000
    // input to the table is a 12-bit number, 2^12 = 4096
    // input to the table is a 12-bit number, 2^14 = 16384
    // input to the table is a 16-bit number, 2^16 = 65536
    // input to the table is a 32-bit number, 2^32 = 4,294,967,296
    ap_uint<8> puppiweight_table[4096];
    _lut_puppiweight_init< ap_uint<8>, 4096 >( puppiweight_table );

    // compute alpha
    const int DR2MAX = 8404; // 0.4 cone
    for (int in = 0; in < NNEUTRALS; ++in) {
        
        if (pfallne[in].hwPt == 0) continue;
        
        // computing the alpha for the analysis
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
      
        // compute the weight
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
                int index = sum_short >= 4096 ? 4096 : sum; // 2^16
                // ap_uint<8> weight_short = puppiweight_table[index];            
                // ap_uint<6> weight = puppiweight_table[1000];            
                // weight = (int) weight_short;
                // weight = 0;
                weight = (int) puppiweight_table[sum_short];
                // printf("sum = %i, sum_short = %i, weight = %i \n", sum, sum_short, (int) puppiweight_table[sum_short]);
            }
        }

        int pT_new = ( pfallne[in].hwPt * weight ) >> 8;
        pfallne[in].hwPtPuppi = (pt_t) pT_new;
        // printf("old pt = %i, new pt = %i \n", (int) pfallne[in].hwPt, (int) pT_new);
    }

}

