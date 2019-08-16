#include "../firmware/data.h"
#include "../firmware/simple_forwardpfalgo.h"
#include "firmware/simple_puppi_forward.h"
#include <cmath>
#include <algorithm>

void simple_puppi_forward_ref(PFNeutralObj pfallne[NNEUTRALS]) {

    // input to the table is a 12-bit number, 2^12 = 4096
    ap_uint<8> puppiweight_table[PUPPI_TABLE_SIZE];
    lut_puppiweight_init< ap_uint<8>, PUPPI_TABLE_SIZE >( puppiweight_table );

    // compute alpha
    const int DR2MAX = 8404; // 0.4 cone
    for (int in = 0; in < NNEUTRALS; ++in) {
        
        //if (pfallne[in].hwPt == 0) continue;
        int sum = 0;
        for (int it = 0; it < NNEUTRALS; ++it) {
            int dr2 = dr2_int(pfallne[it].hwEta, pfallne[it].hwPhi, pfallne[in].hwEta, pfallne[in].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX && it!=in) {
                //std::cout<<"SW:  "<<dr2<<"\t "<<in<<" "<<it<<std::endl;
                ap_uint<9> dr2short = dr2 >> 5; // why?
                int pt2 = pfallne[it].hwPt*pfallne[it].hwPt;
                //int term = (std::min(pt2 >> 5, 131071) << 15)/(dr2short > 0 ? dr2short : ap_uint<9>(1));
                int term = 32768/(dr2short > 0 ? dr2short : ap_uint<9>(1));
                term = std::min(pt2 >> 5, 131071)*term;
                sum += term;
            }
        }
      
        int weight = 0;
        int eToAlpha = sum >> 10;
        if (eToAlpha > 0) { // get the weight if e^alpha is not 0
            if (eToAlpha <= 0) weight = 0; // < e^10 where that is the median
            //else if (eToAlpha > 1174) weight = 256; // e^14 where that is the median
            else if (eToAlpha >= 1174) weight = 256; // e^14 where that is the median
            else{
                // bitshift down the sum by 10 bits to get it in the right range
                int sum_short = eToAlpha;
                // -- compute weight --
                // sum is the current e^alpha
                // med of alpha = 10, rms of alpha = 2
                // pass in sum, alphamed, alpharms into a LUT
                int index = sum_short >= PUPPI_TABLE_SIZE ? PUPPI_TABLE_SIZE : sum_short; // 2^16
                // ap_uint<8> weight_short = puppiweight_table[index];
                weight = (int) puppiweight_table[index];
            }
        }

        int pT_new = pfallne[in].hwPt * weight;
        pT_new =  pT_new >> 8;
        pfallne[in].hwPtPuppi = (pt_t) pT_new;
        //std::cout<<in<<" "<<weight<<" "<<pT_new<<" "<<pfallne[in].hwPtPuppi<<std::endl;
    }

}

