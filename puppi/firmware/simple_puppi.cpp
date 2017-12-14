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

weight_t puppiweight(int iWeight){
  static weight_t _table[PUPPI_TABLE_SIZE];
  lut_puppiweight_init<weight_t,PUPPI_TABLE_SIZE>(_table);
  return _table[iWeight];
}

int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}

void simple_puppi_hw(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne[NNEUTRALS], z0_t Z0) {

    const z0_t DZMAX = 256;
    const int DR2MAX = 8404; // 0.4 cone

    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS pipeline II=1

    ap_uint<32> eToAlphas[NNEUTRALS];
    ap_uint<8> weights[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=eToAlphas complete
    #pragma HLS ARRAY_PARTITION variable=weights complete    

    for (int in = 0; in < NNEUTRALS; ++in){ 
      eToAlphas[in] = 0; weights[in] = 0; 
    }

    for (int in = 0; in < NNEUTRALS; ++in) {
        
        int sum = 0;
        for (int it = 0; it < NTRACK; ++it) {
            // std::cout << "pfch[it].hwPt = " << pfch[it].hwPt << std::endl;

            if ((Z0 - pfch[it].hwZ0 > DZMAX) || (Z0 - pfch[it].hwZ0 < -DZMAX)) continue; // if track is PV
            int dr2 = dr2_int(pfch[it].hwEta, pfch[it].hwPhi, pfallne[in].hwEta, pfallne[in].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
                ap_uint<9> dr2short = dr2 >> 5; // why?
                int pt2 = pfch[it].hwPt*pfch[it].hwPt;
                int term = (std::min(pt2 >> 5, 131071) << 15)/(dr2short > 0 ? dr2short : ap_uint<9>(1));
                sum += term;
            }
        }    
        eToAlphas[in] = sum >> 10;
    }

    for (int in = 0; in < NNEUTRALS; ++in) {
        if (eToAlphas[in] <= 0){
            weights[in] = 0; // < e^10 where that is the median
        } 
        else if (eToAlphas[in] > 1174){
            weights[in] = 256; // e^14 where that is the median + 2 sigma
        }
        else{
            int sum_short = eToAlphas[in];
            int index = sum_short >= PUPPI_TABLE_SIZE ? PUPPI_TABLE_SIZE : sum_short; // 2^16
            weights[in] = puppiweight(index); //(int) puppiweight_table[index];
        }
    	int ptnew = ( pfallne[in].hwPt * weights[in] ) >> 8;
    	pfallne[in].hwPtPuppi = (pt_t) ptnew;
    }

}


