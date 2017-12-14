#include "simple_puppi.h"
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

typedef ap_uint<7> tk2em_dr_t;
typedef ap_uint<8> weight_t;
typedef ap_uint<10> tk2calo_dr_t;
typedef ap_uint<10> em2calo_dr_t;
typedef ap_uint<12> tk2calo_dq_t;
typedef ap_uint<12> mu2trk_dr_t;

// template<class data_T>
// void lut_puppiweight_init3(data_T table_out[4096])
// {
//   for (int ii = 0; ii < 4096; ii++) {
//     data_T real_val = weight_function_float( float(ii) );
//     table_out[ii] = real_val;
//   }
// }

weight_t puppiweight(weight_t iWeight){
  weight_t _table[4096];
  lut_puppiweight_init<weight_t,4096>(_table);
  return _table[iWeight];
}

int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}

void compute_puppi_weight( ap_uint<32> eToAlpha, ap_uint<8> &curweight ){
	
	if (eToAlpha < 22026){
		curweight = 0; // < e^10 where that is the median
	} 
	else if (eToAlpha > 1202604){
		curweight = 256; // e^14 where that is the median + 2 sigma
	}
	else{
		// bitshift down the sum by 10 bits to get it in the right range
		int sum_short = eToAlpha >> 10;
		int index = sum_short >= 4096 ? 4096 : sum_short; // 2^16
		curweight = puppiweight(index);//(int) puppiweight_table[index];
		// weight = puppiweight_table[index];
	}
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
        std::cout << "sum = " << sum << std::endl;
        eToAlphas[in] = sum;
    }

    for (int in = 0; in < NNEUTRALS; ++in) {
    	
    	compute_puppi_weight( eToAlphas[in], weights[in] );
    	std::cout << "eToAlphas[in] = " << eToAlphas[in] << ", weights[in] " << weights[in]  << std::endl;
    	pfallne[in].eToAlpha = eToAlphas[in];
    	int ptnew = ( pfallne[in].hwPt * weights[in] ) >> 8;
    	pfallne[in].hwPtPuppi = (pt_t) ptnew;

    }

}


