#include "simple_puppi.h"
#include <cmath>
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

// need to regulate the range of eToAlpha
float weight_function_float( int eToAlpha ){

    int shifted_eToAlpha = eToAlpha << 10;
    float shifted_eToAlpha_float = float(shifted_eToAlpha);

    float alphaMed = 10.0; // hard-coded for now from phil's studies!
    float alphaRms = 2.0; // hard-coded for now from phil's studies!

    float alpha = log(shifted_eToAlpha_float);
    if (alpha < alphaMed) return 0.; // signed residual

    float sigmas = (alpha - alphaMed)*(alpha - alphaMed) / ( alphaRms * alphaRms );
    // float weight = boost::math::gamma_p(0.5,sigmas/2.);
    // float weight = erf( sqrt(sigmas/2.) ); // samesies as the line above!
    // need a numerical approximation of erf :(
    float arg = sqrt(sigmas/2.);
    float denom = 1 + 0.278393*arg + 0.230389*arg*arg + 0.000972*arg*arg*arg + 0.078108*arg*arg*arg*arg;
    float weight = 1. - ( 1./pow(denom,4) );

    // printf("sigmas = %f, weight = %f, weight2 = %f, weight3 = %f \n", sigmas, weight, weight2, weight3);
    // printf("eToAlpha = %i, shifted = %f, alpha = %f, weight = %f \n", eToAlpha, shifted_eToAlpha_float, alpha, weight);

    return weight; // 8 bit number 

}

// data_T is a 6-bit number
template< class data_T, int N_TABLE > 
void _lut_puppiweight_init( data_T table_out[N_TABLE] ){

    // alpha = pt^2 / dr^2
    // log alpha range is from like 0 to 12
    // alpha range is from 1 to e^12
    // the index is e^alpha, scanning all possible values
    for (int ii = 0; ii < N_TABLE; ii++){
        float curweight = weight_function_float( ii ); // ii = e^alpha
        int weight_integer = int(curweight * 256.);
        table_out[ii] = weight_integer;
        // printf("weight = %f, weight_integer = %i \n", curweight, weight_integer);
    }
}

void simple_puppi_hw(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne[NNEUTRALS], z0_t Z0, 
                      PFChargedObj out_pfch_pv[NTRACK], PFNeutralObj out_pfne_pv[NNEUTRALS]) {

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
                int index = sum_short >= 4096 ? 4096 : sum; // 2^16
                // ap_uint<8> weight_short = puppiweight_table[index];            
                // ap_uint<6> weight = puppiweight_table[1000];            
                // weight = (int) weight_short;
                // weight = 0;
                weight = (int) puppiweight_table[sum_short];
            }
        }

        int pT_new = ( pfallne[in].hwPt * weight ) >> 8;
        pfallne[in].hwPtPuppi = (pt_t) pT_new;
    }

}


