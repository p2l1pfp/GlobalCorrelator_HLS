#ifndef SIMPLE_PUPPI_FORWARD_H
#define SIMPLE_PUPPI_FORWARD_H

#include <cmath>
#include "../../firmware/data.h"

typedef ap_uint<8> weight_t;
// typedef ap_fixed<18,8> weight_t;
#define PUPPI_TABLE_SIZE 1174

//int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2);

inline float erf_approx( float arg ){
    float denom = 1 + 0.278393*arg + 0.230389*arg*arg + 0.000972*arg*arg*arg + 0.078108*arg*arg*arg*arg;
    float weight = 1. - ( 1./pow(denom,4) );
    return weight;
}

// need to regulate the range of eToAlpha
static inline float weight_function_float( float eToAlphaF ){

    float alphaMed = 10.0; // hard-coded for now from phil's studies!
    float alphaRms = 2.0; // hard-coded for now from phil's studies!

    float alpha = log(eToAlphaF+1); // avoid the log(0)!!
    float sigmas = (alpha - alphaMed)*(alpha - alphaMed) / ( alphaRms * alphaRms );
    // float weight = boost::math::gamma_p(0.5,sigmas/2.);
    float weight = erf_approx( sqrt(sigmas/2.) ); // samesies as the line above!
    
    if (alpha < alphaMed) weight = 0.; // signed residual
    return weight*256.; // 8 bit number 

}

template<class data_T, int N_TABLE>
static void lut_puppiweight_init(data_T table_out[N_TABLE])
{
    for (int ii = 0; ii < N_TABLE; ii++) {
        float eToAlpha = float( ii << 10 );
        data_T real_val = (data_T) weight_function_float( eToAlpha );
        // table_out[ii] = (data_T) ii;
        table_out[ii] = real_val;
    }
}

void simple_puppi_forward_ref(PFNeutralObj pfallne[NNEUTRALS]);
//void simple_puppi_forward_hw(PFNeutralObj pfallne[NNEUTRALS], pt_t ptpuppi[NNEUTRALS]);
void simple_puppi_forward_hw(PFNeutralObj pfallne[NNEUTRALS], pt_t ptpuppi[NNEUTRALS], em2calo_dr_t drvals[NPHOTON][NSELCALO]);
void compute_puppi_weight_hw(int index, weight_t &weight);

#endif

#ifndef DRVALSET
#define DRVALSET
#define PFALGO3_DR2MAX_EM_CALO 525
#define PFPUPPI_DR2MAX 8404

//#define PFALGO3_DR2MAX_TK_CALO 27
//#define PFALGO3_DR2MAX_EM_CALO 23
//#define PFALGO3_DR2MAX_TK_MU   45
//#define PFALGO3_DR2MAX_TK_EM   9
#endif
