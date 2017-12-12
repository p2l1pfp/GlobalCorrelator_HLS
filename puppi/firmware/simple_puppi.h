#ifndef SIMPLE_PUPPI_H
#define SIMPLE_PUPPI_H

#include <cmath>
#include "../../firmware/data.h"
int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2);

// need to regulate the range of eToAlpha
inline float weight_function_float( float eToAlphaF ){

    int eToAlpha = (int) eToAlphaF;
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

    return weight*256.; // 8 bit number 

}

template<class data_T, int N_TABLE>
void lut_puppiweight_init(data_T table_out[N_TABLE])
{
    for (int ii = 0; ii < N_TABLE; ii++) {
        data_T real_val = weight_function_float( float(ii) );
        // table_out[ii] = (data_T) ii;
        table_out[ii] = real_val;
        // std::cout << "table_out[ii] = " << table_out[ii] << ", " << ii << std::endl;
    }
}


void simple_puppi_ref(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne[NNEUTRALS], z0_t Z0);
void simple_puppi_hw(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne[NNEUTRALS], z0_t Z0);

#define PFALGO3_DR2MAX_TK_CALO 756
#define PFALGO3_DR2MAX_EM_CALO 525
#define PFALGO3_DR2MAX_TK_MU   2101
#define PFALGO3_DR2MAX_TK_EM   84
#define PFALGO3_TK_MAXINVPT    80

#endif
