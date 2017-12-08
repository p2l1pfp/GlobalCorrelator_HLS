#ifndef SIMPLE_PUPPI_H
#define SIMPLE_PUPPI_H

#include "../../firmware/data.h"
// #include <boost/math/special_functions/gamma.hpp>

//bool match_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, etaphi_t boxSize) ;
//etaphi_t dr_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) ;
int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) ;
//template<int NB> ap_uint<NB>  dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) ;

float weight_function_float( int eToAlpha );
template< class data_T, int N_TABLE > void _lut_puppiweight_init( data_T table_out[N_TABLE] );


void simple_puppi_ref(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne[NNEUTRALS], z0_t Z0, 
                      PFChargedObj out_pfch_pv[NTRACK], PFNeutralObj out_pfne_pv[NNEUTRALS]);

void simple_puppi_hw(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne[NNEUTRALS], z0_t Z0, 
                      PFChargedObj out_pfch_pv[NTRACK], PFNeutralObj out_pfne_pv[NNEUTRALS]);

#define PFALGO3_DR2MAX_TK_CALO 756
#define PFALGO3_DR2MAX_EM_CALO 525
#define PFALGO3_DR2MAX_TK_MU   2101
#define PFALGO3_DR2MAX_TK_EM   84
#define PFALGO3_TK_MAXINVPT    80

#endif
