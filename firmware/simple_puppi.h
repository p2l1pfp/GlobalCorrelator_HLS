#ifndef SIMPLE_PUPPI_H
#define SIMPLE_PUPPI_H

#include "data.h"

//bool match_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, etaphi_t boxSize) ;
//etaphi_t dr_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) ;
int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) ;
//template<int NB> ap_uint<NB>  dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) ;

void simple_puppi_ref(PFChargedObj pfch[NTRACK]);
void simple_puppi_hw(PFChargedObj pfch[NTRACK]);

#define PFALGO3_DR2MAX_TK_CALO 756
#define PFALGO3_DR2MAX_EM_CALO 525
#define PFALGO3_DR2MAX_TK_MU   2101
#define PFALGO3_DR2MAX_TK_EM   84
#define PFALGO3_TK_MAXINVPT    80

#endif
