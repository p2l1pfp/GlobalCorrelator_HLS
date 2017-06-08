#ifndef SIMPLE_PFLOW_SIMPLE_PFLOW_H
#define SIMPLE_PFLOW_SIMPLE_PFLOW_H

#include "data.h"

bool match_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, etaphi_t boxSize) ;
etaphi_t dr_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) ;
int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) ;
ap_uint<12> dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<12> max) ;

void simple_pflow_parallel_ref(CaloObj calo[NCALO], TkObj track[NTRACK], PFObj out[NPF]) ;
void simple_pflow_parallel_hwopt(CaloObj calo[NCALO], TkObj track[NTRACK], PFObj out[NPF]) ;

#endif
