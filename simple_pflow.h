#ifndef SIMPLE_PFLOW_SIMPLE_PFLOW_H
#define SIMPLE_PFLOW_SIMPLE_PFLOW_H

#include "data.h"

bool match_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, etaphi_t boxSize) ;

void simple_pflow_iterative_ref(CaloObj calo[NCALO], TkObj track[NTRACK], PFObj out[NPF]) ;
void simple_pflow_iterative_hwopt(CaloObj calo[NCALO], TkObj track[NTRACK], PFObj out[NPF]) ;

void simple_pflow_parallel_ref(CaloObj calo[NCALO], TkObj track[NTRACK], PFObj out[NPF]) ;
void simple_pflow_parallel_hwopt(CaloObj calo[NCALO], TkObj track[NTRACK], PFObj out[NPF]) ;

#endif
