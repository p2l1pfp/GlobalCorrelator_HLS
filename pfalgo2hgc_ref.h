#ifndef PFALGO2HGC_REF_H
#define PFALGO2HGC_REF_H

#include "firmware/pfalgo2hgc.h"
#include "pfalgo_common_ref.h"

void pfalgo2hgc_ref_set_debug(int debug) ;

void pfalgo2hgc_ref(const HadCaloObj calo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) ;

#endif
