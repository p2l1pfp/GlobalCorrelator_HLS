#ifndef PFALGO3_REF_H
#define PFALGO3_REF_H

#include "firmware/pfalgo3.h"
#include "pfalgo_common_ref.h"

void pfalgo3_ref_set_debug(int debug) ;

void pfalgo3_em_ref(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const bool isMu[NTRACK], bool isEle[NTRACK], PFNeutralObj outpho[NPHOTON], HadCaloObj hadcalo_out[NCALO]) ;
void pfalgo3_ref(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) ;
    // constants
#endif
