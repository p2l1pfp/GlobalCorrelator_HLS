#ifndef FIRMWARE_PFALGO2HGC_H
#define FIRMWARE_PFALGO2HGC_H

#ifndef REG_HGCal
#warning "REG_HGCal is not #defined, but this algorithm has only been tested there"
#endif

#include "pfalgo_common.h"

void pfalgo2hgc(const HadCaloObj calo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) ;


#define PFALGO_DR2MAX_TK_CALO 525
#define PFALGO_TK_MAXINVPT_LOOSE    40
#define PFALGO_TK_MAXINVPT_TIGHT    80

#endif
