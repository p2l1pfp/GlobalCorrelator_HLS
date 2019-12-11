#ifndef FIRMWARE_PFALGO3_H
#define FIRMWARE_PFALGO3_H

#ifndef REG_BARREL
#warning "REG_BARREL not defined in PFALGO3: not validated"
#endif

#include "pfalgo_common.h"

void pfalgo3(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) ;

void pfalgo3_set_debug(bool debug);

#define PFALGO_DR2MAX_TK_CALO 1182
#define PFALGO_DR2MAX_EM_CALO 525
#define PFALGO_DR2MAX_TK_EM   84
#define PFALGO_TK_MAXINVPT_LOOSE    40
#define PFALGO_TK_MAXINVPT_TIGHT    80

#endif
