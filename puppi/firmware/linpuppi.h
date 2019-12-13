#ifndef FIRMWARE_LINPUPPI_H
#define FIRMWARE_LINPUPPI_H

#include <cmath>
#include "../../firmware/data.h"

int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2);

// charged
void linpuppi_chs(z0_t pvZ0, const PFChargedObj pfch[NTRACK], PFChargedObj outallch[NTRACK]) ;

// neutrals, in the tracker
void linpuppiNoCrop(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], PFNeutralObj outallne[NALLNEUTRALS]) ;
void linpuppi(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], PFNeutralObj outselne[NNEUTRALS]) ;

// neutrals, forward
void fwdlinpuppi(const HadCaloObj caloin[NCALO], PFNeutralObj pfselne[NNEUTRALS]);
void fwdlinpuppiNoCrop(const HadCaloObj caloin[NCALO], PFNeutralObj pfallne[NCALO]);

void fwdlinpuppi_set_debug(bool debug);

#define LINPUPPI_ptLSB 0.25
#define LINPUPPI_DR2LSB 1.9e-5
#define LINPUPPI_dzLSB  0.05
#define LINPUPPI_pt2LSB LINPUPPI_ptLSB*LINPUPPI_ptLSB
#define LINPUPPI_pt2DR2_scale LINPUPPI_ptLSB*LINPUPPI_ptLSB/LINPUPPI_DR2LSB

//=================================================
#if defined(REG_Barrel)

#define LINPUPPI_DR2MAX  4727 // 0.3 cone
#define LINPUPPI_DR2MIN   257 // 0.07 cone
#define LINPUPPI_dzCut     10
#define LINPUPPI_ptMax    200 // 50.0/LINPUPPI_ptLSB 

#define LINPUPPI_ptSlopeNe  0.3
#define LINPUPPI_ptSlopePh  0.3
#define LINPUPPI_ptZeroNe   4.0
#define LINPUPPI_ptZeroPh   2.5
#define LINPUPPI_alphaSlope 0.7
#define LINPUPPI_alphaZero  6.0
#define LINPUPPI_alphaCrop  4.0
#define LINPUPPI_priorNe    5.0
#define LINPUPPI_priorPh    1.0

#define LINPUPPI_ptCut        4 // 1.0/LINPUPPI_ptLSB

//=================================================
#elif defined(REG_HGCal) 

#error "Not implemented"
#define LINPUPPI_dzCut     40

//=================================================
#elif defined(REG_HGCalNoTK)

#define LINPUPPI_DR2MAX  4727 // 0.3 cone
#define LINPUPPI_DR2MIN    84 // 0.04 cone
#define LINPUPPI_dzCut     40 // unused
#define LINPUPPI_ptMax    200 // 50.0/LINPUPPI_ptLSB 

#define LINPUPPI_ptSlopeNe  0.3
#define LINPUPPI_ptSlopePh  0.4
#define LINPUPPI_ptZeroNe   9.0
#define LINPUPPI_ptZeroPh   5.0
#define LINPUPPI_alphaSlope 2.2
#define LINPUPPI_alphaZero  9.0
#define LINPUPPI_alphaCrop  4.0
#define LINPUPPI_priorNe    7.0
#define LINPUPPI_priorPh    5.0

#define LINPUPPI_ptCut       16 // 4.0/LINPUPPI_ptLSB

//=================================================
#elif defined(REG_HF)

#define LINPUPPI_DR2MAX  4727 // 0.3 cone
#define LINPUPPI_DR2MIN   525 // 0.1 cone
#define LINPUPPI_dzCut     40 // unused

#define LINPUPPI_ptMax    400 // 100.0/LINPUPPI_ptLSB 

#define LINPUPPI_ptSlopeNe  0.25
#define LINPUPPI_ptSlopePh  0.25
#define LINPUPPI_ptZeroNe   14.
#define LINPUPPI_ptZeroPh   14.
#define LINPUPPI_alphaSlope 0.6
#define LINPUPPI_alphaZero  9.0
#define LINPUPPI_alphaCrop  4.0
#define LINPUPPI_priorNe    6.0
#define LINPUPPI_priorPh    6.0

#define LINPUPPI_ptCut      40  // 10.0/LINPUPPI_ptLSB

#endif

#endif
