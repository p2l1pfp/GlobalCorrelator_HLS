#ifndef FIRMWARE_LINPUPPI_H
#define FIRMWARE_LINPUPPI_H

#include <cmath>
#include "../../firmware/data.h"

int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2);


void fwdlinpuppi(const HadCaloObj caloin[NCALO], PFNeutralObj pfselne[NNEUTRALS]);
void fwdlinpuppiNoCrop(const HadCaloObj caloin[NCALO], PFNeutralObj pfallne[NCALO]);

void fwdlinpuppi_ref(const HadCaloObj caloin[NCALO], PFNeutralObj pfallne[NCALO], PFNeutralObj pfselne[NNEUTRALS], bool debug);
void fwdlinpuppi_flt(const HadCaloObj caloin[NCALO], PFNeutralObj pfallne[NCALO], PFNeutralObj pfselne[NNEUTRALS], bool debug);

void fwdlinpuppi_set_debug(bool debug);

#define LINPUPPI_ptLSB 0.25
#define LINPUPPI_DR2LSB 1.9e-5
#define LINPUPPI_pt2LSB LINPUPPI_ptLSB*LINPUPPI_ptLSB
#define LINPUPPI_pt2DR2_scale LINPUPPI_ptLSB*LINPUPPI_ptLSB/LINPUPPI_DR2LSB

#if defined(REG_HGCalNoTK)

#define LINPUPPI_DR2MAX  4727 // 0.3 cone
#define LINPUPPI_DR2MIN    84 // 0.04 cone
#define LINPUPPI_ptMax    50.0 // 
#define LINPUPPI_ptSlopeNe  0.3
#define LINPUPPI_ptSlopePh  0.4
#define LINPUPPI_ptZeroNe   9.0
#define LINPUPPI_ptZeroPh   5.0
#define LINPUPPI_alphaSlope 2.2
#define LINPUPPI_alphaZero  9.0
#define LINPUPPI_alphaCrop  4.0
#define LINPUPPI_priorNe    7.0
#define LINPUPPI_priorPh    5.0
#define LINPUPPI_ptCut      4.0

#endif

#endif
