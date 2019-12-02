#ifndef FIRMWARE_LINPUPPI_H
#define FIRMWARE_LINPUPPI_H

#include <cmath>
#include "../../firmware/data.h"

int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2);

void fwdlinpuppiSum_hw(const HadCaloObj caloin[NCALO], ap_uint<32> sums[NCALO]);
void fwdlinpuppiSum2Pt_hw(const HadCaloObj caloin[NCALO], const ap_uint<32> sums[NCALO], pt_t puppiPts[NCALO]);
void fwdlinpuppiPt_hw(const HadCaloObj caloin[NCALO], pt_t puppiPts[NCALO]);
void fwdlinpuppiNoCrop_hw(const HadCaloObj caloin[NCALO], PFNeutralObj pfallne[NCALO]);

void fwdlinpuppi_hw(const HadCaloObj caloin[NCALO], PFNeutralObj pfselne[NNEUTRALS]);
void fwdlinpuppi_ref(const HadCaloObj caloin[NCALO], PFNeutralObj pfallne[NCALO], PFNeutralObj pfselne[NNEUTRALS]);
void fwdlinpuppi_flt(const HadCaloObj caloin[NCALO], PFNeutralObj pfallne[NCALO], PFNeutralObj pfselne[NNEUTRALS]);


#endif
