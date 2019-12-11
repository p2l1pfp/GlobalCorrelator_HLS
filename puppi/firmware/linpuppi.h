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

#endif
