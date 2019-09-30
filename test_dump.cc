#include <iostream>
#include <string>
#include <string.h>     /* strlen */
#include <stdlib.h>     /* strtoul */
#include <vector>
#include <stdio.h>
#include <stdint.h>
#include <sstream>

using namespace std;

int main(int argc, char **argv)
{

    FILE *fRegionDump;
    fRegionDump = fopen("dummy.dump", "wb");

    uint32_t run = 1;
    uint32_t lumi = 1000;
    uint64_t event = 1;
    fwrite(&run, sizeof(uint32_t), 1, fRegionDump);
    fwrite(&lumi, sizeof(uint32_t), 1, fRegionDump);
    fwrite(&event, sizeof(uint64_t), 1, fRegionDump);

    uint32_t one = 1;
    fwrite(&one, 4, 1, fRegionDump);//number of regions in event
    
    float etaCenter = 0.;
    float etaMin = -1.5;
    float etaMax = 1.5;
    float phiCenter = 0.;
    float phiHalfWidth = 3.14159;
    float etaExtra = 0.;
    float phiExtra = 0.;
    //these should still be floats, not in hw units
    fwrite(&etaCenter, 4, 1, fRegionDump);
    fwrite(&etaMin,    4, 1, fRegionDump);
    fwrite(&etaMax,    4, 1, fRegionDump);
    fwrite(&phiCenter, 4, 1, fRegionDump);
    fwrite(&phiHalfWidth, 4, 1, fRegionDump);
    fwrite(&etaExtra, 4, 1, fRegionDump);
    fwrite(&phiExtra, 4, 1, fRegionDump);

    uint32_t ncalo = 1;
    fwrite(&ncalo, 4, 1, fRegionDump);//number of calo in event

    for (uint32_t ic = 0; ic < ncalo; ic++) {
      int16_t c_hwPt = 16;
      int16_t c_hwEmPt = 0;
      int16_t c_hwPtErr = 2;
      int16_t c_hwEta = -1;
      int16_t c_hwPhi = 1;
      uint16_t c_hwFlags = 0;
      bool c_isEM = false;
      fwrite(&c_hwPt, 2, 1, fRegionDump);
      fwrite(&c_hwEmPt, 2, 1, fRegionDump);
      fwrite(&c_hwPtErr, 2, 1, fRegionDump);
      fwrite(&c_hwEta, 2, 1, fRegionDump);
      fwrite(&c_hwPhi, 2, 1, fRegionDump);
      fwrite(&c_hwFlags, 2, 1, fRegionDump);
      fwrite(&c_isEM, 1, 1, fRegionDump);
      /*c_hwPt = 128;
      c_hwEmPt = 2;
      c_hwPtErr = 16;
      c_hwEta = -100;
      c_hwPhi = 100;
      c_hwFlags = 0;
      c_isEM = false;
      fwrite(&c_hwPt, 2, 1, fRegionDump);
      fwrite(&c_hwEmPt, 2, 1, fRegionDump);
      fwrite(&c_hwPtErr, 2, 1, fRegionDump);
      fwrite(&c_hwEta, 2, 1, fRegionDump);
      fwrite(&c_hwPhi, 2, 1, fRegionDump);
      fwrite(&c_hwFlags, 2, 1, fRegionDump);
      fwrite(&c_isEM, 1, 1, fRegionDump);*/
    }

    uint32_t nemcalo = 1;
    fwrite(&nemcalo, 4, 1, fRegionDump);//number of emcalo in event

    for (uint32_t ie = 0; ie < nemcalo; ie++) {
      int16_t e_hwPt = 128;
      int16_t e_hwEmPt = 128;
      int16_t e_hwPtErr = 16;
      int16_t e_hwEta = -100;
      int16_t e_hwPhi = 100;
      uint16_t e_hwFlags = 0;
      bool e_isEM = true;
      fwrite(&e_hwPt, 2, 1, fRegionDump);
      fwrite(&e_hwEmPt, 2, 1, fRegionDump);
      fwrite(&e_hwPtErr, 2, 1, fRegionDump);
      fwrite(&e_hwEta, 2, 1, fRegionDump);
      fwrite(&e_hwPhi, 2, 1, fRegionDump);
      fwrite(&e_hwFlags, 2, 1, fRegionDump);
      fwrite(&e_isEM, 1, 1, fRegionDump);
    }
    
    uint32_t ntrack = 1;
    fwrite(&ntrack, 4, 1, fRegionDump);

    for (uint32_t it = 0; it < ntrack; it++) {
      uint16_t t_hwInvpt = 4342;
      int32_t  t_hwVtxEta = -89315;
      int32_t  t_hwVtxPhi = 51695;
      bool     t_hwCharge = 0;
      int16_t  t_hwZ0 = 4;
      uint16_t t_hwChi2 = 102;
      uint16_t t_hwStubs = 6;
      uint16_t t_hwFlags = 1;
      int16_t  t_hwPt = 18;
      int16_t  t_hwPtErr = 0;
      int16_t  t_hwCaloPtErr = 11;
      int16_t  t_hwEta = -205;
      int16_t  t_hwPhi = 68;
      fwrite(&t_hwInvpt, 2, 1, fRegionDump);
      fwrite(&t_hwVtxEta, 4, 1, fRegionDump);
      fwrite(&t_hwVtxPhi, 4, 1, fRegionDump);
      fwrite(&t_hwCharge, 1, 1, fRegionDump);
      fwrite(&t_hwZ0, 2, 1, fRegionDump);
      fwrite(&t_hwChi2, 2, 1, fRegionDump);
      fwrite(&t_hwStubs, 2, 1, fRegionDump);
      fwrite(&t_hwFlags, 2, 1, fRegionDump);
      fwrite(&t_hwPt, 2, 1, fRegionDump);
      fwrite(&t_hwPtErr, 2, 1, fRegionDump);
      fwrite(&t_hwCaloPtErr, 2, 1, fRegionDump);
      fwrite(&t_hwEta, 2, 1, fRegionDump);
      fwrite(&t_hwPhi, 2, 1, fRegionDump);
    }

    uint32_t nmuon = 1;
    fwrite(&nmuon, 4, 1, fRegionDump);//number of muon in event

    for (uint32_t im = 0; im < nmuon; im++) {
        int16_t  m_hwPt=10;   
        int16_t  m_hwEta=-80;
        int16_t  m_hwPhi=80;
        uint16_t m_hwFlags=1;
        bool     m_hwCharge=0;
        fwrite(&m_hwPt, 2, 1, fRegionDump);
        fwrite(&m_hwEta, 2, 1, fRegionDump);
        fwrite(&m_hwPhi, 2, 1, fRegionDump);
        fwrite(&m_hwFlags, 2, 1, fRegionDump);
        fwrite(&m_hwCharge, 1, 1, fRegionDump); 
    }

    float z0 = 0.;
    fwrite(&z0, sizeof(float), 1, fRegionDump);
    fwrite(&z0, sizeof(float), 1, fRegionDump);

    float alpha = 1.;
    fwrite(&alpha, sizeof(float), 1, fRegionDump);
    fwrite(&alpha, sizeof(float), 1, fRegionDump);
    fwrite(&alpha, sizeof(float), 1, fRegionDump);
    fwrite(&alpha, sizeof(float), 1, fRegionDump);

    fclose(fRegionDump);

    return 0;
}
