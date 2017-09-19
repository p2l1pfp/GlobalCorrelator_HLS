#include <cstdio>
#include "firmware/simple_vtx.h"

#define NTEST 50
#define NPV   15
#define Z0_SCALE 20.
#define PT_SCALE 4.0
#define ETAPHI_SCALE (4*180/M_PI)

int main() {
  srand(73); // 73 is also a prime number just to be different than Gio ;p
  float zPVSig = 0.7; 
  float zPU    = 4.; 
  int   nSubSector = NALLTRACK/NSECTOR;
  int   nSubPV     = NPV      /NSECTOR;
  TkObj  track[NSECTOR][NALLTRACK]; 
  VtxObj outvtx; VtxObj outvtxref;
  for (int test = 1; test <= NTEST; ++test) {
    float  gaus = float(rand())/RAND_MAX+float(rand())/RAND_MAX;
    float  zPV  = gaus*zPU;
    for(int j  = 0; j < NSECTOR; j++) { 
      int npucharged = (rand() % (nSubSector-nSubPV)/2) + (nSubSector-nSubPV)/2;
      int npvcharged = (rand() % nSubPV/2             ) + nSubPV/2;
      for (int i = 0; i < nSubSector; ++i) {
	track[j][i].hwPt = 0; track[j][i].hwPtErr = 0; track[j][i].hwEta = 0; track[j][i].hwPhi = 0; track[j][i].hwZ0 = 0;
      }
      for(int i  = 0; i < npucharged; i++) { 
	float pt = 20*(rand()/float(RAND_MAX))+2, eta = (rand()/float(RAND_MAX))*2.0-1.0, phi = (rand()/float(RAND_MAX))*2.0-1.0;
	float gaus = float(rand())/float(RAND_MAX)+float(rand())/float(RAND_MAX);
	float z  = gaus*zPU;
	track[j][i].hwPt    = pt * PT_SCALE;
	track[j][i].hwPtErr = (0.2*pt+4) * PT_SCALE; 
	track[j][i].hwEta = eta * ETAPHI_SCALE;
	track[j][i].hwPhi = phi * ETAPHI_SCALE;
	track[j][i].hwZ0  = z * Z0_SCALE;
      }
      for(int i  = npucharged; i < npucharged+npvcharged; i++) { 
	float pt = 100*(rand()/float(RAND_MAX))+2, eta = (rand()/float(RAND_MAX))*2.0-1.0, phi = (rand()/float(RAND_MAX))*2.0-1.0;
	float  pGaus = 0.5*(rand()/RAND_MAX+rand()/RAND_MAX);
	float  pZPV  = (pGaus-0.5)*zPVSig+zPV;
	track[j][i].hwPt    = pt * PT_SCALE;
	track[j][i].hwPtErr = (0.2*pt+4) * PT_SCALE; 
	track[j][i].hwEta = eta * ETAPHI_SCALE;
	track[j][i].hwPhi = phi * ETAPHI_SCALE;
	track[j][i].hwZ0  = pZPV * Z0_SCALE;
      }
    }
    //simple_vtx_hwopt(track[0],track[1],track[2],track[3],track[4],track[5], &outvtx);
    simple_vtx_hwopt(track[0], &outvtx);
    simple_vtx_ref  (track[0], &outvtxref);
    std::cout << "====>  hardware dZ:" << float(outvtx.hwZ0) << " -- emulated dZ: " << float(outvtxref.hwZ0) << " -- True dZ: " << zPV*Z0_SCALE << std::endl;
  }
  return 0;
}
