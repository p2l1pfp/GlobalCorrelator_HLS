#include "../firmware/data.h"
#include "firmware/simple_vtx.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif
#define Z0_SCALE 20.

void simple_vtx_ref(TkObj track[NALLTRACK], VtxObj *outvtx) { 
  //Fill sum Et and count binned in dZ
  ap_int<5>  vtxbin[NVTXBINS];
  pt_t       sumbin[NVTXBINS];
  for (int ivtx = 0; ivtx < NVTXBINS; ++ivtx) { vtxbin[ivtx] = 0;  sumbin[ivtx] = 0;}
  //  for(int isec=0; isec < NSECTOR;    isec++) { 
  for(int it=0; it < NALLTRACK; it++) { 
    ap_int<5> bin=((track[it].hwZ0))*1./Z0_SCALE*(NVTXBINS/15);
    if(bin < 0)        bin = 0;
    if(bin > NVTXBINS) bin = NVTXBINS;
    pt_t       pt =track[it].hwPt;
    pt > VTXPTMAX ? pt=VTXPTMAX : pt=pt;
    vtxbin[bin]++;
    sumbin[bin]+=pt;
  }
  //Loop through bins and find maximum
  outvtx->hwSumPt = 0;
  outvtx->mult    = 0;
  outvtx->hwZ0    = 0;
  for(int it=0; it < NVTXBINS; it++) { 
    if(sumbin[it] > outvtx->hwSumPt) { 
      outvtx->hwSumPt = sumbin[it];
      outvtx->mult    = vtxbin[it];
      outvtx->hwZ0    = it*Z0_SCALE*(NVTXBINS/15)+Z0_SCALE*0.5;
    }
  }
  outvtx->mult > 3 ? outvtx->hwId = 1 : outvtx->hwId = 0;
}

