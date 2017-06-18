#include "simple_vtx.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

/*
void _lut_zbin_init(ap_uint<16> _table[512]) {
  _table[0] = 0;
  for (int i = 1; i <= 511; ++i) {
    _table[i] = (i / Z0_SCALE)*(NVTXBINS/15.);
    _table[i] > NVTXBINS-1 ? _table[i] = NVTXBINS-1 : _table[i] = _table[i];
  }
}

ap_int<5> z0bin(z0_t iZ0) { 
  //LUT for below
  ap_uint<16> _table[512];
  _lut_zbin_init(_table);
  return _table[iZ0];
}

void _lut_zcbin_init(z0_t _table[NVTXBINS]) {
  _table[0] = 0;
  for (int i = 0; i < NVTXBINS; ++i) {
    _table[i] = (15./NVTXBINS)*(Z0_SCALE*i+Z0_SCALE*0.5);
  }
}
z0_t      convertZ0(z0_t iZ0) { 
  z0_t _table[NVTXBINS];
  _lut_zcbin_init(_table);
  //z0_t  dZ = iZ0*Z0_SCALE*(15/NVTXBINS);
  return _table[iZ0];
}
*/
void fillVtx(TkObj track,  ap_int<5>  vtxbin[NVTXBINS],pt_t       sumbin[NVTXBINS]) { 
  ap_int<5> bin0=z0bin(track.hwZ0);
  pt_t       pt0=track.hwPt;
  pt0 > VTXPTMAX ? pt0=VTXPTMAX     : pt0=pt0;
  (pt0 > 0)      ?  vtxbin[bin0]+=1 : vtxbin[bin0]+=0;
  sumbin[bin0]+=pt0;
}
void simple_vtx_hwopt(TkObj track[NALLTRACK], VtxObj *outvtx) { 
  #pragma HLS RESOURCE variable=track core=RAM_2P_BRAM
  #pragma HLS ARRAY_PARTITION variable=track  complete
  //#pragma HLS INTERFACE axis port=track
  #pragma HLS INTERFACE s_axilite port=return
  #pragma HLS INTERFACE s_axilite port=outvtx
  //#pragma HLS dataflow
  //Fill sum Et and count binned in dZ
  ap_int<5>  vtxbin[NALLTRACK][NVTXBINS];
  pt_t       sumbin[NALLTRACK][NVTXBINS];
  for(int i=0; i < NALLTRACK; ++i) { 
    #pragma HLS unroll
    for(int j=0; j < NVTXBINS; ++j) { 
      vtxbin[i][j]=0;
      sumbin[i][j]=0;
    }
  }
  #pragma HLS ARRAY_PARTITION variable=vtxbin complete dim=1
  #pragma HLS ARRAY_PARTITION variable=sumbin complete dim=1
  //#pragma HLS stream variable=sumbin depth=2
  //#pragma HLS stream variable=vtxbin depth=2
  //#pragma HLS pipeline II=1 rewind
  for(int i=0; i < NALLTRACK; i++) { 
    #pragma HLS unroll
    fillVtx(track[i],vtxbin[i],sumbin[i]);
  }
  #pragma HLS pipeline II=1 rewind
  ap_int<7> pow=1;
  for(int k=0; k < NPOW; ++k) { 
    pow=1<<k;
    for(int i=0; i < NALLTRACK; i+=pow*2) { 
      #pragma HLS unroll
      for(int j=0; j < NVTXBINS; ++j) { 
	vtxbin[i][j]+=vtxbin[i+pow][j];
	sumbin[i][j]+=sumbin[i+pow][j];
      }    
    }
  }
  //Loop through bins and find maximum
  outvtx->hwSumPt = 0;
  for(int it=0; it < NVTXBINS; it++) { 
    if(sumbin[0][it] > outvtx->hwSumPt) { 
      outvtx->hwSumPt = sumbin[0][it];
      outvtx->mult    = vtxbin[0][it];
      outvtx->hwZ0    = it;
    }
  }
  //
  outvtx->hwZ0 = convertZ0(outvtx->hwZ0);
  outvtx->mult > 3 ? outvtx->hwId = 1 : outvtx->hwId = 0;
}

