#ifndef REGIONIZER_H
#define REGIONIZER_H

#include "data.h"
#include "hls_stream.h"

#define N_IN_SECTORS 12
#define N_CLOCKS_SECTOR 18 // number of clock cycles per sector

#define N_OUT_REGIONS_ETA 3
#define N_OUT_REGIONS_PHI 6
#define N_OUT_REGIONS (N_OUT_REGIONS_ETA*N_OUT_REGIONS_PHI)
#define N_SECTORS_PER_PHI_REGION 3 // NOTE: if you change this from 3, you have to also change the code (e.g. merge_sectors)

#define _ETA_025 57    // 0.25    * (4*180/M_PI)
#define _PHI_PIO6 120  // M_PI/6 * (4*180/M_PI)
const int ETA_FID_SIZE = 3*_ETA_025; 
const int ETA_MIN[N_OUT_REGIONS_ETA] = {  -6*_ETA_025, -3*_ETA_025,   _ETA_025 }; 
const int ETA_MAX[N_OUT_REGIONS_ETA] = {   - _ETA_025, +3*_ETA_025, 6*_ETA_025 }; 
const int ETA_SHIFT[N_OUT_REGIONS_ETA] = { +3*_ETA_025, 0, -3*_ETA_025 }; 

// input phi sectors have centers at (iphi+0.5)*phiWidth-M_PI, iphi = 0 .. N_IN_SECTORS-1
// phiWidth = 2*M_PI/N_IN_SECTORS = M_PI/6
// i.e. sector iphi spans in [ 0 - M_PI, iphi * M_PI/6 - M_PI ] in global coordinates
//      local coordinates are always in [ - M_PI/12 , +M_PI/12 ]
// map output phi sectors into input phi sectors
const int IN_SECTOR_OF_REGION[N_OUT_REGIONS_PHI][N_SECTORS_PER_PHI_REGION] = {  
                { 0, 1, 2 }, { 2, 3, 4 }, {4, 5, 6}, { 6, 7, 8 }, {8, 9, 10}, { 10, 11, 0 } 
            };
// by how much I have to shift sector N to make it consistent with the region phi
const int PHI_SHIFT[N_SECTORS_PER_PHI_REGION] = { -_PHI_PIO6, 0, _PHI_PIO6 };
// fiducial size of a region in eta and phi
const int PHI_FID_SIZE = 2*_PHI_PIO6;  

#define NCALO_PER_SECTOR 15
#define NCALO_PER_SECTOR_PER_ETA 8

void regionize_hadcalo(hls::stream<HadCaloObj> fibers[N_IN_SECTORS], HadCaloObj regions[N_OUT_REGIONS][NCALO]) ;
void regionize_hadcalo_ref(hls::stream<HadCaloObj> fibers[N_IN_SECTORS], HadCaloObj regions[N_OUT_REGIONS][NCALO]) ;

#endif
