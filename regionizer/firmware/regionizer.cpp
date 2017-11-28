#include <cassert>
#include "regionizer.h"

#if 0
template<typename T, int N_OBJ_PER_SECTOR> 
void deserialize_ptsort(hls::stream<T> & instream, T out[N_OBJ_PER_SECTOR]) {
    for (unsigned int io = 0; io < N_OBJ_PER_SECTOR; ++io) {
        T in = instream.read();
        for (int i = N_OBJ_PER_SECTOR-1; i >= 0; --i) {
            if (out[i].hwPt < in.hwPt) {
                if (i == 0 || out[i-1].hwPt >= in.hwPt) {
                    out[i] = in;
                } else {
                    out[i] = out[i-1];
                }            
            } 
        }
    }
}
#endif

template<typename T, int N_OBJ_PER_SECTOR_PER_ETA> 
void push_in_sector_eta(const T & in, T out[N_OBJ_PER_SECTOR_PER_ETA], etaphi_t etaMin, etaphi_t etaMax, etaphi_t etaShift) {
    if (etaMin <= in.hwEta && in.hwEta <= etaMax) {
        for (int i = N_OBJ_PER_SECTOR_PER_ETA-1; i >= 0; --i) {
            if (out[i].hwPt < in.hwPt) {
                if (i == 0 || out[i-1].hwPt >= in.hwPt) {
                    out[i] = in;
                    out[i].hwEta += etaShift;
                } else {
                    out[i] = out[i-1];
                }            
            } 
        }
    }
}

template<typename T, int N_OBJ_PER_SECTOR, int N_OBJ_PER_SECTOR_PER_ETA> 
void deserialize_ptsort_eta(hls::stream<T> instream[N_IN_SECTORS], T deserialized[N_IN_SECTORS][N_OUT_REGIONS_ETA][N_OBJ_PER_SECTOR_PER_ETA]) {
    for (unsigned int is = 0; is < N_IN_SECTORS; ++is) {
        for (unsigned int ie = 0; ie < N_OUT_REGIONS_ETA; ++ie) {
            for (unsigned int io = 0; io < N_OBJ_PER_SECTOR_PER_ETA; ++io) {
                clear(deserialized[is][ie][io]);
            }
        }
    }

    // steam the data in, apply eta cuts and sorting
    for (unsigned int io = 0; io < N_OBJ_PER_SECTOR; ++io) {
        for (unsigned int is = 0; is < N_IN_SECTORS; ++is) {
            T in = instream[is].read();
            for (unsigned int ie = 0; ie < N_OUT_REGIONS_ETA; ++ie) {
                push_in_sector_eta<T,N_OBJ_PER_SECTOR_PER_ETA>(in, deserialized[is][ie], ETA_MIN[ie], ETA_MAX[ie], ETA_SHIFT[ie]);
            } 
        }
    }
}

template<typename T, int N_OBJ_PER_SECTOR_PER_ETA, int N_OBJ_PER_REGION> 
void merge_smart(const T list1[N_OBJ_PER_SECTOR_PER_ETA], etaphi_t phiShift1, const T list2[N_OBJ_PER_SECTOR_PER_ETA], etaphi_t phiShift2, const T list3[N_OBJ_PER_SECTOR_PER_ETA], etaphi_t phiShift3, T out[N_OBJ_PER_REGION]) {
    T tmp[N_OBJ_PER_REGION];
    #pragma HLS ARRAY_PARTITION variable=tmp complete

    for (unsigned int i = 0; i < N_OBJ_PER_SECTOR_PER_ETA; ++i) {
        tmp[i] = list1[i];
        if (list1[i].hwPt > 0) tmp[i].hwPhi += phiShift1;
    }
    for (unsigned int i = N_OBJ_PER_SECTOR_PER_ETA; i < N_OBJ_PER_REGION; ++i) {
        clear(tmp[i]);
    }

    for (unsigned int it = 0; it < N_OBJ_PER_SECTOR_PER_ETA; ++it) {
        for (int iout = N_OBJ_PER_REGION-1; iout >= 0; --iout) {
            if (tmp[iout].hwPt < list2[it].hwPt) {
                if (iout == 0 || tmp[iout-1].hwPt >= list2[it].hwPt) {
                    tmp[iout] = list2[it];
                    if (list2[it].hwPt > 0) tmp[iout].hwPhi += phiShift2;
                } else {
                    tmp[iout] = tmp[iout-1];
                }
            }
        }
    }

    for (unsigned int it = 0; it < N_OBJ_PER_SECTOR_PER_ETA; ++it) {
        for (int iout = N_OBJ_PER_REGION-1; iout >= 0; --iout) {
            if (tmp[iout].hwPt < list3[it].hwPt) {
                if (iout == 0 || tmp[iout-1].hwPt >= list3[it].hwPt) {
                    tmp[iout] = list3[it];
                    if (list3[it].hwPt > 0) tmp[iout].hwPhi += phiShift3;
                } else {
                    tmp[iout] = tmp[iout-1];
                }
            }
        }
    }
    for (unsigned int i = 0; i < N_OBJ_PER_REGION; ++i) {
        out[i] = tmp[i];
    }

 }

void merge_hadcalo(HadCaloObj list1[NCALO_PER_SECTOR_PER_ETA], etaphi_t phiShift1, HadCaloObj list2[NCALO_PER_SECTOR_PER_ETA], etaphi_t phiShift2, HadCaloObj list3[NCALO_PER_SECTOR_PER_ETA], etaphi_t phiShift3, HadCaloObj out[NCALO]) {
   #pragma HLS pipeline II=HLS_pipeline_II
   #pragma HLS array_partition variable=list1 complete 
   #pragma HLS array_partition variable=list2 complete 
   #pragma HLS array_partition variable=list3 complete 
   #pragma HLS array_partition variable=out complete
   merge_smart<HadCaloObj,NCALO_PER_SECTOR_PER_ETA,NCALO>(list1, phiShift1, list2, phiShift2, list3, phiShift3, out); 
}

void regionize_hadcalo(hls::stream<HadCaloObj> fibers[N_IN_SECTORS], HadCaloObj regions[N_OUT_REGIONS][NCALO]) {
   HadCaloObj work[N_IN_SECTORS][N_OUT_REGIONS_ETA][NCALO_PER_SECTOR_PER_ETA];
   deserialize_ptsort_eta<HadCaloObj,NCALO_PER_SECTOR,NCALO_PER_SECTOR_PER_ETA>(fibers, work);

   for (int ie = 0; ie < N_OUT_REGIONS_ETA; ++ie) {
       for (int ip = 0; ip < N_OUT_REGIONS_PHI; ++ip) {
           int ir = N_OUT_REGIONS_PHI*ie + ip;
           merge_hadcalo( work[IN_SECTOR_OF_REGION[ip][0]][ie], PHI_SHIFT[0], work[IN_SECTOR_OF_REGION[ip][1]][ie], PHI_SHIFT[1], work[IN_SECTOR_OF_REGION[ip][2]][ie], PHI_SHIFT[2], regions[ir]);
       }
   }
}
