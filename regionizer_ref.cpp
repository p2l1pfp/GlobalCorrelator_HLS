#include "firmware/regionizer.h"
#include "test_utils.h"

template<typename T, int N_OBJ_PER_SECTOR_PER_ETA> 
void push_in_sector_ref(const T & in, T out[N_OBJ_PER_SECTOR_PER_ETA], etaphi_t etaMin, etaphi_t etaMax, etaphi_t etaShift) {
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

template<typename T, int N_OBJ_PER_SECTOR_PER_ETA, int N_OBJ_PER_REGION> 
void merge_sectors_ref(T list1[N_OBJ_PER_SECTOR_PER_ETA], etaphi_t phiShift1, T list2[N_OBJ_PER_SECTOR_PER_ETA], etaphi_t phiShift2, T list3[N_OBJ_PER_SECTOR_PER_ETA], etaphi_t phiShift3, T out[N_OBJ_PER_REGION]) {
    unsigned int i1 = 0, i2 = 0, i3 = 0;
    for (unsigned int i = 0; i < N_OBJ_PER_REGION; ++i) {
        clear(out[i]);
    }
    for (unsigned int i = 0; i < N_OBJ_PER_REGION; ++i) {
        pt_t pt1 = (i1 < N_OBJ_PER_SECTOR_PER_ETA ? list1[i1].hwPt : pt_t(0));
        pt_t pt2 = (i2 < N_OBJ_PER_SECTOR_PER_ETA ? list2[i2].hwPt : pt_t(0));
        pt_t pt3 = (i3 < N_OBJ_PER_SECTOR_PER_ETA ? list3[i3].hwPt : pt_t(0));
        if (pt1 >= pt2 && pt1 >= pt3) {
            if (pt1 == 0) break; // if i'm picking a null, it means I'm out of objects
            out[i] = list1[i1]; 
            out[i].hwPhi += phiShift1;
            i1++;  
        } else if (pt1 < pt2 && pt2 >= pt3)  {
            assert(i2 < N_OBJ_PER_SECTOR_PER_ETA);
            out[i] = list2[i2]; 
            out[i].hwPhi += phiShift2;
            i2++;  
        } else {
            assert(i3 < N_OBJ_PER_SECTOR_PER_ETA);
            out[i] = list3[i3]; 
            out[i].hwPhi += phiShift3;
            i3++;  
        }
    }
}

template<typename T, int N_OBJ_PER_SECTOR, int N_OBJ_PER_SECTOR_PER_ETA, int N_OBJ_PER_REGION> 
void regionize_ref(hls::stream<T> instream[N_IN_SECTORS], T regions[N_OUT_REGIONS][N_OBJ_PER_REGION]) {

    // define and clear work area for deserializers 
    static T obj_sector_eta[N_IN_SECTORS][N_OUT_REGIONS_ETA][N_OBJ_PER_SECTOR_PER_ETA];
    for (unsigned int is = 0; is < N_IN_SECTORS; ++is) {
        for (unsigned int ie = 0; ie < N_OUT_REGIONS_ETA; ++ie) {
            for (unsigned int io = 0; io < N_OBJ_PER_SECTOR_PER_ETA; ++io) {
                clear(obj_sector_eta[is][ie][io]);
            }
        }
    }
    // steam the data in
    for (unsigned int io = 0; io < N_OBJ_PER_SECTOR; ++io) {
        for (unsigned int is = 0; is < N_IN_SECTORS; ++is) {
            T in = instream[is].read();
            //if (in.hwPt > 0) printf("Object %d in sector %d, ieta %+3d (eta %+.3f) bounds [%+3d,%+3d]; [%+3d,%+3d]; [%+3d,%+3d]\n", io, is, int(in.hwEta), in.hwEta*0.25/_ETA_025, ETA_MIN[0], ETA_MAX[0], ETA_MIN[1], ETA_MAX[1], ETA_MIN[2], ETA_MAX[2]);
            for (unsigned int ie = 0; ie < N_OUT_REGIONS_ETA; ++ie) {
                push_in_sector_ref<T,N_OBJ_PER_SECTOR_PER_ETA>(in, obj_sector_eta[is][ie], ETA_MIN[ie], ETA_MAX[ie], ETA_SHIFT[ie]);
                //if (in.hwPt > 0) printf("   now sector %d eta %d (found %u)\n", is, ie, count_nonzero(obj_sector_eta[is][ie], NCALO_PER_SECTOR_PER_ETA));
            }
        }
    }

    //for (unsigned int is = 0; is < N_IN_SECTORS; ++is) {
    //    for (unsigned int ie = 0; ie < N_OUT_REGIONS_ETA; ++ie) {
    //        printf("SECTOR %d ETA %d (FOUND %u)\n", is, ie, count_nonzero(obj_sector_eta[is][ie], NCALO_PER_SECTOR_PER_ETA));
    //    }
    //}

    // clear the output
    for (unsigned int ir = 0; ir < N_OUT_REGIONS; ++ir) {
        for (unsigned int io = 0; io < N_OBJ_PER_REGION; ++io) {
            clear(regions[ir][io]);
        }
    }

    // do the merging of the (sectors x eta) into regions
    for (unsigned int ie = 0; ie < N_OUT_REGIONS_ETA; ++ie) {
        for (unsigned int ip = 0; ip < N_OUT_REGIONS_PHI; ++ip) {
            unsigned int ir = N_OUT_REGIONS_PHI*ie + ip;
            merge_sectors_ref<T,N_OBJ_PER_SECTOR_PER_ETA,N_OBJ_PER_REGION>(
                obj_sector_eta[IN_SECTOR_OF_REGION[ip][0]][ie], PHI_SHIFT[0],
                obj_sector_eta[IN_SECTOR_OF_REGION[ip][1]][ie], PHI_SHIFT[1],
                obj_sector_eta[IN_SECTOR_OF_REGION[ip][2]][ie], PHI_SHIFT[2],
                regions[ir]);
        }
    }
}

void regionize_hadcalo_ref(hls::stream<HadCaloObj> fibers[N_IN_SECTORS], HadCaloObj regions[N_OUT_REGIONS][NCALO]) {
   regionize_ref<HadCaloObj,NCALO_PER_SECTOR,NCALO_PER_SECTOR_PER_ETA,NCALO>(fibers, regions); 
}


