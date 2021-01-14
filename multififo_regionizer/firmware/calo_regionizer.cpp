#include "regionizer.h"
#include "fifos.h"

void calo_route_all_sectors_unpacked(const HadCaloObj calo_in[NCALOSECTORS][NCALOFIBERS], HadCaloObj fifo_in[NCALOSECTORS][NCALOFIFOS], bool fifo_write[NCALOSECTORS][NCALOFIFOS]) {
    #pragma HLS inline
    #pragma HLS array_partition variable=calo_in  complete dim=0
    #pragma HLS array_partition variable=fifo_in  complete dim=0
    #pragma HLS array_partition variable=fifo_write  complete dim=0

    for (int isec = 0; isec < NCALOSECTORS; ++isec) {
        #pragma HLS unroll
        int inxt = (isec == NCALOSECTORS-1 ? 0 : isec + 1);
        // first
        for (int ifib = 0, iout = 0*NCALOFIBERS; ifib < NCALOFIBERS; ++ifib, ++iout) {
            #pragma HLS unroll
            fifo_in[isec][iout]    = calo_in[isec][ifib];
            fifo_write[isec][iout] = calo_in[isec][ifib].hwPt  != 0 && 
                                     calo_in[isec][ifib].hwPhi <= +(PFREGION_PHI_SIZE/2 + PFREGION_PHI_BORDER) && 
                                     calo_in[isec][ifib].hwPhi >= -(PFREGION_PHI_SIZE/2 + PFREGION_PHI_BORDER);
            //if (isec == 0 && ifib == 0) printf("Sector %d, Fiber %d: got pt %d, write %d\n", isec, ifib, calo_in[isec][ifib].hwPt.to_int(), int(fifo_write[isec][iout]));
        }
        // second, from same
        for (int ifib = 0, iout = 1*NCALOFIBERS; ifib < NCALOFIBERS; ++ifib, ++iout) {
            #pragma HLS unroll
            fifo_in[isec][iout]    = phiShifted(calo_in[isec][ifib], -PFREGION_PHI_SIZE);
            fifo_write[isec][iout] = calo_in[isec][ifib].hwPt  != 0 && 
                                     calo_in[isec][ifib].hwPhi >= +(PFREGION_PHI_SIZE/2 - PFREGION_PHI_BORDER);
        }
        // second, from next
        for (int ifib = 0, iout = 2*NCALOFIBERS; ifib < NCALOFIBERS; ++ifib, ++iout) {
            #pragma HLS unroll
            fifo_in[isec][iout]    = phiShifted(calo_in[inxt][ifib], 2*PFREGION_PHI_SIZE);
            fifo_write[isec][iout] = calo_in[inxt][ifib].hwPt  != 0 && 
                                     calo_in[inxt][ifib].hwPhi <= -(3*PFREGION_PHI_SIZE/2 - PFREGION_PHI_BORDER);
        }
        // third, from same
        for (int ifib = 0, iout = 3*NCALOFIBERS; ifib < NCALOFIBERS; ++ifib, ++iout) {
            #pragma HLS unroll
            fifo_in[isec][iout]    = phiShifted(calo_in[isec][ifib], -2*PFREGION_PHI_SIZE);
            fifo_write[isec][iout] = calo_in[isec][ifib].hwPt  != 0 && 
                                     calo_in[isec][ifib].hwPhi >= +(3*PFREGION_PHI_SIZE/2 - PFREGION_PHI_BORDER);
        }
        // third, from next
        for (int ifib = 0, iout = 4*NCALOFIBERS; ifib < NCALOFIBERS; ++ifib, ++iout) {
            #pragma HLS unroll
            fifo_in[isec][iout]    = phiShifted(calo_in[inxt][ifib], PFREGION_PHI_SIZE);
            fifo_write[isec][iout] = calo_in[inxt][ifib].hwPt  != 0 && 
                                     calo_in[inxt][ifib].hwPhi <= -(PFREGION_PHI_SIZE/2 - PFREGION_PHI_BORDER);
        }
    }
}

void calo_route_all_sectors(const PackedCaloObj pcalo_in[NCALOSECTORS][NCALOFIBERS], PackedCaloObj pfifo_in[NCALOSECTORS][NCALOFIFOS], bool pfifo_write[NCALOSECTORS][NCALOFIFOS]) {
    #pragma HLS inline
    #pragma HLS array_partition variable=pcalo_in  complete dim=0
    #pragma HLS array_partition variable=pfifo_in  complete dim=0
    #pragma HLS array_partition variable=pfifo_write  complete dim=0
    
    HadCaloObj calo_in[NCALOSECTORS][NCALOFIBERS], fifo_in[NCALOSECTORS][NCALOFIFOS]; 
    bool fifo_write[NCALOSECTORS][NCALOFIFOS];
    #pragma HLS array_partition variable=calo_in  complete dim=0
    #pragma HLS array_partition variable=fifo_in  complete dim=0
    #pragma HLS array_partition variable=fifo_write  complete dim=0
    for (int isec = 0; isec < NCALOSECTORS; ++isec) {
        #pragma HLS unroll
        for (int ifib = 0; ifib < NCALOFIBERS; ++ifib) {
            #pragma HLS unroll
            l1pf_pattern_unpack_one(pcalo_in[isec][ifib], calo_in[isec][ifib]);
        }
    }
    calo_route_all_sectors_unpacked(calo_in, fifo_in, fifo_write);
    for (int isec = 0; isec < NCALOSECTORS; ++isec) {
        #pragma HLS unroll
        for (int ifib = 0; ifib < NCALOFIFOS; ++ifib) {
            #pragma HLS unroll
            pfifo_in[isec][ifib] = l1pf_pattern_pack_one(fifo_in[isec][ifib]);
            pfifo_write[isec][ifib] = fifo_write[isec][ifib];
        }
    }

}

void calo_router_input_slice(const PackedCaloObj calo_in[NCALOSECTORS][NCALOFIBERS], PackedCaloObj fifo_in[NCALOSECTORS][NCALOFIFOS], bool fifo_write[NCALOSECTORS][NCALOFIFOS]) {
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1
    #pragma HLS array_partition variable=calo_in  complete dim=0
    #pragma HLS array_partition variable=fifo_in  complete dim=0
    #pragma HLS array_partition variable=fifo_write  complete dim=0
    //#pragma HLS interface ap_none port=fifo_in
    //#pragma HLS interface ap_none port=fifo_write

    calo_route_all_sectors(calo_in, fifo_in, fifo_write);
}

void calo_router_fifo_slice(bool newevent, 
                          const PackedCaloObj fifo_in[NCALOSECTORS][NCALOFIFOS], const bool fifo_write[NCALOSECTORS][NCALOFIFOS], const bool fifo_full[NCALOSECTORS][NCALOFIFOS],
                          PackedCaloObj fifo_out[NCALOSECTORS][NCALOFIFOS], bool fifo_out_valid[NCALOSECTORS][NCALOFIFOS], bool fifo_out_roll[NCALOSECTORS][NCALOFIFOS])
{
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1
    #pragma HLS array_partition variable=fifo_in complete dim=0
    #pragma HLS array_partition variable=fifo_write complete dim=0
    #pragma HLS array_partition variable=fifo_full complete dim=0
    #pragma HLS array_partition variable=fifo_out complete dim=0
    #pragma HLS array_partition variable=fifo_out_valid complete dim=0
    #pragma HLS array_partition variable=fifo_out_roll complete dim=0
    //#pragma HLS interface ap_none port=fifo_out
    //#pragma HLS interface ap_none port=fifo_out_valid
    //#pragma HLS interface ap_none port=fifo_out_roll

    static rolling_ram_fifo<PackedCaloObj::width> fifos[NCALOSECTORS*NCALOFIFOS];
    #pragma HLS array_partition variable=fifos complete dim=1 // must be 1D array to avoid unrolling also the RAM

    for (int i = 0; i < NCALOSECTORS; ++i) {
        #pragma HLS unroll
        for (int j = 0; j < NCALOFIFOS; ++j) {
            #pragma HLS unroll
            fifos[i*NCALOFIFOS+j].update(newevent, fifo_in[i][j], fifo_write[i][j], fifo_out[i][j], fifo_out_valid[i][j], fifo_full[i][j], fifo_out_roll[i][j]);
        }
    }
}

void calo_router_merge2_slice(const PackedCaloObj fifo_out[NCALOSECTORS][NCALOFIFOS], const bool fifo_out_valid[NCALOSECTORS][NCALOFIFOS], const bool fifo_out_roll[NCALOSECTORS][NCALOFIFOS],
        const bool merged_full[NCALOSECTORS][NCALOFIFOS/2],
        bool fifo_full[NCALOSECTORS][NCALOFIFOS],
        PackedCaloObj merged_out[NCALOSECTORS][NCALOFIFOS/2], bool merged_out_valid[NCALOSECTORS][NCALOFIFOS/2], bool merged_out_roll[NCALOSECTORS][NCALOFIFOS/2])
{
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1

    #pragma HLS array_partition variable=fifo_full complete dim=0
    #pragma HLS array_partition variable=fifo_out complete dim=0
    #pragma HLS array_partition variable=fifo_out_valid complete dim=0
    #pragma HLS array_partition variable=fifo_out_roll complete dim=0
    #pragma HLS array_partition variable=merged_out complete dim=0
    #pragma HLS array_partition variable=merged_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged_out_roll  complete dim=0
    #pragma HLS array_partition variable=merged_full  complete dim=0

    //#pragma HLS interface ap_none port=fifo_full
    //#pragma HLS interface ap_none port=merged_out
    //#pragma HLS interface ap_none port=merged_out_valid
    //#pragma HLS interface ap_none port=merged_out_roll

    static fifo_merge2_full<PackedCaloObj::width> merger[NCALOSECTORS*NCALOFIFOS/2];
    #pragma HLS array_partition variable=mergers complete dim=1 

    for (int i = 0; i < NCALOSECTORS; ++i) {
        #pragma HLS unroll
        for (int j = 0; j < NCALOFIFOS/2; ++j) {
            //if (i == 0) merger[i*(NCALOFIFOS/2)+j].debug_ = j+1;
            #pragma HLS unroll
            merger[i*(NCALOFIFOS/2)+j].update(fifo_out_roll[i][2*j],
                                          fifo_out[i][2*j], fifo_out[i][2*j+1], 
                                          fifo_out_valid[i][2*j], fifo_out_valid[i][2*j+1], 
                                          merged_full[i][j],  
                                          merged_out[i][j], 
                                          merged_out_valid[i][j],
                                          fifo_full[i][2*j], fifo_full[i][2*j+1], 
                                          merged_out_roll[i][j]);
        }
    }

}

void calo_router_merge4_slice(const PackedCaloObj merged2_out[NCALOSECTORS][NCALOFIFOS/2], const bool merged2_out_valid[NCALOSECTORS][NCALOFIFOS/2], const bool merged2_out_roll[NCALOSECTORS][NCALOFIFOS/2],
        const bool merged_full[NCALOSECTORS][NCALOFIFOS/4],
        bool merged2_full[NCALOSECTORS][NCALOFIFOS/2],
        PackedCaloObj merged_out[NCALOSECTORS][NCALOFIFOS/4], bool merged_out_valid[NCALOSECTORS][NCALOFIFOS/4], bool merged_out_roll[NCALOSECTORS][NCALOFIFOS/4])
{
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1

    #pragma HLS array_partition variable=merged2_full complete dim=0
    #pragma HLS array_partition variable=merged2_out complete dim=0
    #pragma HLS array_partition variable=merged2_out_valid complete dim=0
    #pragma HLS array_partition variable=merged2_out_roll complete dim=0
    #pragma HLS array_partition variable=merged_out complete dim=0
    #pragma HLS array_partition variable=merged_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged_out_roll  complete dim=0
    #pragma HLS array_partition variable=merged_full  complete dim=0

    //#pragma HLS interface ap_none port=merged2_full
    //#pragma HLS interface ap_none port=merged_out
    //#pragma HLS interface ap_none port=merged_out_valid
    //#pragma HLS interface ap_none port=merged_out_roll

    static fifo_merge2_full<PackedCaloObj::width> merger[NCALOSECTORS*NCALOFIFOS/4];
    #pragma HLS array_partition variable=mergers complete dim=1 

    for (int i = 0; i < NCALOSECTORS; ++i) {
        #pragma HLS unroll
        for (int j = 0; j < NCALOFIFOS/4; ++j) {
            //if (i == 0) merger[i*(NCALOFIFOS/2)+j].debug_ = j+1;
            #pragma HLS unroll
            merger[i*(NCALOFIFOS/4)+j].update(merged2_out_roll[i][2*j],
                                          merged2_out[i][2*j], merged2_out[i][2*j+1], 
                                          merged2_out_valid[i][2*j], merged2_out_valid[i][2*j+1], 
                                          merged_full[i][j],  
                                          merged_out[i][j], 
                                          merged_out_valid[i][j],
                                          merged2_full[i][2*j], merged2_full[i][2*j+1], 
                                          merged_out_roll[i][j]);
        }
    }

}


void calo_router_merge_slice(const PackedCaloObj merged4_out[NCALOSECTORS][NCALOFIFOS/4], const bool merged4_out_valid[NCALOSECTORS][NCALOFIFOS/4], const bool merged4_out_roll[NCALOSECTORS][NCALOFIFOS/4],
        bool merged4_full[NCALOSECTORS][NCALOFIFOS/4],
        PackedCaloObj merged_out[NPFREGIONS], bool merged_out_valid[NPFREGIONS], bool merged_out_roll[NPFREGIONS])
{
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1

    #pragma HLS array_partition variable=merged4_out complete dim=0
    #pragma HLS array_partition variable=merged4_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged4_out_roll  complete dim=0
    #pragma HLS array_partition variable=merged4_full  complete dim=0
    #pragma HLS array_partition variable=merged_out complete dim=0
    #pragma HLS array_partition variable=merged_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged_out_roll  complete dim=0

    //#pragma HLS interface ap_none port=merged4_full
    //#pragma HLS interface ap_none port=merged_out
    //#pragma HLS interface ap_none port=merged_out_valid
    //#pragma HLS interface ap_none port=merged_out_roll

    static fifo_merge2_simple<PackedCaloObj::width> merger[NCALOSECTORS*2];
    #pragma HLS array_partition variable=mergers complete dim=1 

    for (int i = 0; i < NCALOSECTORS; ++i) {
        #pragma HLS unroll
        // region with no extra merge (4 -> 2 -> 1 -> 1)
        merged_out[3*i+0]       = merged4_out[i][0];
        merged_out_valid[3*i+0] = merged4_out_valid[i][0];
        merged_out_roll[3*i+0]  = merged4_out_roll[i][0];
        merged4_full[i][0] = false;
        // two regions with extra merge ( 8 -> 4 -> 2 -> 1 )
        for (int j = 0; j <= 1; ++j) {
            #pragma HLS unroll
            merger[2*i+j].update(merged4_out_roll[i][2*j+1],
                                 merged4_out[i][2*j+1],       merged4_out[i][2*j+2], 
                                 merged4_out_valid[i][2*j+1], merged4_out_valid[i][2*j+2], 
                                 merged_out[3*i+j+1], 
                                 merged_out_valid[3*i+j+1],
                                 merged4_full[i][2*j+1],      merged4_full[i][2*j+2], 
                                 merged_out_roll[3*i+j+1]);
;
        }
    }
}


void calo_router_nomerge_output_slice(const PackedCaloObj merged_out[NCALOSECTORS][NCALOFIFOS], const bool merged_out_valid[NCALOSECTORS][NCALOFIFOS],
                            PackedCaloObj calo_out[NCALOSECTORS*NCALOFIFOS])
{
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1

    #pragma HLS array_partition variable=merged_out complete dim=0
    #pragma HLS array_partition variable=merged_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged_out_roll  complete dim=0
    #pragma HLS array_partition variable=calo_out complete
    //#pragma HLS interface ap_none port=calo_out

    for (int i = 0; i < NCALOSECTORS; ++i) {
        for (int j = 0; j < NCALOFIFOS; ++j) {
            if (merged_out_valid[i][j]) {
                calo_out[i*NCALOFIFOS+j] = merged_out[i][j];
            } else {
                clear(calo_out[i*NCALOFIFOS+j]);
            }
        }
    }
}

void calo_router_output_slice(const PackedCaloObj merged_out[NPFREGIONS], const bool merged_out_valid[NPFREGIONS], const bool merged_out_roll[NPFREGIONS],
                            PackedCaloObj calo_out[NCALOOUT], bool & newevent_out)
{
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1

    #pragma HLS array_partition variable=merged_out complete dim=0
    #pragma HLS array_partition variable=merged_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged_out_roll  complete dim=0
    #pragma HLS array_partition variable=calo_out complete
    //#pragma HLS interface ap_none port=calo_out

    for (int i = 0; i < NPFREGIONS; ++i) {
        if (merged_out_valid[i]) {
            calo_out[i] = merged_out[i];
        } else {
            clear(calo_out[i]);
        }
    }

    newevent_out = merged_out_roll[0];
}


bool calo_router(bool newevent, const PackedCaloObj calo_in[NCALOSECTORS][NCALOFIBERS], PackedCaloObj calo_out[NCALOOUT], bool & newevent_out) {
    #pragma HLS pipeline II=1 enable_flush
    #pragma HLS array_partition variable=calo_in  complete dim=0
    #pragma HLS array_partition variable=calo_out complete
    //#pragma HLS interface ap_none port=calo_out

    static PackedCaloObj fifo_in[NCALOSECTORS][NCALOFIFOS]; 
    static bool fifo_write[NCALOSECTORS][NCALOFIFOS], newevent_in;
    #pragma HLS array_partition variable=fifo_in  complete dim=0
    #pragma HLS array_partition variable=fifo_write  complete dim=0

#if defined(ROUTER_NOMERGE)
    calo_router_nomerge_output_slice(fifo_in, fifo_write, calo_out); newevent_out = newevent_in;
    calo_router_input_slice(calo_in, fifo_in, fifo_write); newevent_in = newevent;
#else
    static PackedCaloObj fifo_out[NCALOSECTORS][NCALOFIFOS], merged2_out[NCALOSECTORS][NCALOFIFOS/2], merged4_out[NCALOSECTORS][NCALOFIFOS/4], merged_out[NPFREGIONS];
    static bool fifo_out_valid[NCALOSECTORS][NCALOFIFOS], fifo_out_roll[NCALOSECTORS][NCALOFIFOS], fifo_full[NCALOSECTORS][NCALOFIFOS];
    static bool merged2_out_valid[NCALOSECTORS][NCALOFIFOS/2], merged2_out_roll[NCALOSECTORS][NCALOFIFOS/2], merged2_full[NCALOSECTORS][NCALOFIFOS/2];
    static bool merged4_out_valid[NCALOSECTORS][NCALOFIFOS/4], merged4_out_roll[NCALOSECTORS][NCALOFIFOS/4], merged4_full[NCALOSECTORS][NCALOFIFOS/4];
    static bool merged_out_valid[NPFREGIONS], merged_out_roll[NPFREGIONS];
    #pragma HLS array_partition variable=fifo_full  complete dim=0§
    #pragma HLS array_partition variable=fifo_out complete dim=0
    #pragma HLS array_partition variable=fifo_out_valid complete dim=0
    #pragma HLS array_partition variable=fifo_out_roll complete dim=0
    #pragma HLS array_partition variable=merged2_out complete dim=0
    #pragma HLS array_partition variable=merged2_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged2_out_roll complete dim=0
    #pragma HLS array_partition variable=merged2_full  complete dim=0
    #pragma HLS array_partition variable=merged4_out complete dim=0
    #pragma HLS array_partition variable=merged4_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged4_out_roll complete dim=0
    #pragma HLS array_partition variable=merged4_full  complete dim=0
    #pragma HLS array_partition variable=merged_out complete dim=0
    #pragma HLS array_partition variable=merged_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged_out_roll complete dim=0

    bool fifo_full_new[NCALOSECTORS][NCALOFIFOS], merged2_full_new[NCALOSECTORS][NCALOFIFOS/2], merged4_full_new[NCALOSECTORS][NCALOFIFOS/4];
    #pragma HLS array_partition variable=fifo_full_new  complete dim=0§
    #pragma HLS array_partition variable=merged2_full_new  complete dim=0
    #pragma HLS array_partition variable=merged4_full_new  complete dim=0
   
    calo_router_output_slice(merged_out, merged_out_valid, merged_out_roll, calo_out, newevent_out);
    calo_router_merge_slice(merged4_out, merged4_out_valid, merged4_out_roll, merged4_full_new, merged_out, merged_out_valid, merged_out_roll);
    calo_router_merge4_slice(merged2_out, merged2_out_valid, merged2_out_roll, merged4_full, merged2_full_new, merged4_out, merged4_out_valid, merged4_out_roll);
    calo_router_merge2_slice(fifo_out, fifo_out_valid, fifo_out_roll, merged2_full, fifo_full_new, merged2_out, merged2_out_valid, merged2_out_roll);
    calo_router_fifo_slice(newevent_in, fifo_in, fifo_write, fifo_full, fifo_out, fifo_out_valid, fifo_out_roll);
    calo_router_input_slice(calo_in, fifo_in, fifo_write); newevent_in = newevent;

    for (int is = 0, i = 0; is < NCALOSECTORS; ++is) {
        for (int f = 0; f < NCALOFIFOS; ++f) fifo_full[is][f] = fifo_full_new[is][f];
        for (int f = 0; f < NCALOFIFOS/2; ++f) merged2_full[is][f] = merged2_full_new[is][f];
        for (int f = 0; f < NCALOFIFOS/4; ++f) merged4_full[is][f] = merged4_full_new[is][f];
    }
#endif

    return true;
}






