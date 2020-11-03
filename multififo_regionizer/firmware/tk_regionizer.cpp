#include "regionizer.h"
#include "fifos.h"

void tk_route_link2fifo_unpacked(const TkObj & in, TkObj & center, bool & write_center, TkObj & after, bool & write_after, TkObj & before, bool & write_before) {
    #pragma HSL inline
    bool valid = (in.hwPt != 0);
    bool in_next = in.hwPhi >= +(PFREGION_PHI_SIZE/2-PFREGION_PHI_BORDER);
    bool in_prev = in.hwPhi <= -(PFREGION_PHI_SIZE/2-PFREGION_PHI_BORDER);
    write_center = valid;              center = in; 
    write_after  = valid && (in_next);  after = phiShifted(in, -PFREGION_PHI_SIZE); 
    write_before = valid && (in_prev); before = phiShifted(in, +PFREGION_PHI_SIZE);
}

void tk_route_link2fifo(const PackedTkObj & pin, 
        PackedTkObj & pcenter, bool & write_pcenter, 
        PackedTkObj & pafter,  bool & write_pafter, 
        PackedTkObj & pbefore, bool & write_pbefore) {
    #pragma HSL inline
    TkObj in, center, before, after; 
    bool  write_center, write_before, write_after;
    l1pf_pattern_unpack_one(pin, in); 
    tk_route_link2fifo_unpacked(in, center, write_center, after, write_after, before, write_before);
    pcenter = l1pf_pattern_pack_one(center); write_pcenter = write_center;
    pbefore = l1pf_pattern_pack_one(before); write_pbefore = write_before;
    pafter  = l1pf_pattern_pack_one(after);  write_pafter  = write_after;
}

void tk_route_all_sectors(const PackedTkObj tracks_in[NTKSECTORS][NTKFIBERS], PackedTkObj fifo_in[NTKSECTORS][NTKFIFOS], bool fifo_write[NTKSECTORS][NTKFIFOS]) {
    #pragma HLS inline
    #pragma HLS array_partition variable=tracks_in  complete dim=0
    #pragma HLS array_partition variable=fifo_in  complete dim=0
    #pragma HLS array_partition variable=fifo_write  complete dim=0
    #pragma HLS interface ap_none port=fifo_in
    #pragma HLS interface ap_none port=fifo_write

#if NTKSECTORS == 9
    tk_route_link2fifo(tracks_in[0][0], fifo_in[0][0], fifo_write[0][0], fifo_in[1][2], fifo_write[1][2], fifo_in[8][4], fifo_write[8][4]);
    tk_route_link2fifo(tracks_in[0][1], fifo_in[0][1], fifo_write[0][1], fifo_in[1][3], fifo_write[1][3], fifo_in[8][5], fifo_write[8][5]);
    tk_route_link2fifo(tracks_in[1][0], fifo_in[1][0], fifo_write[1][0], fifo_in[2][2], fifo_write[2][2], fifo_in[0][4], fifo_write[0][4]);
    tk_route_link2fifo(tracks_in[1][1], fifo_in[1][1], fifo_write[1][1], fifo_in[2][3], fifo_write[2][3], fifo_in[0][5], fifo_write[0][5]);
    tk_route_link2fifo(tracks_in[2][0], fifo_in[2][0], fifo_write[2][0], fifo_in[3][2], fifo_write[3][2], fifo_in[1][4], fifo_write[1][4]);
    tk_route_link2fifo(tracks_in[2][1], fifo_in[2][1], fifo_write[2][1], fifo_in[3][3], fifo_write[3][3], fifo_in[1][5], fifo_write[1][5]);
    tk_route_link2fifo(tracks_in[3][0], fifo_in[3][0], fifo_write[3][0], fifo_in[4][2], fifo_write[4][2], fifo_in[2][4], fifo_write[2][4]);
    tk_route_link2fifo(tracks_in[3][1], fifo_in[3][1], fifo_write[3][1], fifo_in[4][3], fifo_write[4][3], fifo_in[2][5], fifo_write[2][5]);
    tk_route_link2fifo(tracks_in[4][0], fifo_in[4][0], fifo_write[4][0], fifo_in[5][2], fifo_write[5][2], fifo_in[3][4], fifo_write[3][4]);
    tk_route_link2fifo(tracks_in[4][1], fifo_in[4][1], fifo_write[4][1], fifo_in[5][3], fifo_write[5][3], fifo_in[3][5], fifo_write[3][5]);
    tk_route_link2fifo(tracks_in[5][0], fifo_in[5][0], fifo_write[5][0], fifo_in[6][2], fifo_write[6][2], fifo_in[4][4], fifo_write[4][4]);
    tk_route_link2fifo(tracks_in[5][1], fifo_in[5][1], fifo_write[5][1], fifo_in[6][3], fifo_write[6][3], fifo_in[4][5], fifo_write[4][5]);
    tk_route_link2fifo(tracks_in[6][0], fifo_in[6][0], fifo_write[6][0], fifo_in[7][2], fifo_write[7][2], fifo_in[5][4], fifo_write[5][4]);
    tk_route_link2fifo(tracks_in[6][1], fifo_in[6][1], fifo_write[6][1], fifo_in[7][3], fifo_write[7][3], fifo_in[5][5], fifo_write[5][5]);
    tk_route_link2fifo(tracks_in[7][0], fifo_in[7][0], fifo_write[7][0], fifo_in[8][2], fifo_write[8][2], fifo_in[6][4], fifo_write[6][4]);
    tk_route_link2fifo(tracks_in[7][1], fifo_in[7][1], fifo_write[7][1], fifo_in[8][3], fifo_write[8][3], fifo_in[6][5], fifo_write[6][5]);
    tk_route_link2fifo(tracks_in[8][0], fifo_in[8][0], fifo_write[8][0], fifo_in[0][2], fifo_write[0][2], fifo_in[7][4], fifo_write[7][4]);
    tk_route_link2fifo(tracks_in[8][1], fifo_in[8][1], fifo_write[8][1], fifo_in[0][3], fifo_write[0][3], fifo_in[7][5], fifo_write[7][5]);
#elif NTKSECTORS == 3
    tk_route_link2fifo(tracks_in[0][0], fifo_in[0][0], fifo_write[0][0], fifo_in[1][2], fifo_write[1][2], fifo_in[2][4], fifo_write[2][4]);
    tk_route_link2fifo(tracks_in[0][1], fifo_in[0][1], fifo_write[0][1], fifo_in[1][3], fifo_write[1][3], fifo_in[2][5], fifo_write[2][5]);
    tk_route_link2fifo(tracks_in[1][0], fifo_in[1][0], fifo_write[1][0], fifo_in[2][2], fifo_write[2][2], fifo_in[0][4], fifo_write[0][4]);
    tk_route_link2fifo(tracks_in[1][1], fifo_in[1][1], fifo_write[1][1], fifo_in[2][3], fifo_write[2][3], fifo_in[0][5], fifo_write[0][5]);
    tk_route_link2fifo(tracks_in[2][0], fifo_in[2][0], fifo_write[2][0], fifo_in[0][2], fifo_write[0][2], fifo_in[1][4], fifo_write[1][4]);
    tk_route_link2fifo(tracks_in[2][1], fifo_in[2][1], fifo_write[2][1], fifo_in[0][3], fifo_write[0][3], fifo_in[1][5], fifo_write[1][5]);
#else
    #error "Unsupported number of sectors"
#endif

}

void tk_router_input_slice(const PackedTkObj tracks_in[NTKSECTORS][NTKFIBERS], PackedTkObj fifo_in[NTKSECTORS][NTKFIFOS], bool fifo_write[NTKSECTORS][NTKFIFOS]) {
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1
    #pragma HLS array_partition variable=tracks_in  complete dim=0
    #pragma HLS array_partition variable=fifo_in  complete dim=0
    #pragma HLS array_partition variable=fifo_write  complete dim=0
    #pragma HLS array_partition variable=fifo_in complete dim=0
    #pragma HLS array_partition variable=fifo_write complete dim=0
    //#pragma HLS interface ap_none port=fifo_in
    //#pragma HLS interface ap_none port=fifo_write

    tk_route_all_sectors(tracks_in, fifo_in, fifo_write);
}

void tk_router_fifo_slice(bool newevent, 
                          const PackedTkObj fifo_in[NTKSECTORS][NTKFIFOS], const bool fifo_write[NTKSECTORS][NTKFIFOS], const bool fifo_full[NTKSECTORS][NTKFIFOS],
                          PackedTkObj fifo_out[NTKSECTORS][NTKFIFOS], bool fifo_out_valid[NTKSECTORS][NTKFIFOS], bool fifo_out_roll[NTKSECTORS][NTKFIFOS])
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

    static rolling_ram_fifo<PackedTkObj::width> fifos[NTKSECTORS*NTKFIFOS];
    #pragma HLS array_partition variable=fifos complete dim=1 // must be 1D array to avoid unrolling also the RAM

    for (int i = 0; i < NTKSECTORS; ++i) {
        #pragma HLS unroll
        for (int j = 0; j < NTKFIFOS; ++j) {
            #pragma HLS unroll
            fifos[i*NTKFIFOS+j].update(newevent, fifo_in[i][j], fifo_write[i][j], fifo_out[i][j], fifo_out_valid[i][j], fifo_full[i][j], fifo_out_roll[i][j]);
        }
    }
}


void tk_router_merge2_slice(const PackedTkObj fifo_out[NTKSECTORS][NTKFIFOS], const bool fifo_out_valid[NTKSECTORS][NTKFIFOS], const bool fifo_out_roll[NTKSECTORS][NTKFIFOS],
        const bool merged_full[NTKSECTORS][NTKFIFOS/2],
        bool fifo_full[NTKSECTORS][NTKFIFOS],
        PackedTkObj merged_out[NTKSECTORS][NTKFIFOS/2], bool merged_out_valid[NTKSECTORS][NTKFIFOS/2], bool merged_out_roll[NTKSECTORS][NTKFIFOS/2])
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

    static fifo_merge2_full<PackedTkObj::width> merger[NTKSECTORS*NTKFIFOS/2];
    #pragma HLS array_partition variable=mergers complete dim=1 

    for (int i = 0; i < NTKSECTORS; ++i) {
        #pragma HLS unroll
        for (int j = 0; j < NTKFIFOS/2; ++j) {
            #pragma HLS unroll
            merger[i*(NTKFIFOS/2)+j].update(fifo_out_roll[i][2*j],
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

void tk_router_merge3_slice(const PackedTkObj merged_out[NTKSECTORS][NTKFIFOS/2], const bool merged_out_valid[NTKSECTORS][NTKFIFOS/2], const bool merged_out_roll[NTKSECTORS][NTKFIFOS/2],
        bool merged_full[NTKSECTORS][NTKFIFOS/2],
        PackedTkObj merged3_out[NTKSECTORS], bool merged3_out_valid[NTKSECTORS], bool merged3_out_roll[NTKSECTORS])
{
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1

    #pragma HLS array_partition variable=merged_out complete dim=0
    #pragma HLS array_partition variable=merged_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged_out_roll  complete dim=0
    #pragma HLS array_partition variable=merged_full  complete dim=0
    #pragma HLS array_partition variable=merged3_out complete dim=0
    #pragma HLS array_partition variable=merged3_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged3_out_roll  complete dim=0

    //#pragma HLS interface ap_none port=merged_full
    //#pragma HLS interface ap_none port=merged3_out
    //#pragma HLS interface ap_none port=merged3_out_valid
    //#pragma HLS interface ap_none port=merged3_out_roll

    static fifo_merge3<PackedTkObj::width> merger[NTKSECTORS];
    #pragma HLS array_partition variable=mergers complete dim=1 

    for (int i = 0; i < NTKSECTORS; ++i) {
        //if (i == 0) merger[i].debug_ = i+1;
        #pragma HLS unroll
        merger[i].update(merged_out_roll[i][0],
                                          merged_out[i][0], merged_out[i][1], merged_out[i][2],
                                          merged_out_valid[i][0], merged_out_valid[i][1], merged_out_valid[i][2],
                                          merged3_out[i], 
                                          merged3_out_valid[i],
                                          merged_full[i][0], merged_full[i][1], merged_full[i][2],
                                          merged3_out_roll[i]);
    }
}

void tk_router_nomerge_output_slice(const PackedTkObj fifo_out[NTKSECTORS][NTKFIBERS], const bool fifo_out_valid[NTKSECTORS][NTKFIBERS],
                            PackedTkObj tracks_out[NTKSECTORS*NTKFIBERS])
{
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1

    #pragma HLS array_partition variable=fifo_out complete dim=0
    #pragma HLS array_partition variable=fifo_out_valid  complete dim=0
    #pragma HLS array_partition variable=fifo_out_roll  complete dim=0
    #pragma HLS array_partition variable=tracks_out complete
    //#pragma HLS interface ap_none port=tracks_out

    for (int i = 0; i < NTKSECTORS; ++i) {
        #pragma HLS unroll
        for (int j = 0; j < NTKFIFOS; ++j) {
            #pragma HLS unroll
            if (fifo_out_valid[i]) {
                tracks_out[i*NTKFIFOS+j] = fifo_out[i];
            } else {
                clear(tracks_out[i*NTKFIFOS+j]);
            }
        }
    }
}

void tk_router_output_slice(const PackedTkObj merged3_out[NTKSECTORS], const bool merged3_out_valid[NTKSECTORS], const bool merged3_out_roll[NTKSECTORS],
                            PackedTkObj tracks_out[NTKSECTORS], bool & newevent_out)
{
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1

    #pragma HLS array_partition variable=merged3_out complete dim=0
    #pragma HLS array_partition variable=merged3_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged3_out_roll  complete dim=0
    #pragma HLS array_partition variable=tracks_out complete
    //#pragma HLS interface ap_none port=tracks_out

    for (int i = 0; i < NTKSECTORS; ++i) {
        if (merged3_out_valid[i]) {
            tracks_out[i] = merged3_out[i];
        } else {
            clear(tracks_out[i]);
        }
    }

    newevent_out = merged3_out_roll[0];
}


bool tk_router(bool newevent, const PackedTkObj tracks_in[NTKSECTORS][NTKFIBERS], PackedTkObj tracks_out[NTKOUT], bool & newevent_out) {
    #pragma HLS pipeline II=1 enable_flush
    #pragma HLS array_partition variable=tracks_in  complete dim=0
    #pragma HLS array_partition variable=tracks_out complete
    //#pragma HLS interface ap_none port=tracks_out

    // all the layers must run in parallel taking as input the output from the previous
    // clock cycle, that is saved in the static variables, and not anything new that is
    // produced during this call of the function.
    // because of this, I call them in reverse order, and write the new "full" flags into a new array,
    // and copy that into the static at the end of this function

    static PackedTkObj fifo_in[NTKSECTORS][NTKFIFOS]; 
    static bool fifo_write[NTKSECTORS][NTKFIFOS], newevent_in;
    #pragma HLS array_partition variable=fifo_in  complete dim=0
    #pragma HLS array_partition variable=fifo_write  complete dim=0

#if defined(ROUTER_NOMERGE)
    tk_router_nomerge_output_slice(fifo_in, fifo_write, tracks_out); newevent_out = newevent_in;
    tk_router_input_slice(tracks_in, fifo_in, fifo_write); newevent_in = newevent;
#else
    static PackedTkObj fifo_out[NTKSECTORS][NTKFIFOS]; 
    static bool fifo_out_valid[NTKSECTORS][NTKFIFOS], fifo_out_roll[NTKSECTORS][NTKFIFOS], fifo_full[NTKSECTORS][NTKFIFOS];
    bool        fifo_full_new[NTKSECTORS][NTKFIFOS];
    #pragma HLS array_partition variable=fifo_out complete dim=0
    #pragma HLS array_partition variable=fifo_out_valid complete dim=0
    #pragma HLS array_partition variable=fifo_out_roll complete dim=0
    #pragma HLS array_partition variable=fifo_full  complete dim=0§
    #pragma HLS array_partition variable=fifo_full_new  complete dim=0§

    static PackedTkObj merged_out[NTKSECTORS][NTKFIFOS/2]; 
    static bool merged_out_valid[NTKSECTORS][NTKFIFOS/2], merged_out_roll[NTKSECTORS][NTKFIFOS/2], merged_full[NTKSECTORS][NTKFIFOS/2];
    bool        merged_full_new[NTKSECTORS][NTKFIFOS/2];
    #pragma HLS array_partition variable=merged_out complete dim=0
    #pragma HLS array_partition variable=merged_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged_out_roll complete dim=0
    #pragma HLS array_partition variable=merged_full  complete dim=0
    #pragma HLS array_partition variable=merged_full_new  complete dim=0

    static PackedTkObj merged3_out[NTKSECTORS];
    static bool merged3_out_valid[NTKSECTORS], merged3_out_roll[NTKSECTORS];
    #pragma HLS array_partition variable=merged3_out complete dim=0
    #pragma HLS array_partition variable=merged3_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged3_out_roll complete dim=0

    tk_router_output_slice(merged3_out, merged3_out_valid, merged3_out_roll, tracks_out, newevent_out);
    tk_router_merge3_slice(merged_out, merged_out_valid, merged_out_roll, merged_full_new, merged3_out, merged3_out_valid, merged3_out_roll);
    tk_router_merge2_slice(fifo_out, fifo_out_valid, fifo_out_roll, merged_full, fifo_full_new, merged_out, merged_out_valid, merged_out_roll);
    tk_router_fifo_slice(newevent_in, fifo_in, fifo_write, fifo_full, fifo_out, fifo_out_valid, fifo_out_roll);
    tk_router_input_slice(tracks_in, fifo_in, fifo_write); newevent_in = newevent;

    for (int is = 0, i = 0; is < NTKSECTORS; ++is) {
        for (int f = 0; f < NTKFIFOS; ++f) fifo_full[is][f] = fifo_full_new[is][f];
        for (int f = 0; f < NTKFIFOS/2; ++f) merged_full[is][f] = merged_full_new[is][f];
    }
#endif

    return true;
}
