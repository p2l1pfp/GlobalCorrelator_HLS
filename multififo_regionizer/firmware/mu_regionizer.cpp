#include "regionizer.h"
#include "fifos.h"


void mu_router_input_slice(const glbeta_t etaCenter, const PackedMuObj mu_in[NMUFIBERS], PackedMuObj fifo_in[NPFREGIONS][NMUFIBERS], bool fifo_write[NPFREGIONS][NMUFIBERS]) {
    #pragma HLS pipeline II=1 
    #pragma HLS array_partition variable=mu_in  complete dim=0
    #pragma HLS array_partition variable=fifo_in  complete dim=0
    #pragma HLS array_partition variable=fifo_write  complete dim=0
    //#pragma HLS interface ap_none port=fifo_in
    //#pragma HLS interface ap_none port=fifo_write

    const ap_int<glbphi_t::width+1> INT_MPI = -(PFREGION_PHI_SIZE * 9)/2;
    const ap_int<glbphi_t::width+1> INT_2PI =  (PFREGION_PHI_SIZE * 9);
    const ap_int<glbphi_t::width+1> PHI_HALFSIZE = (PFREGION_PHI_SIZE / 2 + PFREGION_PHI_BORDER);
    const ap_int<glbeta_t::width+1> ETA_HALFSIZE = (PFREGION_ETA_SIZE / 2 + PFREGION_ETA_BORDER);
    for (unsigned int f = 0; f < NMUFIBERS; ++f) {
        #pragma HLS unroll
        GlbMuObj gmu; 
        l1pf_pattern_unpack_one(mu_in[f], gmu);
        for (unsigned int i = 0; i < NPFREGIONS; ++i) {
            #pragma HLS unroll
            ap_int<glbeta_t::width+1> local_eta = gmu.hwEta - etaCenter;
            ap_int<glbphi_t::width+1> local_phi = gmu.hwPhi - ap_int<glbphi_t::width+1>(i * PFREGION_PHI_SIZE);
            if (local_phi < INT_MPI)  local_phi += INT_2PI;
#ifndef __SYNTHESIS__
            //if (gmu.hwPt != 0) printf("hwd push mu ipt %4d  glb eta %+4d phi %+4d -> local  eta %+4d phi %+4d \n",
            //         gmu.hwPt.to_int(), gmu.hwEta.to_int(),  gmu.hwPhi.to_int(), local_eta.to_int(), local_phi.to_int());
#endif
            bool write = gmu.hwPt != 0 && 
                         (local_eta >= -ETA_HALFSIZE && local_eta <= ETA_HALFSIZE) &&
                         (local_phi >= -PHI_HALFSIZE && local_phi <= PHI_HALFSIZE);
            MuObj lmu;
            lmu.hwPt = gmu.hwPt; lmu.hwPtErr = gmu.hwPtErr;
            lmu.hwEta = local_eta; lmu.hwPhi = local_phi;
            fifo_in[i][f]    = l1pf_pattern_pack_one(lmu);
            fifo_write[i][f] = write;
        }
    }
}

void mu_router_fifo_slice(bool newevent, 
                          const PackedMuObj fifo_in[NPFREGIONS][NMUFIBERS], const bool fifo_write[NPFREGIONS][NMUFIBERS], const bool fifo_full[NPFREGIONS][NMUFIBERS],
                          PackedMuObj fifo_out[NPFREGIONS][NMUFIBERS], bool fifo_out_valid[NPFREGIONS][NMUFIBERS], bool fifo_out_roll[NPFREGIONS][NMUFIBERS])
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

    static rolling_ram_fifo<PackedMuObj::width> fifos[NPFREGIONS*NMUFIBERS];
    #pragma HLS array_partition variable=fifos complete dim=1 // must be 1D array to avoid unrolling also the RAM

    for (int i = 0; i < NPFREGIONS; ++i) {
        #pragma HLS unroll
        for (int j = 0; j < NMUFIBERS; ++j) {
            #pragma HLS unroll
            fifos[i*NMUFIBERS+j].update(newevent, fifo_in[i][j], fifo_write[i][j], fifo_out[i][j], fifo_out_valid[i][j], fifo_full[i][j], fifo_out_roll[i][j]);
        }
    }
}


void mu_router_merge_slice(const PackedMuObj fifo_out[NPFREGIONS][NMUFIBERS], const bool fifo_out_valid[NPFREGIONS][NMUFIBERS], const bool fifo_out_roll[NPFREGIONS][NMUFIBERS],
        bool fifo_full[NPFREGIONS][NMUFIBERS],
        PackedMuObj merged_out[NPFREGIONS], bool merged_out_valid[NPFREGIONS], bool merged_out_roll[NPFREGIONS])
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
    //#pragma HLS interface ap_none port=fifo_full
    //#pragma HLS interface ap_none port=merged_out
    //#pragma HLS interface ap_none port=merged_out_valid
    //#pragma HLS interface ap_none port=merged_out_roll

    static fifo_merge2_simple<PackedMuObj::width> merger[NPFREGIONS];
    #pragma HLS array_partition variable=mergers complete dim=1 

    for (int i = 0; i < NPFREGIONS; ++i) {
        #pragma HLS unroll
        merger[i].update(fifo_out_roll[i][0],
                         fifo_out[i][0], fifo_out[i][1], 
                         fifo_out_valid[i][0], fifo_out_valid[i][1], 
                         merged_out[i], 
                         merged_out_valid[i],
                         fifo_full[i][0], fifo_full[i][1], 
                         merged_out_roll[i]);
    }
}

void mu_router_nomerge_output_slice(const PackedTkObj fifo_out[NPFREGIONS][NMUFIBERS], const bool fifo_out_valid[NPFREGIONS][NMUFIBERS],
                            PackedTkObj mu_out[NPFREGIONS*NMUFIBERS])
{
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1

    #pragma HLS array_partition variable=fifo_out complete dim=0
    #pragma HLS array_partition variable=fifo_out_valid  complete dim=0
    #pragma HLS array_partition variable=fifo_out_roll  complete dim=0
    #pragma HLS array_partition variable=mu_out complete
    //#pragma HLS interface ap_none port=mu_out

    for (int i = 0; i < NPFREGIONS; ++i) {
        #pragma HLS unroll
        for (int j = 0; j < NMUFIBERS; ++j) {
            #pragma HLS unroll
            if (fifo_out_valid[i]) {
                mu_out[i*NMUFIBERS+j] = fifo_out[i];
            } else {
                clear(mu_out[i*NMUFIBERS+j]);
            }
        }
    }
}

void mu_router_output_slice(const PackedMuObj merged_out[NPFREGIONS], const bool merged_out_valid[NPFREGIONS], const bool merged_out_roll[NPFREGIONS],
                            PackedMuObj mu_out[NPFREGIONS], bool & newevent_out)
{
    #pragma HLS pipeline II=1 
    #pragma HLS latency min=1 max=1

    #pragma HLS array_partition variable=merged_out complete dim=0
    #pragma HLS array_partition variable=merged_out_valid complete dim=0
    #pragma HLS array_partition variable=merged_out_roll  complete dim=0
    #pragma HLS array_partition variable=mu_out complete
    //#pragma HLS interface ap_none port=mu_out

    for (int i = 0; i < NPFREGIONS; ++i) {
        if (merged_out_valid[i]) {
            mu_out[i] = merged_out[i];
        } else {
            clear(mu_out[i]);
        }
    }

    newevent_out = merged_out_roll[0];
}


bool mu_router(bool newevent, const glbeta_t etaCenter, const PackedMuObj mu_in[NMUFIBERS], PackedMuObj mu_out[NMUOUT], bool & newevent_out) {
    #pragma HLS pipeline II=1 enable_flush
    #pragma HLS array_partition variable=mu_in  complete dim=0
    #pragma HLS array_partition variable=mu_out complete
    //#pragma HLS interface ap_none port=mu_out

    // all the layers must run in parallel taking as input the output from the previous
    // clock cycle, that is saved in the static variables, and not anything new that is
    // produced during this call of the function.
    // because of this, I call them in reverse order, and write the new "full" flags into a new array,
    // and copy that into the static at the end of this function

    static PackedMuObj fifo_in[NPFREGIONS][NMUFIBERS]; 
    static bool fifo_write[NPFREGIONS][NMUFIBERS], newevent_in;
    #pragma HLS array_partition variable=fifo_in  complete dim=0
    #pragma HLS array_partition variable=fifo_write  complete dim=0

#if defined(ROUTER_NOMERGE)
    mu_router_nomerge_output_slice(fifo_in, fifo_write, mu_out); newevent_out = newevent_in;
    mu_router_input_slice(etaCenter, mu_in, fifo_in, fifo_write); newevent_in = newevent;
#else
    static PackedMuObj fifo_out[NPFREGIONS][NMUFIBERS]; 
    static bool fifo_out_valid[NPFREGIONS][NMUFIBERS], fifo_out_roll[NPFREGIONS][NMUFIBERS], fifo_full[NPFREGIONS][NMUFIBERS];
    bool        fifo_full_new[NPFREGIONS][NMUFIBERS];
    #pragma HLS array_partition variable=fifo_out complete dim=0
    #pragma HLS array_partition variable=fifo_out_valid complete dim=0
    #pragma HLS array_partition variable=fifo_out_roll complete dim=0
    #pragma HLS array_partition variable=fifo_full  complete dim=0§
    #pragma HLS array_partition variable=fifo_full_new  complete dim=0§

    static PackedMuObj merged_out[NPFREGIONS]; 
    static bool merged_out_valid[NPFREGIONS], merged_out_roll[NPFREGIONS];
    #pragma HLS array_partition variable=merged_out complete dim=0
    #pragma HLS array_partition variable=merged_out_valid  complete dim=0
    #pragma HLS array_partition variable=merged_out_roll complete dim=0

    mu_router_output_slice(merged_out, merged_out_valid, merged_out_roll, mu_out, newevent_out);
    mu_router_merge_slice(fifo_out, fifo_out_valid, fifo_out_roll, fifo_full_new, merged_out, merged_out_valid, merged_out_roll);
    mu_router_fifo_slice(newevent_in, fifo_in, fifo_write, fifo_full, fifo_out, fifo_out_valid, fifo_out_roll);
    mu_router_input_slice(etaCenter, mu_in, fifo_in, fifo_write); newevent_in = newevent;

    for (int is = 0, i = 0; is < NPFREGIONS; ++is) {
        for (int f = 0; f < NMUFIBERS; ++f) fifo_full[is][f] = fifo_full_new[is][f];
    }
#endif

    return true;
}
