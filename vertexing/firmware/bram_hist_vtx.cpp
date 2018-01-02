#include "bram_hist_vtx.h"

void bhv_add_track(zbin_vt zbin, pt_t tkpt, ptsum_t hist[BHV_NBINS]) {
    #pragma HLS pipeline II=2
    #pragma HLS interface ap_memory port=hist
    #pragma HLS resource  variable=hist core=ram_1p
    if (zbin.valid) {
        pt_t pt = (tkpt >> 1);
        if (pt > BHV_MAXPT) pt = BHV_MAXPT;
        int sum = int(hist[zbin.bin])+pt;
        hist[zbin.bin] = (sum & (BHV_MAXBIN+1)) ? BHV_MAXBIN : sum;
    }
}

zbin_t bhv_find_pv(twoptsums_t hist[BHV_NSECTORS][BHV_NHALFBINS], pt_t *sumpt) {
    #pragma HLS pipeline II=36
    #pragma HLS array_partition variable=hist complete dim=1
    #pragma HLS interface ap_memory port=hist
    #pragma HLS resource  variable=hist core=ram_1p
    zbin_t ibest = 0; int sbest = 0;
    for (unsigned int i = 0; i < BHV_NHALFBINS; ++i) {
        int sum1 = 0, sum2 = 0;
        for (unsigned int is = 0; is < BHV_NSECTORS; ++is) {
            ptsum_t s1 = hist[is][i]( 8, 0);
            ptsum_t s2 = hist[is][i](17, 9);
            sum1 += s1;
            sum2 += s2;
            //printf("HW  is %2d bin %3d   val %8d sum %8d\n", is, 2*i+0, int(s2), sum2);
            //printf("HW  is %2d bin %3d   val %8d sum %8d\n", is, 2*i+1, int(s1), sum1);
        }
        //printf("HW  bin %3d  sum %8d\n", 2*i+0, sum2);
        //printf("HW  bin %3d  sum %8d\n", 2*i+1, sum1);
        if (sum1 > sum2) {
            if (sum1 > sbest) { ibest = 2*i+1; sbest = sum1; }
        } else {
            if (sum2 > sbest) { ibest = 2*i+0; sbest = sum2; }
        }
    }
    *sumpt = sbest;
    return ibest;
}

bool dummy(z0_t z0) { return (z0 > 0); }
