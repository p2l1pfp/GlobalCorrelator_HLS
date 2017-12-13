#include "firmware/bram_hist_vtx.h"
#include <iostream>

void bhv_find_pv_ref(TkObj tracks[BHV_NSECTORS][BHV_NTRACKS], zbin_t & pvbin, z0_t & pv, int &pvsum) {
    ptsum_t histos[BHV_NSECTORS][BHV_NBINS];
    for (unsigned int is = 0; is < BHV_NSECTORS; ++is) {
        for (unsigned int b = 0; b < BHV_NBINS; ++b) histos[is][b] = 0;
        for (unsigned int it = 0; it < BHV_NTRACKS; ++it) {
            pt_t tkpt = (tracks[is][it].hwPt >> 1);
            if (tkpt > 0) {
                if (tkpt > BHV_MAXPT) tkpt = BHV_MAXPT;
                zbin_vt bin = fetch_bin_ref(tracks[is][it].hwZ0);
                assert(bin.bin >= 0 && bin.bin < BHV_NBINS);
                if (bin.valid) {
                    int newsum = int(histos[is][bin.bin]) + tkpt;
                    histos[is][bin.bin] = ptsum_t(newsum > BHV_MAXBIN ? BHV_MAXBIN : newsum);
                }
            }
        }
    }
    int ibest = 0, sbest = 0;
    for (unsigned int b = 0; b < BHV_NBINS; ++b) {
        int sum = 0;
        for (unsigned int is = 0; is < BHV_NSECTORS; ++is) {
            sum += histos[is][b];
            //printf("REF is %2d bin %3d   val %8d sum %8d\n", is, b, int(histos[is][b]), sum);
        }
        //printf("REF bin %3d  sum %8d\n", b, sum);
        if (sum > sbest) { ibest = b; sbest = sum; }
    }
    pvbin = ibest;
    pvsum = sbest;
    pv = bin_center_ref(ibest);
}
