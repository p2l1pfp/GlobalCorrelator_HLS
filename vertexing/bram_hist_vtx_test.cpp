
#include <cstdio>
#include "firmware/bram_hist_vtx.h"
#include "../utils/DiscretePFInputs_IO.h"
#include "../firmware/mp7pf_encoding.h"
#include "../utils/pattern_serializer.h"
#define NTEST 250

int main() {
    DiscretePFInputs inputs("barrel_alltracks_sectors_1x12_TTbar_PU140.dump");
    //DiscretePFInputs inputs("barrel_sectors_1x12_TTbar_PU140.dump");

    TkObj track_in[BHV_NSECTORS][BHV_NTRACKS], track_tmp[2*BHV_NTRACKS];

    TkObj track_in_transposed[BHV_NTRACKS][BHV_NSECTORS];
    MP7DataWord mp7_in[MP7_NCHANN];
    MP7DataWord mp7_out[MP7_NCHANN];

    unsigned int ngood = 0, ntot = 0, ncmssw_good = 0, nagree = 0, nerr = 0;
    double resol = 0;
    double cmssw_resol = 0;
    // run multiple tests

    FILE *test_in  = fopen("bhv_dump_in.txt","w");
    FILE *test_out = fopen("bhv_dump_out.txt","w");
    unsigned int test_in_frame = 0, test_out_frame = 0;
    MP7PatternSerializer serMP7_in("mp7_input.txt",2,1);  
    MP7PatternSerializer serMP7_in_debug("mp7_input_debug.txt",1,0);  
    MP7PatternSerializer serMP7_out("mp7_output.txt",1,0,3); 
    MP7PatternSerializer serMP7_out_debug("mp7_output_debug.txt",1,0); 

    for (int test = 1; test <= NTEST; ++test) {
        // read the event
        if (!inputs.nextEvent()) break;
        if (inputs.event().regions.size() != N_IN_SECTORS) { printf("ERROR: Mismatching number of input regions: %lu\n", inputs.event().regions.size()); return 2; }
        for (int is = 0; is < N_IN_SECTORS; ++is) {
            const Region & r = inputs.event().regions[is];
            dpf2fw::convert<2*BHV_NTRACKS>(r.track, track_tmp); 
            for (int i = 0, t = 0; i < BHV_NTRACKS; ++i) {
                track_in[2*is+0][i] = track_tmp[t++]; 
                track_in[2*is+1][i] = track_tmp[t++]; 
            }
        }

        z0_t   pv_gen = round(inputs.event().genZ0*l1tpf_int::InputTrack::Z0_SCALE);
        z0_t   pv_cmssw = round(inputs.event().z0*l1tpf_int::InputTrack::Z0_SCALE);
        zbin_t pvbin_ref;
        z0_t   pv_ref;
        int    ptsum_ref;
        bhv_find_pv_ref(track_in, pvbin_ref, pv_ref, ptsum_ref);

        for (unsigned int it = 0; it < 18; ++it) {
            for (int go = 1; go >= 0; --go) {
                fprintf(test_in,"Frame %4d : %1d  %3d  %6d\n", test_in_frame++, go, int(fetch_bin_ref(track_in[0][it].hwZ0).bin), int(track_in[0][it].hwPt>>1));
            }
        }
        for (int is = 0; is < BHV_NSECTORS; ++is) { 
            for (int io = 0; io < BHV_NTRACKS; ++io) track_in_transposed[io][is] = track_in[is][io];
        } 
        for (unsigned int ic = 0; ic < N_CLOCKS/2; ++ic) {
            for (unsigned int i = 0; i < MP7_NCHANN; ++i) mp7_in[i] = 0; // clear
            if (ic < BHV_NTRACKS) mp7_pack<BHV_NSECTORS,0>(track_in_transposed[ic], mp7_in);
            serMP7_in(mp7_in);  
        } 
        for (unsigned int ic = 0; ic < N_CLOCKS; ++ic) {
            for (unsigned int i = 0; i < MP7_NCHANN; ++i) mp7_in[i] = 0; // clear
            for (int is = 0; is < BHV_NSECTORS; ++is) {
                zbin_vt bin = fetch_bin_ref(track_in[is][ic/2].hwZ0);
                mp7_in[2*is+0] = bin.bin + (0x10000000 * bin.valid);
                mp7_in[2*is+1] = track_in[is][ic/2].hwPt >> 1;
            } 
            serMP7_in_debug(mp7_in);  
        }

        twoptsums_t hists[BHV_NSECTORS][BHV_NHALFBINS];
        ptsum_t hist[BHV_NBINS];
        for (int is = 0; is < BHV_NSECTORS; ++is) {
             for (unsigned int b = 0; b < BHV_NBINS; ++b) hist[b] = 0;
             for (unsigned int it = 0; it < BHV_NTRACKS; ++it) {
                bhv_add_track(fetch_bin_ref(track_in[is][it].hwZ0), track_in[is][it].hwPt, hist);
                if (it == 17 && is == 0) {
                    ptsum_t max = 0; zbin_t bmax = BHV_NBINS;
                    for (unsigned int b = 0; b < BHV_NBINS; b += 2) {
                       fprintf(test_out,"Frame %4d : %3d   %6d   %6d\n", test_out_frame++, int(b), int(hist[b]), int(hist[b+1]));
                    }
                }
             }
             for (unsigned int b = 0, tb = 0; b < BHV_NBINS; b += 2, ++tb) {
                hists[is][tb] = ( hist[b], hist[b+1] );
             }
        }
        pt_t ptsum_hw;
        zbin_t pvbin_hw = bhv_find_pv(hists, &ptsum_hw);
        z0_t pv_hw = bin_center_ref(pvbin_hw);

        for (unsigned int i = 0; i < MP7_NCHANN; ++i) mp7_out[i] = 0; 
        for (unsigned int ic = 0; ic < N_CLOCKS; ++ic) {
            bool done = (ic == 0);
            mp7_out[0] = done ? 1 : 0;
            mp7_out[1] = done ? ptsum_hw : pt_t(0);
            mp7_out[2] = done ? pvbin_hw : zbin_t(0);
            serMP7_out(mp7_out);
        }
        assert(N_CLOCKS == BHV_NHALFBINS);
        for (unsigned int i = 0; i < MP7_NCHANN; ++i) mp7_out[i] = 0; 
        for (unsigned int ic = 0; ic < N_CLOCKS; ++ic) {
            for (unsigned int is = 0, i = 3; is < BHV_NSECTORS; ++is, i += 2) {
                mp7_out[i+0] = 0x10000 * (ic == 0) + (ic<<1) + ((test % 2) << 24);
                mp7_out[i+1] = (int(hists[is][ic](8,0)) << 16) | int(hists[is][ic](17,9));
            }
            serMP7_out_debug(mp7_out);
        }
        
        printf("GEN PV %+4d    CMSSW PV %+4d  bin %3d :  REF %+4d  bin %3d, ptsum %8d, diff %+4d :  HW %+4d  bin %3d, ptsum %8d diff %+4d  :  MATCH %+1d\n", int(pv_gen), int(pv_cmssw), int(fetch_bin_ref(pv_cmssw).bin), int(pv_ref), int(pvbin_ref), ptsum_ref, int(pv_ref-pv_gen), int(pv_hw), int(pvbin_hw), int(ptsum_hw), int(pv_hw-pv_gen), int(pv_hw-pv_ref));
        if (abs(int(pv_ref-pv_gen)) <= 10) { ngood++; resol += std::pow(double(pv_gen - pv_ref), 2); }
        ntot++;
        if (abs(int(pv_gen-pv_cmssw)) <= 10) { ncmssw_good++; cmssw_resol += std::pow(double(pv_gen - pv_cmssw), 2); }

        if (abs(int(pv_ref-pv_cmssw)) <= 10) nagree++;
        if (pv_ref != pv_hw || ptsum_hw != ptsum_ref) nerr++;
        
    }

    fclose(test_in);
    fclose(test_out);
    printf("Good matches: CMSSW %4d/%4d = %.3f  REF %4d/%4d = %.3f  (REF vs CMSSW: %4d/%4d = %.3f). ERRORS: %d\n", 
        ncmssw_good, ntot, float(ncmssw_good)/(ntot), ngood, ntot, float(ngood)/(ntot), nagree, ntot, float(nagree)/(ntot), nerr);
    printf("Resolution: CMSSW: %5.3f mm   REF %5.3f mm\n", 0.5*sqrt(cmssw_resol/(ncmssw_good)), 0.5*sqrt(resol/(ngood)));
    return 0;
}
