#include "firmware/regionizer.h"
#include "firmware/obj_unpackers.h"
#include "../utils/pattern_serializer.h"

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>

#define TLEN REGIONIZERNCLOCKS 

bool readEventTk(FILE *file, std::vector<TkObj> inputs[NTKSECTORS][NTKFIBERS], uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventCalo(FILE *file, std::vector<HadCaloObj> inputs[NCALOSECTORS][NCALOFIBERS], bool zside, uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventMu(FILE *file, std::vector<GlbMuObj> inputs[NMUFIBERS], uint32_t &run, uint32_t &lumi, uint64_t &event) ;
 

int main(int argc, char **argv) {
    std::string sample = "TTbar_PU200";
    FILE *fMC_calo  = fopen(("caloDump_hgcal."+sample+".txt").c_str(), "r");
    FILE *fMC_tk  = fopen(("trackDump_hgcalPos."+sample+".txt").c_str(), "r");
    FILE *fMC_mu  = fopen(("muonDump_all."+sample+".txt").c_str(), "r");
    if (!fMC_calo || !fMC_tk || !fMC_mu) {
        printf("Couldn't open input files\n");
        return 2;
    }
    const glbeta_t etaCenter = 2*PFREGION_ETA_SIZE; // eta = +2.0

    PatternSerializer serPatternsIn("input-emp.txt"), serPatternsOut("output-emp.txt"), serPatternsRef("output-emp-ref.txt");
    assert(NTKFIBERS == 2 && NCALOFIBERS == 4 && NMUFIBERS == 2);
    assert(PACKING_NCHANN >= NTKSECTORS*NTKFIBERS/2*3 + NCALOSECTORS*NCALOFIBERS*3 + NMUFIBERS/2*3);
    ap_uint<64> all_channels_in[PACKING_NCHANN], all_channels_out[PACKING_NCHANN], all_channels_ref[PACKING_NCHANN];
    bool        all_valid_in[PACKING_NCHANN], all_valid_out[PACKING_NCHANN], all_valid_ref[PACKING_NCHANN];

    for (unsigned int i = 0; i < PACKING_NCHANN; ++i) {
        all_channels_in[i] = 0; all_channels_out[i] = 0; all_channels_ref[i] = 0;   
        all_valid_in[i] = 0; all_valid_out[i] = 0; all_valid_ref[i] = 0;   
    }

    bool ok = true; unsigned int frame = 0;

    for (int itest = 0; itest < 10; ++itest) {
        std::vector<TkObj>      tk_inputs[NTKSECTORS][NTKFIBERS];
        std::vector<HadCaloObj> calo_inputs[NCALOSECTORS][NCALOFIBERS];
        std::vector<GlbMuObj>   mu_inputs[NMUFIBERS];
        std::vector<std::pair<z0_t,pt_t>> vtx_inputs;

        uint32_t run = 0, lumi = 0; uint64_t event = 0;
        if (!readEventTk(fMC_tk, tk_inputs, run, lumi, event) || 
            !readEventCalo(fMC_calo, calo_inputs, /*zside=*/true, run, lumi, event) ||
            !readEventMu(fMC_mu, mu_inputs, run, lumi, event)) break;

        std::vector<GlbMuObj> mu_flat;
        for (int i = 0, nmu = mu_inputs[0].size() + mu_inputs[1].size(); i < nmu; ++i) {
            mu_flat.push_back(mu_inputs[i%2][i/2]);
        }

        for (int i = 0; i < TLEN; ++i, ++frame) {
            unsigned int ilink = 0, iref = 0;

            for (int s = 0; s < NTKSECTORS; ++s) {
                ap_uint<64> tk[NTKFIBERS];
                bool       tkv[NTKFIBERS];
                for (int f = 0; f < NTKFIBERS; ++f) {
                    tk[f] = 0; tkv[f] = 0;
                    if (i < TLEN-1 && i < int(tk_inputs[s][f].size())) { // emp protocol, must leave one null frame at the end
                        tk[f]  = l1pf_pattern_pack_one(tk_inputs[s][f][i]);
                        tkv[f] = 1;
                    }
                }

                ap_uint<96> w1 = tk[0], w2 = tk[1];
                all_channels_in[ilink+0](63, 0) = w1(95, 32); all_valid_in[ilink+0] = tkv[0];
                all_channels_in[ilink+1](63,32) = w1(31,  0); all_valid_in[ilink+1] = tkv[0] || tkv[1];
                all_channels_in[ilink+1](31, 0) = w2(95, 64);
                all_channels_in[ilink+2](63, 0) = w2(63,  0); all_valid_in[ilink+2] = tkv[1];
                ilink += 3;

                all_channels_ref[iref] = tk[0]; all_valid_ref[iref++] = tkv[0];
                all_channels_ref[iref] = tk[1]; all_valid_ref[iref++] = tkv[1];
            }

            for (int s = 0; s < NCALOSECTORS; ++s) {
                for (int f = 0; f < NCALOFIBERS; ++f) {
                    ap_uint<64> calo = 0; bool cv = 0;
                    if (i < TLEN-1 && i < int(calo_inputs[s][f].size())) { // emp protocol, must leave one null frame at the end
                        calo = l1pf_pattern_pack_one(calo_inputs[s][f][i]);
                        cv   = 1;
                    }
                    all_channels_in[ilink] =    0; all_valid_in[ilink++] = cv; 
                    all_channels_in[ilink] = calo; all_valid_in[ilink++] = cv; 
                    all_channels_in[ilink] =    0; all_valid_in[ilink++] =  0;

                    all_channels_ref[iref] = calo; all_valid_ref[iref++] = cv;
                }
            }

            ap_uint<64> mu[NMUFIBERS]; static ap_uint<64> mu_queue;
            unsigned int imu = i*3/2, nmu = mu_flat.size();
            if (i % 2 == 0) {
                bool muv = imu < nmu, muv2 = imu+1 < nmu;
                ap_uint<64> mu = (muv ? l1pf_pattern_pack_one(mu_flat[imu]) : ap_uint<64>(0));
                all_channels_in[ilink] =  0; all_valid_in[ilink++] = muv;
                all_channels_in[ilink] = mu; all_valid_in[ilink++] = muv;
                all_channels_in[ilink] =  0; all_valid_in[ilink++] = muv2;

                all_channels_ref[iref] = mu; all_valid_ref[iref++] = muv; 
                all_channels_ref[iref] = 0;  all_valid_ref[iref++] = 0;
            } else {
                bool muv1 = imu < nmu, muv2 = imu+1 < nmu;
                ap_uint<64> mu1 = (muv1 ? l1pf_pattern_pack_one(mu_flat[imu  ]) : ap_uint<64>(0));
                ap_uint<64> mu2 = (muv2 ? l1pf_pattern_pack_one(mu_flat[imu+1]) : ap_uint<64>(0));
                all_channels_in[ilink] = mu1; all_valid_in[ilink++] = muv1;
                all_channels_in[ilink] =   0; all_valid_in[ilink++] = muv2;
                all_channels_in[ilink] = mu2; all_valid_in[ilink++] = muv2;

                all_channels_ref[iref] = mu1; all_valid_ref[iref++] = muv1; 
                all_channels_ref[iref] = mu2; all_valid_ref[iref++] = muv2; 
            }

            // run the unpacking
            ilink = 0; unsigned int iout = 0;
            for (int s = 0; s < NTKSECTORS; ++s) {
                unpack_track_3to2(all_channels_in[ilink+0], all_valid_in[ilink+0],
                                  all_channels_in[ilink+1], all_valid_in[ilink+1],
                                  all_channels_in[ilink+2], all_valid_in[ilink+2],
                                  all_channels_out[iout+0], all_valid_out[iout+0],
                                  all_channels_out[iout+1], all_valid_out[iout+1]);
                ilink += 3; iout += 2;
            }

            for (int s = 0; s < NCALOSECTORS; ++s) {
                for (int f = 0; f < NCALOFIBERS; ++f) {
                    unpack_hgcal_3to1(all_channels_in[ilink+0], all_valid_in[ilink+0],
                                      all_channels_in[ilink+1], all_valid_in[ilink+1],
                                      all_channels_in[ilink+2], all_valid_in[ilink+2],
                                      all_channels_out[iout+0], all_valid_out[iout+0]);
                    ilink += 3; iout += 1;
                }
            }

            unpack_mu_3to12(all_channels_in[ilink+0], all_valid_in[ilink+0],
                            all_channels_in[ilink+1], all_valid_in[ilink+1],
                            all_channels_in[ilink+2], all_valid_in[ilink+2],
                            all_channels_out[iout+0], all_valid_out[iout+0],
                            all_channels_out[iout+1], all_valid_out[iout+1]);
            ilink += 3; iout += 2;


            serPatternsIn(all_channels_in, all_valid_in);
            serPatternsOut(all_channels_out, all_valid_out);
            serPatternsRef(all_channels_ref, all_valid_ref);

            if (itest == 0 && i == 0) printf("unpacking of %d channels to %d\n", ilink, iout);
            for (int j = 0; j < iout; ++j) {
                if (all_channels_out[j] != all_channels_ref[j] || all_valid_out[j] != all_valid_ref[j]) {
                    printf("Mismatch at event %d, frame %2d, channel %3d:  %1dv%016llx vs %1dv%016llx\n", 
                                itest, i, j,
                                int(all_valid_ref[j]), all_channels_ref[j].to_uint64(),
                                int(all_valid_out[j]), all_channels_out[j].to_uint64());
                    ok = false;
                }
            }
    
            if (!ok) break;
        }
        if (!ok) break;
        printf("event %d ok\n", itest);
    } 

    fclose(fMC_tk);
    fclose(fMC_calo);
    fclose(fMC_mu);
    return ok ? 0 : 1;
}
