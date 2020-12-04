#include "firmware/regionizer.h"
#include "../utils/pattern_serializer.h"
#include "../utils/test_utils.h"
#include "firmware/obj_unpackers.h"
#include "utils/obj_packers.h"

#include "utils/tmux18_utils.h"
#include "utils/readMC.h"

#include "tdemux/tdemux_ref.h"
#include "regionizer_ref.h"
#include "../ref/pfalgo2hgc_ref.h"
#include "../puppi/linpuppi_ref.h"

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>

#define TLEN REGIONIZERNCLOCKS 


int main(int argc, char **argv) {
    pfalgo_config pfcfg(NTRACK,NCALO,NMU, NSELCALO,
                        PFALGO_DR2MAX_TK_MU, PFALGO_DR2MAX_TK_CALO,
                        PFALGO_TK_MAXINVPT_LOOSE, PFALGO_TK_MAXINVPT_TIGHT);
    linpuppi_config pucfg(NTRACK, NALLNEUTRALS, NNEUTRALS,
                          LINPUPPI_DR2MIN, LINPUPPI_DR2MAX, LINPUPPI_ptMax, LINPUPPI_dzCut,
                          LINPUPPI_etaCut, LINPUPPI_invertEta,
                          LINPUPPI_ptSlopeNe, LINPUPPI_ptSlopeNe_1, LINPUPPI_ptSlopePh, LINPUPPI_ptSlopePh_1, 
                          LINPUPPI_ptZeroNe, LINPUPPI_ptZeroNe_1, LINPUPPI_ptZeroPh, LINPUPPI_ptZeroPh_1, 
                          LINPUPPI_alphaSlope, LINPUPPI_alphaSlope_1, LINPUPPI_alphaZero, LINPUPPI_alphaZero_1, LINPUPPI_alphaCrop, LINPUPPI_alphaCrop_1, 
                          LINPUPPI_priorNe, LINPUPPI_priorNe_1, LINPUPPI_priorPh, LINPUPPI_priorPh_1,
                          LINPUPPI_ptCut, LINPUPPI_ptCut_1);

    std::string sample = "TTbar_PU200";
    FILE *fMC_calo  = fopen(("caloDump_hgcal."+sample+".txt").c_str(), "r");
    FILE *fMC_tk  = fopen(("trackDump_hgcalPos."+sample+".txt").c_str(), "r");
    FILE *fMC_mu  = fopen(("muonDump_all."+sample+".txt").c_str(), "r");
    FILE *fMC_vtx = fopen(("vertexDump_all."+sample+".txt").c_str(), "r");
    if (!fMC_calo || !fMC_tk || !fMC_mu || !fMC_vtx) {
        printf("Couldn't open input files\n");
        return 2;
    }
    const glbeta_t etaCenter = 2*PFREGION_ETA_SIZE; // eta = +2.0

    PatternSerializer serPatternsTM("input-emp.txt"), serPatternsIn("input-emp-tm6.txt");
    PatternSerializer serPatternsTDemux("input-emp-tdemux.txt"), serPatternsDecode("input-emp-decoded.txt");
    PatternSerializer serPatternsReg("output-emp-regionized-ref.txt"), serPatternsPf("output-emp-pf-ref.txt");
    PatternSerializer serPatternsPuppi("output-emp-puppi-ref.txt"), serPatternsPuppiSort("output-emp-puppisort-ref.txt");;
    assert(PACKING_NCHANN >= NTKSECTORS*3 + 3*NCALOSECTORS*NCALOFIBERS + 3 + 1);
    assert(PACKING_NCHANN >= NTKOUT + NCALOOUT + NMUOUT);
    ap_uint<64> all_channels_tmux[PACKING_NCHANN], all_channels_in[PACKING_NCHANN], all_channels_regionized[PACKING_NCHANN];
    ap_uint<64> all_channels_pf[PACKING_NCHANN], all_channels_puppi[PACKING_NCHANN], all_channels_puppisort[PACKING_NCHANN];
    ap_uint<64> all_channels_tdemux[PACKING_NCHANN], all_channels_decode[PACKING_NCHANN];
    bool all_valids_tmux[PACKING_NCHANN], all_valids_tdemux[PACKING_NCHANN], all_valids_decode[PACKING_NCHANN];
    for (unsigned int i = 0; i < PACKING_NCHANN; ++i) {
        all_channels_tmux[i] = 0; all_channels_tdemux[i] = 0;  all_channels_decode[i] = 0; 
        all_channels_in[i] = 0; all_channels_regionized[i] = 0; all_channels_pf[i] = 0; all_channels_puppi[i] = 0;  all_channels_puppisort[i] = 0;  
        all_valids_tmux[i] = 0;  all_valids_tdemux[i] = 0; all_valids_decode[i] = 0; 

    }
    serPatternsIn(all_channels_in, false); // prepend one null frame at the beginning
    serPatternsTM(all_channels_tmux, all_valids_tmux); // prepend one null frame at the beginning
    serPatternsTDemux(all_channels_tdemux, false); // prepend one null frame at the beginning
    serPatternsDecode(all_channels_decode, false); // prepend one null frame at the beginning


    // TMux encoders
    TM18LinkMultiplet<ap_uint<64>,TLEN> tk_tmuxer(NTKSECTORS), calo_tmuxer(NCALOSECTORS*NCALOFIBERS), mu_tmuxer(1);
    // TMux decoders, for testing
    TDemuxRef tk_tdemuxer[NTKSECTORS], calo_tdemuxer[NCALOSECTORS][NCALOFIBERS], mu_tdemuxer;
    // make a delay queue for the PV to realign it to the first frame
    DelayQueue pv_delayer(TLEN*2+1); // latency of the TDemuxRef (measured from first valid frame in to first valid frame out)

    int frame = 0; 
    bool ok = true; 
    z0_t pvZ0_prev = 0; // we have 1 event of delay in the reference regionizer, so we need to use the PV from 54 clocks before
    for (int itest = 0; itest < 50; ++itest) {
        std::vector<TkObj>      tk_inputs[NTKSECTORS];
        std::vector<HadCaloObj> calo_inputs[NCALOSECTORS*NCALOFIBERS];
        std::vector<GlbMuObj>   mu_inputs;
        std::vector<std::pair<z0_t,pt_t>> vtx_inputs;

        uint32_t run = 0, lumi = 0; uint64_t event = 0;
        if (!readEventTk(fMC_tk, tk_inputs, run, lumi, event) || 
            !readEventCalo(fMC_calo, calo_inputs, /*zside=*/true, run, lumi, event) ||
            !readEventMu(fMC_mu, mu_inputs, run, lumi, event) ||
            !readEventVtx(fMC_vtx, vtx_inputs, run, lumi, event)) {
                printf("Reached end of input file.\n");
            break;
        }

        // enqueue frames (outside of the frame loop, since it takes 3*TLEN and not TLEN)
        tk_tmuxer.push_links(itest, tk_inputs, pack_tracks);
        calo_tmuxer.push_links(itest, calo_inputs, pack_hgcal);
        mu_tmuxer.push_link(itest, mu_inputs, pack_muons);

        z0_t vtxZ0 = vtx_inputs.empty() ? z0_t(0) : vtx_inputs.front().first;
        //if (itest == 0) printf("Vertexis at z0 = %d\n", vtxZ0.to_int());

        for (int i = 0; i < TLEN; ++i, ++frame) {
            TkObj tk_links_in[NTKSECTORS][NTKFIBERS];
            PackedTkObj tk_links64_in[NTKSECTORS][NTKFIBERS];

            unsigned int ilink = 0;

            for (int s = 0; s < NTKSECTORS; ++s) {
                for (int f = 0; f < NTKFIBERS; ++f) {
                    int itk = 2*i+f;
                    clear(tk_links_in[s][f]);
                    if (itk < TLEN-1 && itk < int(tk_inputs[s].size())) { // emp protocol, must leave one null frame at the end
                        tk_links_in[s][f]  = tk_inputs[s][itk];
                    }
                    tk_links64_in[s][f] = l1pf_pattern_pack_one(tk_links_in[s][f]);
                    all_channels_in[ilink++] = tk_links64_in[s][f];
                }
            }

            HadCaloObj    calo_links_in[NCALOSECTORS][NCALOFIBERS];
            PackedCaloObj calo_links64_in[NCALOSECTORS][NCALOFIBERS];
            for (int s = 0; s < NCALOSECTORS; ++s) {
                for (int f = 0; f < NCALOFIBERS; ++f) {
                    clear(calo_links_in[s][f]);
                    if (i < TLEN-1 && i < int(calo_inputs[s*NCALOFIBERS+f].size())) { // emp protocol, must leave one null frame at the end
                        calo_links_in[s][f]  = calo_inputs[s*NCALOFIBERS+f][i];
                    }
                    calo_links64_in[s][f] = l1pf_pattern_pack_one(calo_links_in[s][f]);
                    all_channels_in[ilink++] = calo_links64_in[s][f];
                }
            }

            GlbMuObj    mu_links_in[NMUFIBERS];
            PackedMuObj mu_links64_in[NMUFIBERS];
            for (int f = 0; f < NMUFIBERS; ++f) {
                int imu = 3*i/2+f;
                clear(mu_links_in[f]);
                if ((f == 0 || i%2 == 1) && imu < TLEN-1 && imu < int(mu_inputs.size())) { // emp protocol, must leave one null frame at the end
                    mu_links_in[f]  = mu_inputs[imu];
                }
                mu_links64_in[f] = l1pf_pattern_pack_one(mu_links_in[f]);
                all_channels_in[ilink++] = mu_links64_in[f];
            }

            all_channels_in[ilink++] = (i < int(vtx_inputs.size())) ? vtx_inputs[i].first : z0_t(0);

            // pop out frames from the tmuxer for printing. for each sector, we put the 3 links for 3 set of events next to each other
            unsigned int calo_offs = NTKSECTORS*3, mu_offs = calo_offs + NCALOSECTORS * NCALOFIBERS * 3, vtx_offs = mu_offs + 3;
            tk_tmuxer.pop_frame(all_channels_tmux, all_valids_tmux); 
            calo_tmuxer.pop_frame(all_channels_tmux, all_valids_tmux, calo_offs);
            mu_tmuxer.pop_frame(all_channels_tmux, all_valids_tmux, mu_offs);
            // the vertex is not TMUXed so we just add it at the end
            all_channels_tmux[vtx_offs] = (i < int(vtx_inputs.size())) ? vtx_inputs[i].first : z0_t(0);
            all_valids_tmux[vtx_offs] = (i < TLEN-1);

            // now let's run the time demultiplexer
            bool newEvt = (i == 0 && itest == 0);
            for (int s = 0; s < NTKSECTORS; ++s) {
                tk_tdemuxer[s]( newEvt, &all_channels_tmux  [3*s], &all_valids_tmux  [3*s],
                                        &all_channels_tdemux[3*s], &all_valids_tdemux[3*s]);
            }
            for (int s = 0; s < NCALOSECTORS; ++s) {
                for (int f = 0; f < NCALOFIBERS; ++f) {
                    ilink = calo_offs + 3*(s*NCALOFIBERS + f);
                    calo_tdemuxer[s][f]( newEvt, &all_channels_tmux  [ilink], &all_valids_tmux  [ilink],
                                                 &all_channels_tdemux[ilink], &all_valids_tdemux[ilink]);
                }
            }
            mu_tdemuxer(newEvt, &all_channels_tmux  [mu_offs], &all_valids_tmux  [mu_offs],
                                &all_channels_tdemux[mu_offs], &all_valids_tdemux[mu_offs]);
            // note: the PV is delayed to realign it to the other demuxed channels
            pv_delayer(all_channels_tmux  [vtx_offs], all_valids_tmux  [vtx_offs],
                       all_channels_tdemux[vtx_offs], all_valids_tdemux[vtx_offs]);

            // and now we unpack to 64 bit format
            ilink = 0; unsigned int iout = 0;
            for (int s = 0; s < NTKSECTORS; ++s) {
                unpack_track_3to2(all_channels_tdemux[ilink+0], all_valids_tdemux[ilink+0],
                                  all_channels_tdemux[ilink+1], all_valids_tdemux[ilink+1],
                                  all_channels_tdemux[ilink+2], all_valids_tdemux[ilink+2],
                                  all_channels_decode[iout+0], all_valids_decode[iout+0],
                                  all_channels_decode[iout+1], all_valids_decode[iout+1]);
                ilink += 3; iout += 2;
            }
            for (int s = 0; s < NCALOSECTORS; ++s) {
                for (int f = 0; f < NCALOFIBERS; ++f) {
                    unpack_hgcal_3to1(all_channels_tdemux[ilink+0], all_valids_tdemux[ilink+0],
                                      all_channels_tdemux[ilink+1], all_valids_tdemux[ilink+1],
                                      all_channels_tdemux[ilink+2], all_valids_tdemux[ilink+2],
                                      all_channels_decode[iout+0], all_valids_decode[iout+0]);
                    ilink += 3; iout += 1;
                }
            }
            unpack_mu_3to12(all_channels_tdemux[ilink+0], all_valids_tdemux[ilink+0],
                            all_channels_tdemux[ilink+1], all_valids_tdemux[ilink+1],
                            all_channels_tdemux[ilink+2], all_valids_tdemux[ilink+2],
                            all_channels_decode[iout+0], all_valids_decode[iout+0],
                            all_channels_decode[iout+1], all_valids_decode[iout+1]); 
            if (all_valids_decode[iout+0]) all_valids_decode[iout+1] = 1; // for our purposes, mark both valid if the first is valid
            ilink += 3; iout += 2;
            // the vertex is trivial
            all_channels_decode[iout] = all_channels_tdemux[ilink];
            all_valids_decode  [iout] = all_valids_tdemux  [ilink];
            ilink++; iout++;
            // done unpacking


            TkObj        tk_links_ref[NTKOUT];
            HadCaloObj calo_links_ref[NCALOOUT];
            MuObj        mu_links_ref[NMUOUT]; 

            bool newevt_ref = (i == 0);
            bool   tk_ref_good =   tk_router_ref(i == 0,   tk_links_in, tk_links_ref);
            bool calo_ref_good = calo_router_ref(i == 0, calo_links_in, calo_links_ref);
            bool   mu_ref_good =   mu_router_ref(i == 0, etaCenter, mu_links_in, mu_links_ref);

            l1pf_pattern_pack<NTKOUT,0>(tk_links_ref, all_channels_regionized);
            l1pf_pattern_pack<NCALOOUT,NTKOUT>(calo_links_ref, all_channels_regionized);
            l1pf_pattern_pack<NMUOUT,NTKOUT+NCALOOUT>(mu_links_ref, all_channels_regionized);

            if ((itest > 0) && (i % PFLOWII == 0) && (i/PFLOWII < NPFREGIONS)) {
                ////  ok we can run PF and puppi
                // PF objects
                PFChargedObj pfch[NTRACK], pfmu[NMU]; PFNeutralObj pfallne[NALLNEUTRALS];
                assert(NTKOUT == NTRACK && NCALOOUT == NCALO && NMUOUT == NMU);
                if (itest <= 5) printf("Will run PF event %d, region %d\n", itest-1, i/PFLOWII);
                pfalgo2hgc_ref_set_debug(itest <= 5);
                pfalgo2hgc_ref(pfcfg, calo_links_ref, tk_links_ref, mu_links_ref, pfch, pfallne, pfmu); 
                pfalgo2hgc_pack_out(pfch, pfallne, pfmu, all_channels_pf);
                // Puppi objects
                if (itest <= 5) printf("Will run Puppi with z0 = %d in event %d, region %d\n", pvZ0_prev.to_int(), itest-1, i/PFLOWII);
                PuppiObj outallch[NTRACK];
                PuppiObj outallne_nocut[NALLNEUTRALS], outallne[NALLNEUTRALS], outselne[NNEUTRALS]; 
                PuppiObj outpresort[NTRACK+NALLNEUTRALS];
                linpuppi_ref(pucfg, tk_links_ref, pvZ0_prev, pfallne, outallne_nocut, outallne, outselne, itest <= 1);
                linpuppi_chs_ref(pucfg, pvZ0_prev, pfch, outallch, itest <= 1);
                std::copy(outallch, outallch+NTRACK, outpresort);
                std::copy(outallne, outallne+NALLNEUTRALS, outpresort+NTRACK);
                PuppiObj outsorted[NPUPPIFINALSORTED];
                puppisort_and_crop_ref(NTRACK+NALLNEUTRALS, NPUPPIFINALSORTED, outpresort, outsorted);
                l1pf_pattern_pack<NTRACK+NALLNEUTRALS,0>(outpresort, all_channels_puppi);
                l1pf_pattern_pack<NPUPPIFINALSORTED,0>(outsorted, all_channels_puppisort);
            }

            serPatternsTM(all_channels_tmux, all_valids_tmux);
            if (frame >= 108) {
                serPatternsTDemux(all_channels_tdemux, all_valids_tdemux);
                serPatternsDecode(all_channels_decode, all_valids_decode);
            }
            serPatternsIn(all_channels_in, (i < TLEN-1));
            serPatternsReg(all_channels_regionized);
            serPatternsPf(all_channels_pf);
            serPatternsPuppi(all_channels_puppi);
            serPatternsPuppiSort(all_channels_puppisort);

            if (i == TLEN-1) pvZ0_prev = vtx_inputs.empty() ? z0_t(0) : vtx_inputs.front().first;
        }
    } 

    fclose(fMC_tk);
    fclose(fMC_calo);
    fclose(fMC_mu);
    fclose(fMC_vtx);
    return ok ? 0 : 1;
}
