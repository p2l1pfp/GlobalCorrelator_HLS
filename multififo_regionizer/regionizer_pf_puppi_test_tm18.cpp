#include "firmware/regionizer.h"
#include "../utils/pattern_serializer.h"
#include "../utils/test_utils.h"
#include "firmware/obj_unpackers.h"
#include "tdemux/tdemux_ref.h"

#include "../ref/pfalgo2hgc_ref.h"
#include "../puppi/linpuppi_ref.h"

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <queue>
#include <string>

#define TLEN REGIONIZERNCLOCKS 

bool tk_router_ref(bool newevent, const TkObj tracks_in[NTKSECTORS][NTKFIBERS], TkObj tracks_out[NTKOUT]) ;
bool calo_router_ref(bool newevent, const HadCaloObj calo_in[NCALOSECTORS][NCALOFIBERS], HadCaloObj calo_out[NCALOOUT]) ;
bool mu_router_ref(bool newevent, const glbeta_t etaCenter, const GlbMuObj mu_in[NMUFIBERS], MuObj mu_out[NMUOUT]) ;
bool readEventTkTM18(FILE *file, std::vector<TkObj> inputs[NTKSECTORS], uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventCalo(FILE *file, std::vector<HadCaloObj> inputs[NCALOSECTORS][NCALOFIBERS], bool zside, uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventMuTM18(FILE *file, std::vector<GlbMuObj> inputs, uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventVtx(FILE *file, std::vector<std::pair<z0_t,pt_t>> & inputs, uint32_t &irun, uint32_t &ilumi, uint64_t &ievent) ;

template<typename T>
struct TM18LinkTriplet {
    std::queue<std::pair<T,bool>> links[3];
    
    TM18LinkTriplet() { 
        for (int i = 0; i <  TLEN; ++i)  links[1].emplace(T(0), false);
        for (int i = 0; i < 2*TLEN; ++i) links[2].emplace(T(0), false);
    };


    template<typename C>
    void push_event(unsigned int iev, const C & objs) {
        auto & q = links[iev % 3];
        unsigned int n = std::min<unsigned int>(3*TLEN-3, objs.size()); // let's leave 3 empty frames at the end, so they get reassebled as 1 row of nulls
        //printf("Writing %u objects on link %u for iev %u\n", n, iev%3, iev);
        for (unsigned int i = 0; i < n; ++i) {
            q.emplace(objs[i], true);
        }

        for (unsigned int i = n; i < 3*TLEN; ++i) {
            q.emplace(T(0), i < 3*TLEN-3);
        }
    }
    template<typename TC, typename BC>
    void pop_frame(TC & values, BC & valids, unsigned int start=0, unsigned int stride=1) {
        for (int i = 0; i < 3; ++i) {
            if (links[i].empty()) { 
                printf("ERROR: link %d is empty (start = %u, stride = %u)\n", i, start, stride); fflush(stdout); 
                continue;
            }
            assert(!links[i].empty());
            auto obj = links[i].front();
            values[start+i*stride] = obj.first;
            valids[start+i*stride] = obj.second;
            links[i].pop();
        }
    }
};

#ifdef TRIVIAL_ENCODING_64
template<typename T>
std::vector<ap_uint<64>> encode_objs(const std::vector<T> & objs) {
    std::vector<ap_uint<64>> ret;
    for (unsigned int i = 0, n = objs.size(); i < n; ++i) {
        ret.push_back(l1pf_pattern_pack_one(objs[i]));
    }
    return ret;
}
#else
std::vector<ap_uint<64>> encode_objs(const std::vector<TkObj> & tracks) {
    std::vector<ap_uint<64>> ret;
    for (unsigned int i = 0, n = tracks.size(); i < n; ++i) {
        // simulate 96 bit objects
        ap_uint<96> packedtk = l1pf_pattern_pack_one(tracks[i]);
        if (i % 2 == 0) {
            ret.emplace_back(packedtk(95,32));
            ret.emplace_back((packedtk(31,0), ap_uint<32>(0)));
        } else {
            ret.back()(31,0) = packedtk(95,64);
            ret.emplace_back(packedtk(63,0));
        }
    }
    return ret;
}
std::vector<ap_uint<64>> encode_objs(const std::vector<HadCaloObj> & calo) {
    std::vector<ap_uint<64>> ret;
    for (unsigned int i = 0, n = calo.size(); i < n; ++i) {
        // simulate 128 bit objects on a 16 G link
        ap_uint<128> packed = l1pf_pattern_pack_one(calo[i]);
        ret.emplace_back(packed(127,64));
        ret.emplace_back(packed( 63, 0));
        ret.push_back(0); // third zero frame, that will have the strobe bit off
    }
    return ret;

}
std::vector<ap_uint<64>> encode_objs(const std::vector<GlbMuObj> & mu) {
    std::vector<ap_uint<64>> ret;
    for (unsigned int i = 0, n = mu.size(); i < n; ++i) {
        // simulate 128 bit objects on a 16 G link
        ap_uint<128> packed = l1pf_pattern_pack_one(mu[i]);
        ret.emplace_back(packed(127,64));
        ret.emplace_back(packed( 63, 0));
    }
    return ret;
}
#endif

class DelayQueue {
    public:
        DelayQueue(unsigned int n) : n_(n), data_(n, 0), ptr_(0) {}
        ap_uint<65> operator()(const ap_uint<65> & in) {
            ap_uint<65> ret = data_[ptr_];
            data_[ptr_] = in;
            ptr_++; 
            if (ptr_ == n_) ptr_ = 0;
            return ret;
        }
        void operator()(const ap_uint<64> & in,  const bool & in_valid,
                              ap_uint<64> & out,       bool & out_valid) {
            ap_uint<65> in65  = in; in65[64] = in_valid;
            ap_uint<65> out65 = (*this)(in65);
            out = out65(63,0); out_valid = out65[64];
        }
    private:
        unsigned int n_, ptr_;
        std::vector<ap_uint<65>> data_;
};



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
    PatternSerializer serPatternsReg("output-emp-regionized-ref.txt"), serPatternsPf("output-emp-pf-ref.txt"), serPatternsPuppi("output-emp-puppi-ref.txt");
    assert(PACKING_NCHANN >= NTKSECTORS*3 + 3*NCALOSECTORS*NCALOFIBERS + 3 + 1);
    assert(PACKING_NCHANN >= NTKOUT + NCALOOUT + NMUOUT);
    ap_uint<64> all_channels_tmux[PACKING_NCHANN], all_channels_in[PACKING_NCHANN], all_channels_regionized[PACKING_NCHANN], all_channels_pf[PACKING_NCHANN], all_channels_puppi[PACKING_NCHANN];
    ap_uint<64> all_channels_tdemux[PACKING_NCHANN], all_channels_decode[PACKING_NCHANN];
    bool all_valids_tmux[PACKING_NCHANN], all_valids_tdemux[PACKING_NCHANN], all_valids_decode[PACKING_NCHANN];
    for (unsigned int i = 0; i < PACKING_NCHANN; ++i) {
        all_channels_tmux[i] = 0; all_channels_tdemux[i] = 0;  all_channels_decode[i] = 0; 
        all_channels_in[i] = 0; all_channels_regionized[i] = 0; all_channels_pf[i] = 0; all_channels_puppi[i] = 0;  
        all_valids_tmux[i] = 0;  all_valids_tdemux[i] = 0; all_valids_decode[i] = 0; 

    }
    serPatternsIn(all_channels_in, false); // prepend one null frame at the beginning
    serPatternsTM(all_channels_tmux, all_valids_tmux); // prepend one null frame at the beginning
    serPatternsTDemux(all_channels_tdemux, false); // prepend one null frame at the beginning
    serPatternsDecode(all_channels_decode, false); // prepend one null frame at the beginning


    // TMux encoders
    TM18LinkTriplet<ap_uint<64>> tk_tmuxer[NTKSECTORS], calo_tmuxer[NCALOSECTORS][NCALOFIBERS], mu_tmuxer;
    // TMux decoders, for testing
    TDemuxRef tk_tdemuxer[NTKSECTORS], calo_tdemuxer[NCALOSECTORS][NCALOFIBERS], mu_tdemuxer;
    // make a delay queue for the PV to realign it to the first frame
    DelayQueue pv_delayer(TLEN*2+1); // latency of the TDemuxRef (measured from first valid frame in to first valid frame out)

    int frame = 0; 
    bool ok = true; 
    z0_t pvZ0_prev = 0; // we have 1 event of delay in the reference regionizer, so we need to use the PV from 54 clocks before
    for (int itest = 0; itest < 10; ++itest) {
        std::vector<TkObj>      tk_inputs[NTKSECTORS];
        std::vector<HadCaloObj> calo_inputs[NCALOSECTORS][NCALOFIBERS];
        std::vector<GlbMuObj>   mu_inputs;
        std::vector<std::pair<z0_t,pt_t>> vtx_inputs;

        uint32_t run = 0, lumi = 0; uint64_t event = 0;
        if (!readEventTkTM18(fMC_tk, tk_inputs, run, lumi, event) || 
            !readEventCalo(fMC_calo, calo_inputs, /*zside=*/true, run, lumi, event) ||
            !readEventMuTM18(fMC_mu, mu_inputs, run, lumi, event) ||
            !readEventVtx(fMC_vtx, vtx_inputs, run, lumi, event)) break;

        // enqueue frames (outside of the frame loop, since it takes 3*TLEN and not TLEN)
        for (int s = 0; s < NTKSECTORS; ++s) {
            tk_tmuxer[s].push_event(itest, encode_objs(tk_inputs[s]));
        }
        for (int s = 0; s < NCALOSECTORS; ++s) {
            for (int f = 0; f < NCALOFIBERS; ++f) {
                calo_tmuxer[s][f].push_event(itest, encode_objs(calo_inputs[s][f]));
            }
        }
        mu_tmuxer.push_event(itest, encode_objs(mu_inputs));

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
                    if (i < TLEN-1 && i < int(calo_inputs[s][f].size())) { // emp protocol, must leave one null frame at the end
                        calo_links_in[s][f]  = calo_inputs[s][f][i];
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
                if ((f == 0 || i%2 == 0) && imu < TLEN-1 && imu < int(mu_inputs.size())) { // emp protocol, must leave one null frame at the end
                    mu_links_in[f]  = mu_inputs[imu];
                }
                mu_links64_in[f] = l1pf_pattern_pack_one(mu_links_in[f]);
                all_channels_in[ilink++] = mu_links64_in[f];
            }

            all_channels_in[ilink++] = (i < int(vtx_inputs.size())) ? vtx_inputs[i].first : z0_t(0);

            // pop out frames from the tmuxer for printing. for each sector, we put the 3 links for 3 set of events next to each other
            unsigned int calo_offs = NTKSECTORS*3, mu_offs = calo_offs + NCALOSECTORS * NCALOFIBERS * 3, vtx_offs = mu_offs + 3;
            for (int s = 0; s < NTKSECTORS; ++s) {
                tk_tmuxer[s].pop_frame(all_channels_tmux, all_valids_tmux, 3*s); //s, NTKSECTORS); 
            }
            for (int s = 0; s < NCALOSECTORS; ++s) {
                for (int f = 0; f < NCALOFIBERS; ++f) {
                    calo_tmuxer[s][f].pop_frame(all_channels_tmux, all_valids_tmux, calo_offs+3*(s*NCALOFIBERS+f)); // s*NCALOFIBERS+f, NCALOSECTORS*NCALOFIBERS);
                }
            }
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
                PFChargedObj outallch[NTRACK];
                PFNeutralObj outallne_nocut[NALLNEUTRALS], outallne[NALLNEUTRALS], outselne[NNEUTRALS]; 
                linpuppi_ref(pucfg, tk_links_ref, pvZ0_prev, pfallne, outallne_nocut, outallne, outselne, itest <= 1);
                linpuppi_chs_ref(pucfg, pvZ0_prev, pfch, outallch, itest <= 1);
                l1pf_pattern_pack<NTRACK,0>(outallch, all_channels_puppi);
                l1pf_pattern_pack<NALLNEUTRALS,NTRACK>(outallne, all_channels_puppi);
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

            if (i == TLEN-1) pvZ0_prev = vtx_inputs.empty() ? z0_t(0) : vtx_inputs.front().first;
        }
    } 

    fclose(fMC_tk);
    fclose(fMC_calo);
    fclose(fMC_mu);
    fclose(fMC_vtx);
    return ok ? 0 : 1;
}
