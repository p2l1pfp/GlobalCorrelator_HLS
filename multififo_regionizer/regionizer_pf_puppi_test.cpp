#include "firmware/regionizer.h"
#include "../utils/pattern_serializer.h"
#include "../utils/test_utils.h"

#include "../ref/pfalgo2hgc_ref.h"
#include "../puppi/linpuppi_ref.h"

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>

#define TLEN REGIONIZERNCLOCKS 

bool tk_router_ref(bool newevent, const TkObj tracks_in[NTKSECTORS][NTKFIBERS], TkObj tracks_out[NTKOUT]) ;
bool calo_router_ref(bool newevent, const HadCaloObj calo_in[NCALOSECTORS][NCALOFIBERS], HadCaloObj calo_out[NCALOOUT]) ;
bool mu_router_ref(bool newevent, const glbeta_t etaCenter, const GlbMuObj mu_in[NMUFIBERS], MuObj mu_out[NMUOUT]) ;
bool readEventTk(FILE *file, std::vector<TkObj> inputs[NTKSECTORS][NTKFIBERS], uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventCalo(FILE *file, std::vector<HadCaloObj> inputs[NCALOSECTORS][NCALOFIBERS], bool zside, uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventMu(FILE *file, std::vector<GlbMuObj> inputs[NMUFIBERS], uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventVtx(FILE *file, std::vector<std::pair<z0_t,pt_t>> & inputs, uint32_t &irun, uint32_t &ilumi, uint64_t &ievent) ;
 

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

    std::string sample = "TTbar_PU0";
    FILE *fMC_calo  = fopen(("caloDump_hgcal."+sample+".txt").c_str(), "r");
    FILE *fMC_tk  = fopen(("trackDump_hgcalPos."+sample+".txt").c_str(), "r");
    FILE *fMC_mu  = fopen(("muonDump_all."+sample+".txt").c_str(), "r");
    FILE *fMC_vtx = fopen(("vertexDump_all."+sample+".txt").c_str(), "r");
    if (!fMC_calo || !fMC_tk || !fMC_mu || !fMC_vtx) {
        printf("Couldn't open input files\n");
        return 2;
    }
    const glbeta_t etaCenter = 2*PFREGION_ETA_SIZE; // eta = +2.0

    PatternSerializer serPatternsIn("input-emp.txt"), serPatternsReg("output-emp-regionized-ref.txt"), serPatternsPf("output-emp-pf-ref.txt"), serPatternsPuppi("output-emp-puppi-ref.txt");
    assert(PACKING_NCHANN >= NTKSECTORS*NTKFIBERS + NCALOSECTORS*NCALOFIBERS + NMUFIBERS + 1);
    assert(PACKING_NCHANN >= NTKOUT + NCALOOUT + NMUOUT);
    ap_uint<64> all_channels_in[PACKING_NCHANN], all_channels_regionized[PACKING_NCHANN], all_channels_pf[PACKING_NCHANN], all_channels_puppi[PACKING_NCHANN];
    for (unsigned int i = 0; i < PACKING_NCHANN; ++i) {
        all_channels_in[i] = 0; all_channels_regionized[i] = 0; all_channels_pf[i] = 0; all_channels_puppi[i] = 0;  
    }
    serPatternsIn(all_channels_in, false); // prepend one null frame at the beginning

    int frame = 0; 

    bool ok = true; 

    z0_t pvZ0_prev = 0; // we have 1 event of delay in the reference regionizer, so we need to use the PV from 54 clocks before

    for (int itest = 0; itest < 10; ++itest) {
        std::vector<TkObj>      tk_inputs[NTKSECTORS][NTKFIBERS];
        std::vector<HadCaloObj> calo_inputs[NCALOSECTORS][NCALOFIBERS];
        std::vector<GlbMuObj>   mu_inputs[NMUFIBERS];
        std::vector<std::pair<z0_t,pt_t>> vtx_inputs;

        uint32_t run = 0, lumi = 0; uint64_t event = 0;
        if (!readEventTk(fMC_tk, tk_inputs, run, lumi, event) || 
            !readEventCalo(fMC_calo, calo_inputs, /*zside=*/true, run, lumi, event) ||
            !readEventMu(fMC_mu, mu_inputs, run, lumi, event) ||
            !readEventVtx(fMC_vtx, vtx_inputs, run, lumi, event)) break;

        z0_t vtxZ0 = vtx_inputs.empty() ? z0_t(0) : vtx_inputs.front().first;
        //if (itest == 0) printf("Vertexis at z0 = %d\n", vtxZ0.to_int());

        for (int i = 0; i < TLEN; ++i, ++frame) {
            TkObj tk_links_in[NTKSECTORS][NTKFIBERS];
            PackedTkObj tk_links64_in[NTKSECTORS][NTKFIBERS];

            unsigned int ilink = 0;

            for (int s = 0; s < NTKSECTORS; ++s) {
                for (int f = 0; f < NTKFIBERS; ++f) {
                    clear(tk_links_in[s][f]);
                    if (i < TLEN-1 && i < int(tk_inputs[s][f].size())) { // emp protocol, must leave one null frame at the end
                        tk_links_in[s][f]  = tk_inputs[s][f][i];
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
                clear(mu_links_in[f]);
                if (i < TLEN-1 && i < int(mu_inputs[f].size())) { // emp protocol, must leave one null frame at the end
                    mu_links_in[f]  = mu_inputs[f][i];
                }
                mu_links64_in[f] = l1pf_pattern_pack_one(mu_links_in[f]);
                all_channels_in[ilink++] = mu_links64_in[f];
            }

            all_channels_in[ilink++] = (i < int(vtx_inputs.size())) ? vtx_inputs[i].first : z0_t(0);

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
                linpuppi_ref(pucfg, tk_links_ref, pvZ0_prev, pfallne, outallne_nocut, outallne, outselne, itest <= 5);
                linpuppi_chs_ref(pucfg, pvZ0_prev, pfch, outallch, itest <= 5);
                l1pf_pattern_pack<NTRACK,0>(outallch, all_channels_puppi);
                l1pf_pattern_pack<NALLNEUTRALS,NTRACK>(outallne, all_channels_puppi);
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
