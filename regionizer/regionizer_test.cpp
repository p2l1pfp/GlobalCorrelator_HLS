#include <cstdio>
#include "firmware/regionizer.h"
#include "../firmware/mp7pf_encoding.h"
#include "../utils/random_inputs.h"
#include "../utils/DiscretePFInputs_IO.h"
#include "../utils/pattern_serializer.h"
#include "../utils/test_utils.h"

#define NTEST 50

template<unsigned int N, typename T>
bool fill_stream(hls::stream<T> & stream, T data[N], int step=1, int offset=0, const char *name="unknown", int isector=-999) {
    if (!stream.empty()) { printf("ERROR: %s stream for sector %d is not empty\n", name, isector); return false; }
    for (unsigned int i = offset; i < N; i += step) { stream.write(data[i]); }
    return true;
}

void dump_o(FILE *f, const HadCaloObj & obj) { 
    // NOTE: no leading + on positive numbers, VHDL doesn't like it
    fprintf(f, "       %4d % 4d % 4d %3d %1d", int(obj.hwPt), int(obj.hwEta), int(obj.hwPhi), int(obj.hwEmPt), obj.hwIsEM); 
}
void dump_z(FILE *f, const HadCaloObj &) { // the second argument is only to resolve overloading
    HadCaloObj dummy; clear(dummy); dump_o(f,dummy);
}
void dump_o(FILE *f, const EmCaloObj & obj) { 
    // NOTE: no leading + on positive numbers, VHDL doesn't like it
    fprintf(f, "       %4d % 4d % 4d %3d", int(obj.hwPt), int(obj.hwEta), int(obj.hwPhi), int(obj.hwPtErr)); 
}
void dump_z(FILE *f, const EmCaloObj &) { // the second argument is only to resolve overloading
    EmCaloObj dummy; clear(dummy); dump_o(f,dummy);
}

void dump_o(FILE *f, const TkObj & obj) { 
    // NOTE: no leading + on positive numbers, VHDL doesn't like it
    fprintf(f, "       %4d % 4d % 4d %4d %4d", int(obj.hwPt), int(obj.hwEta), int(obj.hwPhi), int(obj.hwPtErr), int(obj.hwZ0)); 
}
void dump_z(FILE *f, const TkObj &) { // the second argument is only to resolve overloading
    TkObj dummy; clear(dummy); dump_o(f,dummy);
}
void dump_o(FILE *f, const MuObj & obj) { 
    // NOTE: no leading + on positive numbers, VHDL doesn't like it
    fprintf(f, "       %4d % 4d % 4d %4d", int(obj.hwPt), int(obj.hwEta), int(obj.hwPhi), int(obj.hwPtErr)); 
}
void dump_z(FILE *f, const MuObj &) { // the second argument is only to resolve overloading
    MuObj dummy; clear(dummy); dump_o(f,dummy);
}

int main() {

    // input format: could be random or coming from simulation
    //RandomPFInputs inputs(37); // 37 is a good random number
    DiscretePFInputs inputs("barrel_sectors_1x12_TTbar_PU140.dump");
    HumanReadablePatternSerializer debug("-"); // this will print on stdout, we'll use it for errors

    printf(" --- Configuration --- \n");
    printf(" Sectors: %d \n", N_IN_SECTORS);
    printf("    max N(Calo), input:     %2d \n", NCALO_PER_SECTOR);
    printf("    max N(Calo), eta slice: %2d \n", NCALO_PER_SECTOR_PER_ETA);
    printf("    max N(EmCalo), input:     %2d \n", NEMCALO_PER_SECTOR);
    printf("    max N(EmCalo), eta slice: %2d \n", NEMCALO_PER_SECTOR_PER_ETA);
    printf("    max N(Track), input:     %2d \n", NTRACK_PER_SECTOR);
    printf("    max N(Track), eta slice: %2d \n", NTRACK_PER_SECTOR_PER_ETA);
    printf(" Regions: %d (%d x %d )\n", N_OUT_REGIONS, N_OUT_REGIONS_ETA, N_OUT_REGIONS_PHI);
    printf("    max N(Calo): %2d \n", NCALO);
    printf("    max N(EmCalo): %2d \n", NEMCALO);
    printf("    max N(Track): %2d \n", NTRACK);

    HadCaloObj calo_in[N_IN_SECTORS][NCALO_PER_SECTOR];
    hls::stream<HadCaloObj> calo_fibers[N_IN_SECTORS];
    hls::stream<HadCaloObj> calo_fibers_ref[N_IN_SECTORS];
    HadCaloObj calo_regions[N_OUT_REGIONS][NCALO]; 
    HadCaloObj calo_regions_ref[N_OUT_REGIONS][NCALO]; 

    EmCaloObj emcalo_in[N_IN_SECTORS][NEMCALO_PER_SECTOR];
    hls::stream<EmCaloObj> emcalo_fibers[N_IN_SECTORS];
    hls::stream<EmCaloObj> emcalo_fibers_ref[N_IN_SECTORS];
    EmCaloObj emcalo_regions[N_OUT_REGIONS][NEMCALO]; 
    EmCaloObj emcalo_regions_ref[N_OUT_REGIONS][NEMCALO]; 

    TkObj track_in[N_IN_SECTORS][NTRACK_PER_SECTOR];
    hls::stream<TkObj> track_fibers[2*N_IN_SECTORS]; // two fibers per sector
    hls::stream<TkObj> track_fibers_ref[2*N_IN_SECTORS];
    TkObj track_regions[N_OUT_REGIONS][NTRACK]; 
    TkObj track_regions_ref[N_OUT_REGIONS][NTRACK]; 

    MuObj mu_in_cmssw[N_IN_SECTORS][NMU]; // 12-fold 
    MuObj mu_in[N_MUON_SECTORS][NMU]; // 4-fold
    hls::stream<MuObj> mu_fibers[N_MUON_SECTORS]; 
    hls::stream<MuObj> mu_fibers_ref[N_MUON_SECTORS];
    MuObj mu_regions[N_OUT_REGIONS][NMU]; 
    MuObj mu_regions_ref[N_OUT_REGIONS][NMU]; 


    FILE *f_in  = fopen("dump_in.txt","w");
    FILE *f_out = fopen("dump_out.txt","w");
    int frame_in = 0, frame_out = 0;

#ifdef MP7
    HadCaloObj calo_in_transposed[NCALO_PER_SECTOR][N_IN_SECTORS];
    EmCaloObj emcalo_in_transposed[NEMCALO_PER_SECTOR][N_IN_SECTORS];
    TkObj track_in_transposed[NTRACK_PER_SECTOR/2][2*N_IN_SECTORS];
    MuObj mu_in_transposed[NMU][N_MUON_SECTORS];
    MP7PatternSerializer serMP7_in( "mp7_input.txt",2,1);  
    MP7PatternSerializer serMP7_in_calo( "mp7_input_calo.txt",2,1);  
    MP7PatternSerializer serMP7_in_track( "mp7_input_track.txt",2,1);  
    MP7PatternSerializer serMP7_in_mu( "mp7_input_mu.txt",2,1);  
    MP7PatternSerializer serMP7_out("mp7_output.txt",2,1); 
    MP7PatternSerializer serMP7_out_calo("mp7_output_calo.txt",2,1); 
    MP7PatternSerializer serMP7_out_track("mp7_output_track.txt",2,1); 
    MP7PatternSerializer serMP7_out_mu("mp7_output_mu.txt",2,1); 
    MP7PatternSerializer serMP7_out_calo_nomux("mp7_output_calo_nomux.txt",1,0); 
    MP7PatternSerializer serMP7_out_track_nomux("mp7_output_track_nomux.txt",1,0); 
    MP7DataWord mp7_in[MP7_NCHANN];
    MP7DataWord mp7_out[MP7_NCHANN];
#endif

    // -----------------------------------------
    // run multiple tests
    for (int test = 1; test <= NTEST; ++test) {
        // read the event
        if (!inputs.nextEvent()) break;
        if (inputs.event().regions.size() != N_IN_SECTORS) { printf("ERROR: Mismatching number of input regions: %lu\n", inputs.event().regions.size()); return 2; }
        //fill in the streams
        for (int is = 0; is < N_IN_SECTORS; ++is) {
            const Region & r = inputs.event().regions[is];
            // CALO
            dpf2fw::convert<NCALO_PER_SECTOR>(r.calo, calo_in[is]); 
            for (unsigned int i = 0; i < NCALO_PER_SECTOR; ++i) assert(calo_in[is][i].hwPt >= 0);
            if (!fill_stream<NCALO_PER_SECTOR>(calo_fibers[is], calo_in[is], 1, 0, "calo stream", is)) return 3;
            if (!fill_stream<NCALO_PER_SECTOR>(calo_fibers_ref[is], calo_in[is], 1, 0, "calo ref stream", is)) return 3;
            // EMCALO
            dpf2fw::convert<NEMCALO_PER_SECTOR>(r.emcalo, emcalo_in[is]); 
            if (!fill_stream<NEMCALO_PER_SECTOR>(emcalo_fibers[is], emcalo_in[is], 1, 0, "emcalo stream", is)) return 3;
            if (!fill_stream<NEMCALO_PER_SECTOR>(emcalo_fibers_ref[is], emcalo_in[is], 1, 0, "emcalo ref stream", is)) return 3;
            // TRACK
            dpf2fw::convert<NTRACK_PER_SECTOR>(r.track, track_in[is]); 
            for (unsigned int i = 0; i < 2; ++i) {
                if (!fill_stream<NTRACK_PER_SECTOR>(track_fibers[2*is+i],     track_in[is], 2, i, "track stream ",    2*is+i)) return 3;
                if (!fill_stream<NTRACK_PER_SECTOR>(track_fibers_ref[2*is+i], track_in[is], 2, i, "track ref stream", 2*is+i)) return 3;
            }
            // MUON (12-fold CMSSW input)
            dpf2fw::convert<NMU>(r.muon, mu_in_cmssw[is]); 
            //for (int i = 0; i < NMU; ++i) if (mu_in_cmssw[is][i].hwPt > 0) printf("Muon 12-fold %2d/%d of hwPt %6d local phi %+6d  global phi %+8d  sector %+6.4f\n", is, i,
            //    int(mu_in_cmssw[is][i].hwPt), int(mu_in_cmssw[is][i].hwPhi), int(mu_in_cmssw[is][i].hwPhi) + is * _PHI_PIO6 + _PHI_PIO6/2, (int(mu_in_cmssw[is][i].hwPhi) + is * _PHI_PIO6 + _PHI_PIO6/2)/float(_PHI_PIO6));
        }
        // MUONS need a dedicate handling to fake a 4-fold readout instead of a 12-fold
        merge_muon_in(mu_in_cmssw, mu_in); // 12-fold to 4-fold
        for (int is = 0; is < N_MUON_SECTORS; ++is) {
            //for (int i = 0; i < NMU; ++i) if (mu_in[is][i].hwPt > 0) printf("Muon  4-fold %2d/%d of hwPt %6d local phi %+6d  global phi %+8d\n", is, i,
            //    int(mu_in[is][i].hwPt), int(mu_in[is][i].hwPhi), int(mu_in[is][i].hwPhi) + is * 3*_PHI_PIO6 + 3*_PHI_PIO6/2);
            if (!fill_stream<NMU>(mu_fibers[is], mu_in[is], 1, 0, "mu stream", is)) return 3;
            if (!fill_stream<NMU>(mu_fibers_ref[is], mu_in[is], 1, 0, "mu ref stream", is)) return 3;
        }
        // dump inputs
        for (unsigned int ic = 0; ic < N_CLOCKS; ++ic) {
            // takes 2 clocks to send one input; so we just duplicate the lines for now
            int iobj = ic/2; bool send = (ic % 2 == 0);
            fprintf(f_in,"Frame %04d : %2d %2d", ++frame_in, test, iobj);
            for (int is = 0; is < N_IN_SECTORS; ++is) {
                if (iobj < NCALO_PER_SECTOR && send) dump_o(f_in, calo_in[is][iobj]);
                else                                 dump_z(f_in, calo_in[is][0]);
            }
            for (int is = 0; is < N_IN_SECTORS; ++is) {
                if (iobj < NEMCALO_PER_SECTOR && send) dump_o(f_in, emcalo_in[is][iobj]);
                else                                   dump_z(f_in, emcalo_in[is][0]);
            }
            iobj = ic;
            for (int is = 0; is < N_IN_SECTORS; ++is) {
                if (iobj+0 < NTRACK_PER_SECTOR/2 && send) dump_o(f_in, track_in[is][iobj+0]);
                else                                      dump_z(f_in, track_in[is][0]);
                if (iobj+1 < NTRACK_PER_SECTOR/2 && send) dump_o(f_in, track_in[is][iobj+1]);
                else                                      dump_z(f_in, track_in[is][0]);
            }
            iobj = ic/2; 
            for (int is = 0; is < N_MUON_SECTORS; ++is) {
                if (iobj < NMU && send) dump_o(f_in, mu_in[is][iobj]);
                else                    dump_z(f_in, mu_in[is][0]);
            }
            fprintf(f_in,"\n");
        }
#ifdef MP7
        for (int is = 0; is < N_IN_SECTORS; ++is) { 
            for (int io = 0; io < NCALO_PER_SECTOR; ++io) calo_in_transposed[io][is] = calo_in[is][io];
            for (int io = 0; io < NEMCALO_PER_SECTOR; ++io) emcalo_in_transposed[io][is] = emcalo_in[is][io];
            for (int io = 0; io < NTRACK_PER_SECTOR; ++io) track_in_transposed[io/2][2*is+(io%2)] = track_in[is][io];
        } 
        for (int is = 0; is < N_MUON_SECTORS; ++is) {
            for (int io = 0; io < NMU; ++io) mu_in_transposed[io][is] = mu_in[is][io];
        }
        
        for (unsigned int ic = 0; ic < N_CLOCKS/2; ++ic) {
            for (unsigned int i = 0; i < MP7_NCHANN; ++i) mp7_in[i] = 0; // clear
            if (ic < NCALO_PER_SECTOR) mp7_pack<N_IN_SECTORS,0>(calo_in_transposed[ic], mp7_in);
            if (ic < NEMCALO_PER_SECTOR) mp7_pack<N_IN_SECTORS,2*N_IN_SECTORS>(emcalo_in_transposed[ic], mp7_in);
            if (ic < NTRACK_PER_SECTOR/2) mp7_pack<2*N_IN_SECTORS,4*N_IN_SECTORS>(track_in_transposed[ic], mp7_in);
            if (ic < NMU) mp7_pack<N_MUON_SECTORS,8*N_IN_SECTORS>(mu_in_transposed[ic], mp7_in);
            serMP7_in(mp7_in);  
            // also make calo-only and track-only dumps for development
            for (unsigned int i = 0; i < MP7_NCHANN; ++i) mp7_in[i] = 0; // clear
            if (ic < NCALO_PER_SECTOR) mp7_pack<N_IN_SECTORS,0>(calo_in_transposed[ic], mp7_in);
            serMP7_in_calo(mp7_in);  
            for (unsigned int i = 0; i < MP7_NCHANN; ++i) mp7_in[i] = 0; // clear
            if (ic < NTRACK_PER_SECTOR/2) mp7_pack<2*N_IN_SECTORS,0>(track_in_transposed[ic], mp7_in);
            serMP7_in_track(mp7_in);  
            for (unsigned int i = 0; i < MP7_NCHANN; ++i) mp7_in[i] = 0; // clear
            if (ic < NMU) mp7_pack<N_MUON_SECTORS,0>(mu_in_transposed[ic], mp7_in);
            serMP7_in_mu(mp7_in);  
        }
#endif

        // run ref
        regionize_hadcalo(calo_fibers, calo_regions);
        regionize_hadcalo_ref(calo_fibers_ref, calo_regions_ref);
        regionize_emcalo_ref(emcalo_fibers, emcalo_regions);
        regionize_emcalo_ref(emcalo_fibers_ref, emcalo_regions_ref);
        regionize_track_ref(track_fibers, track_regions); // FIXME: I know both are _ref
        regionize_track_ref(track_fibers_ref, track_regions_ref);
        regionize_muon_ref(mu_fibers, mu_regions);
        regionize_muon_ref(mu_fibers_ref, mu_regions_ref);

        for (unsigned int ic = 0; ic < N_CLOCKS; ++ic) {
            fprintf(f_out,"Frame %04d :", ++frame_out);
            for (int i = 0; i < NCALO; ++i) {
                if (ic/2 < N_OUT_REGIONS) dump_o(f_out, calo_regions[ic/2][i]);
                else                      dump_z(f_out, calo_regions[0 ][i]);
            }
            for (int i = 0; i < NEMCALO; ++i) {
                if (ic/2 < N_OUT_REGIONS) dump_o(f_out, emcalo_regions[ic/2][i]);
                else                      dump_z(f_out, emcalo_regions[0 ][i]);
            }
            for (int i = 0; i < NTRACK; ++i) {
                if (ic/2 < N_OUT_REGIONS) dump_o(f_out, track_regions[ic/2][i]);
                else                      dump_z(f_out, track_regions[0 ][i]);
            }
            for (int i = 0; i < NMU; ++i) {
                if (ic/2 < N_OUT_REGIONS) dump_o(f_out, mu_regions[ic/2][i]);
                else                      dump_z(f_out, mu_regions[0 ][i]);
            }
            fprintf(f_out,"       %d\n", 1);
            //for (int i = 0; i < NMU; ++i) if (mu_regions[ic][i].hwPt > 0 && ic < N_OUT_REGIONS) printf("Muon  regio  %2d/%d of hwPt %6d local phi %+6d\n", ic, i,
            //    int(mu_regions[ic][i].hwPt), int(mu_regions[ic][i].hwPhi));
        }

#ifdef MP7
        for (unsigned int ic = 0; ic < N_CLOCKS/2; ++ic) {
            for (unsigned int i = 0; i < MP7_NCHANN; ++i) mp7_out[i] = 0; // clear
            if (ic < N_OUT_REGIONS) mp7_pack<NCALO,0>(calo_regions[ic], mp7_out);
            if (ic < N_OUT_REGIONS) mp7_pack<NEMCALO,2*NCALO>(emcalo_regions[ic], mp7_out);
            if (ic < N_OUT_REGIONS) mp7_pack<NTRACK,2*(NCALO+NEMCALO)>(track_regions[ic], mp7_out);
            if (ic < N_OUT_REGIONS) mp7_pack<NMU,2*(NCALO+NEMCALO+NTRACK)>(mu_regions[ic], mp7_out);
            serMP7_out(mp7_out);  
            // also make calo-only and track-only dumps for development
            for (unsigned int i = 0; i < MP7_NCHANN; ++i) mp7_out[i] = 0; // clear
            if (ic < N_OUT_REGIONS) mp7_pack<NCALO,0>(calo_regions[ic], mp7_out);
            serMP7_out_calo(mp7_out); serMP7_out_calo_nomux(mp7_out); serMP7_out_calo_nomux(mp7_out); // no-mux must be done twice
            for (unsigned int i = 0; i < MP7_NCHANN; ++i) mp7_out[i] = 0; // clear
            if (ic < N_OUT_REGIONS) mp7_pack<NTRACK,0>(track_regions[ic], mp7_out);
            serMP7_out_track(mp7_out); serMP7_out_track_nomux(mp7_out); serMP7_out_track_nomux(mp7_out);  // no-mux must be done twice
            for (unsigned int i = 0; i < MP7_NCHANN; ++i) mp7_out[i] = 0; // clear
            if (ic < N_OUT_REGIONS) mp7_pack<NMU,0>(mu_regions[ic], mp7_out);
            serMP7_out_mu(mp7_out); 
        }
#endif

 
        // -----------------------------------------
        // validation against the reference algorithm
        int errors = 0;
        for (int ir = 0; ir < N_OUT_REGIONS; ++ir) {
            for (int i = 0; i < NCALO; ++i) {
                if (!had_equals(calo_regions_ref[ir][i], calo_regions[ir][i], "regionized hadcalo", ir*100+i)) { errors++; break; }
            }
            for (int i = 0; i < NEMCALO; ++i) {
                if (!em_equals(emcalo_regions_ref[ir][i], emcalo_regions[ir][i], "regionized emcalo", ir*100+i)) { errors++; break; }
            }
            for (int i = 0; i < NTRACK; ++i) {
                if (!track_equals(track_regions_ref[ir][i], track_regions[ir][i], "regionized track", ir*100+i)) { errors++; break; }
            }
            for (int i = 0; i < NMU; ++i) {
                if (!mu_equals(mu_regions_ref[ir][i], mu_regions[ir][i], "regionized muon", ir*100+i)) { errors++; break; }
            }
        }
 
        if (errors != 0) {
            printf("Error in computing test %d (%d)\n", test, errors);
            for (int is = 0; is < N_IN_SECTORS; ++is) {
                printf("INPUT SECTOR %d (FOUND %u): \n", is, count_nonzero(calo_in[is], NCALO_PER_SECTOR));
                debug.dump_hadcalo(calo_in[is], NCALO_PER_SECTOR);
            }
            for (int ir = 0; ir < N_OUT_REGIONS; ++ir) {
                printf("OUTPUT REGION %d (ETA %d PHI %d) (REF; FOUND %u): \n", ir, ir / N_OUT_REGIONS_PHI, ir % N_OUT_REGIONS_PHI, count_nonzero(calo_regions_ref[ir], NCALO));
                debug.dump_hadcalo(calo_regions_ref[ir], NCALO);
                printf("OUTPUT REGION %d (TEST): \n", ir);
                debug.dump_hadcalo(calo_regions[ir], NCALO);
            }
            fclose(f_in); fclose(f_out);
            return 1;
        } else {
            printf("Passed test %d\n", test);
        }

    }
    fclose(f_in); fclose(f_out);
    return 0;
}
