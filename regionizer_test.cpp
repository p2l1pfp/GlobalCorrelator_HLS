#include <cstdio>
#include "firmware/regionizer.h"
#include "random_inputs.h"
#include "DiscretePFInputs_IO.h"
#include "pattern_serializer.h"
#include "test_utils.h"

#define NTEST 5

template<unsigned int N, typename T>
bool fill_stream(hls::stream<T> & stream, T data[N], const char *name="unknown", int isector=-999) {
    if (!stream.empty()) { printf("ERROR: %s stream for sector %d is not empty\n", name, isector); return false; }
    for (unsigned int i = 0; i < N; ++i) stream.write(data[i]);
    return true;
}

void dump_o(FILE *f, const HadCaloObj & obj) { 
    // NOTE: no leading + on positive numbers, VHDL doesn't like it
    fprintf(f, "       %4d % 4d % 4d %3d %1d", int(obj.hwPt), int(obj.hwEta), int(obj.hwPhi), int(obj.hwEmPt), obj.hwIsEM); 
}
void dump_z(FILE *f, const HadCaloObj &) { // the second argument is only to resolve overloading
    HadCaloObj dummy; clear(dummy); dump_o(f,dummy);
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
    printf(" Regions: %d (%d x %d )\n", N_OUT_REGIONS, N_OUT_REGIONS_ETA, N_OUT_REGIONS_PHI);
    printf("    max N(Calo): %2d \n", NCALO);

    HadCaloObj calo_in[N_IN_SECTORS][NCALO_PER_SECTOR];
    hls::stream<HadCaloObj> calo_fibers[N_IN_SECTORS];
    hls::stream<HadCaloObj> calo_fibers_ref[N_IN_SECTORS];
    HadCaloObj calo_regions[N_OUT_REGIONS][NCALO]; 
    HadCaloObj calo_regions_ref[N_OUT_REGIONS][NCALO]; 

    FILE *f_in  = fopen("dump_in.txt","w");
    FILE *f_out = fopen("dump_out.txt","w");
    int frame_in = 0, frame_out = 0;
    // -----------------------------------------
    // run multiple tests
    for (int test = 1; test <= NTEST; ++test) {
        // read the event
        if (!inputs.nextEvent()) break;
        if (inputs.event().regions.size() != N_IN_SECTORS) { printf("ERROR: Mismatching number of input regions: %lu\n", inputs.event().regions.size()); return 2; }
        //fill in the streams
        for (int is = 0; is < N_IN_SECTORS; ++is) {
            const Region & r = inputs.event().regions[is];
            dpf2fw::convert<NCALO_PER_SECTOR>(r.calo, calo_in[is]); 
            for (unsigned int i = 0; i < NCALO_PER_SECTOR; ++i) assert(calo_in[is][i].hwPt >= 0);
            if (!fill_stream<NCALO_PER_SECTOR>(calo_fibers[is], calo_in[is], "calo stream", is)) return 3;
            if (!fill_stream<NCALO_PER_SECTOR>(calo_fibers_ref[is], calo_in[is], "calo ref stream", is)) return 3;
        }
        // dump inputs
        for (unsigned int ic = 0; ic < N_CLOCKS_SECTOR; ++ic) {
            for (unsigned int id = 0; id < 2; ++id) { // takes 2 clocks to send one input; so we just duplicate the lines for now
                fprintf(f_in,"Frame %04d : %2d %2d", ++frame_in, test, ic);
                for (int is = 0; is < N_IN_SECTORS; ++is) {
                    if (ic < NCALO_PER_SECTOR && !id) dump_o(f_in, calo_in[is][ic]);
                    else                              dump_z(f_in, calo_in[is][0]);
                }
                fprintf(f_in,"\n");
            }
        }

        // run ref
        regionize_hadcalo(calo_fibers, calo_regions);
        regionize_hadcalo_ref(calo_fibers_ref, calo_regions_ref);

        for (unsigned int ic = 0; ic < N_CLOCKS_SECTOR; ++ic) {
            fprintf(f_out,"Frame %04d :", ++frame_out);
            for (int i = 0; i < NCALO; ++i) {
                if (ic < N_OUT_REGIONS) dump_o(f_out, calo_regions[ic][i]);
                else                    dump_z(f_out, calo_regions[0 ][i]);
            }
            fprintf(f_out,"       %d\n", 1);
        }
 
        // -----------------------------------------
        // validation against the reference algorithm
        int errors = 0;
        for (int ir = 0; ir < N_OUT_REGIONS; ++ir) {
            for (int i = 0; i < NCALO; ++i) {
                if (!had_equals(calo_regions_ref[ir][i], calo_regions[ir][i], "regionized had calo", ir*100+i)) { errors++; break; }
            }
        }
 
        int n_expected[N_OUT_REGIONS];
        for (int ir = 0; ir < N_OUT_REGIONS; ++ir) { n_expected[ir] = 0; }
        for (int io = 0; io < NCALO_PER_SECTOR; ++io) {
            for (int is = 0; is < N_IN_SECTORS; ++is) {
                int phi0s = (1 + 2*is - 12)*_PHI_PIO6/2; // pi/12 + is*pi/6 - pi
                if (calo_in[is][io].hwPt == 0) continue;
                for (unsigned int ie = 0; ie < N_OUT_REGIONS_ETA; ++ie) {
                    if (!(ETA_MIN[ie] <= int(calo_in[is][io].hwEta) && int(calo_in[is][io].hwEta) <= ETA_MAX[ie])) continue;
                    for (int ip = 0; ip < N_OUT_REGIONS_PHI; ++ip) {
                        unsigned int ir = N_OUT_REGIONS_PHI*ie + ip;
                        int phi0r = (3 + 4*ip - 12)*_PHI_PIO6/2; // pi/4 + ip*pi/3 - pi // NOTE the offset is chosen so that the boundaries of the region, including the padding, align with the sectors borders
                        int dphi = int(calo_in[is][io].hwPhi) + phi0s - phi0r; 
                        while (dphi >  6*_PHI_PIO6) dphi -= 12*_PHI_PIO6;
                        while (dphi < -6*_PHI_PIO6) dphi += 12*_PHI_PIO6;
                        bool by_cabling = false;
                        for (unsigned int ic = 0; ic < 3; ++ic) if (IN_SECTOR_OF_REGION[ip][ic] == is) by_cabling = true;
                        bool by_dphi = (std::abs(dphi) <= 3*_PHI_PIO6/2);
                        if (by_cabling) n_expected[ir]++;
                        if (by_cabling != by_dphi) {
                            printf("LOGIC ERROR in phi mapping for region %d (iphi %d) : by cabling %d, by dphi %d\n", ir, ip, by_cabling, by_dphi); 
                            printf("object local  iphi in sector: %+6d\n", int(calo_in[is][io].hwPhi)); 
                            printf("object global iphi:           %+6d\n", int(calo_in[is][io].hwPhi) + phi0s); 
                            printf("region global iphi:           %+6d\n", phi0r); 
                            printf("                pi:           %+6d\n", 6*_PHI_PIO6); 
                            printf("              2*pi:           %+6d\n", 12*_PHI_PIO6); 
                            printf("raw     local dphi in region: %+6d\n", int(calo_in[is][io].hwPhi) + phi0s - phi0r); 
                            printf("wrapped local dphi in region: %+6d\n", dphi); 
                            printf("region half size  :           %+6d\n", 3*_PHI_PIO6/2); 
                            printf("object global phi * PI/12:  %+.3f   [-12 to 12 range]\n", (int(calo_in[is][io].hwPhi) + phi0s)/float(_PHI_PIO6/2)); 
                            return 37;
                        }
                        //if (by_cabling) printf("Object %d in sector %d, ieta %+3d (eta %+.3f) extected in region %d (eta %d phi %d)\n", io, is, int(calo_in[is][io].hwEta), calo_in[is][io].hwEta*0.25/_ETA_025, ir, ie, ip);
                    }
                }
            }
        }
        for (int ir = 0; ir < N_OUT_REGIONS; ++ir) { 
            if (std::min<int>(n_expected[ir],NCALO) != count_nonzero(calo_regions_ref[ir], NCALO)) errors++; 
        }

        if (errors != 0) {
            printf("Error in computing test %d (%d)\n", test, errors);
            for (int is = 0; is < N_IN_SECTORS; ++is) {
                printf("INPUT SECTOR %d (FOUND %u): \n", is, count_nonzero(calo_in[is], NCALO_PER_SECTOR));
                debug.dump_hadcalo(calo_in[is], NCALO_PER_SECTOR);
            }
            for (int ir = 0; ir < N_OUT_REGIONS; ++ir) {
                printf("OUTPUT REGION %d (ETA %d PHI %d) (REF; EXPECTED: %u; FOUND %u): \n", ir, ir / N_OUT_REGIONS_PHI, ir % N_OUT_REGIONS_PHI, n_expected[ir], count_nonzero(calo_regions_ref[ir], NCALO));
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
