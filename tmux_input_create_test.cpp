#include <cstdio>
#include <iomanip>
#include "firmware/simple_fullpfalgo.h"
#include "vertexing/firmware/simple_vtx.h"
#include "puppi/firmware/simple_puppi.h"
#include "utils/random_inputs.h"
#include "utils/DiscretePFInputs_IO.h"
#include "utils/pattern_serializer.h"
#include "utils/test_utils.h"

#define NTEST 36
#define NLINKS_APX_GEN0 48
#define NFRAMES_APX_GEN0 3
#define TMUX_IN 18
#define TMUX_OUT 6
int mp7DataLength = 2*(NTRACK+NCALO+NEMCALO+NMU);

int main() {

    // input format: could be random or coming from simulation
    //RandomPFInputs inputs(37); // 37 is a good random number
    DiscretePFInputs inputs("regions_TTbar_PU140.dump");
    
    // input TP objects
    HadCaloObj calo[NCALO]; EmCaloObj emcalo[NEMCALO]; TkObj track[NTRACK]; z0_t hwZPV;
    HadCaloObj calo_subem[NCALO], calo_subem_ref[NCALO]; 
    MuObj mu[NMU];

    // output PF objects
    PFChargedObj outch[NTRACK], outch_ref[NTRACK];
    PFNeutralObj outpho[NPHOTON], outpho_ref[NPHOTON];
    PFNeutralObj outne[NSELCALO], outne_ref[NSELCALO];
    PFChargedObj outmupf[NMU], outmupf_ref[NMU];

    //printf("NTRACK = %i, NEMCALO = %i, NCALO = %i, NMU = %i, MP7_NCHANN = %i \n", NTRACK, NEMCALO, NCALO, NMU, MP7_NCHANN);

    const int listLength = NFRAMES_APX_GEN0*TMUX_OUT*(((NTEST-1)/TMUX_OUT)+1);
    std::string datawords [NLINKS_APX_GEN0][listLength];
    for (int ia = 0; ia < NLINKS_APX_GEN0; ia++){
        for (int ib = 0; ib < listLength; ib++){
            datawords[ia][ib] = "0x0000000000000000";
        }
    }

    // -----------------------------------------
    // run multiple tests

    int link_ctr = 0;
    int link_off = 0;
    int offset = 0;
    for (int test = 0; test < NTEST; ++test) {
        
        offset = NFRAMES_APX_GEN0*TMUX_OUT*(test/TMUX_OUT);
        //std::cout<<offset<<std::endl;


        // initialize TP objects
        for (int i = 0; i < NTRACK; ++i) {
            track[i].hwPt = 0; track[i].hwPtErr = 0; track[i].hwEta = 0; track[i].hwPhi = 0; track[i].hwZ0 = 0; 
        }
        for (int i = 0; i < NCALO; ++i) {
            calo[i].hwPt = 0; calo[i].hwEmPt = 0; calo[i].hwEta = 0; calo[i].hwPhi = 0; calo[i].hwIsEM = 0; 
        }
        for (int i = 0; i < NEMCALO; ++i) {
            emcalo[i].hwPt = 0; emcalo[i].hwPtErr = 0;  emcalo[i].hwEta = 0; emcalo[i].hwPhi = 0;
        }
        for (int i = 0; i < NMU; ++i) {
            mu[i].hwPt = 0; mu[i].hwPtErr = 0; mu[i].hwEta = 0; mu[i].hwPhi = 0;
        }

        // get the inputs from the input object
        if (!inputs.nextRegion(calo, emcalo, track, mu, hwZPV)) break;

        VtxObj curvtx;    
        simple_vtx_ref(track,&curvtx);
        //printf("Vertex Z   %i\n",(int)(curvtx.hwZ0));

        MP7DataWord data_in[MP7_NCHANN], data_out[MP7_NCHANN];
        mp7wrapped_pack_in(emcalo, calo, track, mu, data_in);
        // for (unsigned int in = 0; in < MP7_NCHANN; in++){
        //     printf("data_in[%i] = %i \n", in, (int) data_in[in]);
        // }

        std::stringstream stream1;
        stream1 << "0x";
        stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[0]);
        stream1 << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(curvtx.hwZ0.range(9,0))) << 14) << "00";
        datawords[link_off+link_ctr][offset] = stream1.str();

        // put the data on the link number = link_ctr;
        stream1.str("");
        int id = 1;
        int index = 1;
        while (id < mp7DataLength) {
            stream1.str("");
            stream1 << "0x";
            if (index%NFRAMES_APX_GEN0==0) {
                stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]) << "00000000";
                id+=1;
            }
            else {
                stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id+1]);
                stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]);
                id+=2;
            }
            datawords[link_off+link_ctr][offset+index] = stream1.str();
            index++;
            // std::cout << stream1.str() << std::endl;
            // std::cout << datawords[link_ctr][offset+(id/2)] << std::endl;
            //printf("test = %i, link ctr = %i, clock = %i, offset = %i \n", test, link_ctr, offset+((id+1)/2), offset);
        }
        std::cout<<"link_ctr = "<<link_ctr<<" link_off = "<<link_off<<" offset = "<<offset<<" index = "<<index<<std::endl;
        
        
        link_ctr++;
        if (link_ctr >= TMUX_OUT) {
            link_off += link_ctr;
            link_ctr = 0;
        }
        if (link_off>= TMUX_IN) link_off = 0;

        //std::cout<<"-----------"<<std::endl;
    }

    for (int ib = 0; ib < listLength; ib++){
        // std::cout << ib << " ";
        std::cout << "0x" << std::setfill('0') << std::setw(4) << std::hex << ib << "  ";
        for (int ia = 0; ia < NLINKS_APX_GEN0; ia++){
            //datawords[ia][ib] = "0x0000000000000000";
            std::cout << datawords[ia][ib] << "   ";
        }
        std::cout << std::endl;
    }

    return 0;
}
