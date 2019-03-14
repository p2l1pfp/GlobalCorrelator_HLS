#include <cstdio>
#include <iomanip>
#include "firmware/simple_fullpfalgo.h"
#include "vertexing/firmware/simple_vtx.h"
#include "puppi/firmware/simple_puppi.h"
#include "utils/random_inputs.h"
#include "utils/DiscretePFInputs_IO.h"
#include "utils/pattern_serializer.h"
#include "utils/test_utils.h"

#define NTEST 2
#define NLINKS_APX_GEN0 48
#define NFRAMES_APX_GEN0 3
#define TMUX_IN 18
#define TMUX_OUT 6

#define MAXETA_INT 345
#define MAXPHI_INT 720

int mp7DataLength = 2*(NTRACK+NCALO+NEMCALO+NMU);

int main() {

    // input format: could be random or coming from simulation
    //RandomPFInputs inputs(37); // 37 is a good random number
    //DiscretePFInputs inputs("regions_TTbar_PU140.dump");
    DiscretePFInputs inputs("barrel_sectors_1x1_TTbar_PU140.dump");
    
    // input TP objects
    HadCaloObj calo[NCALO_TMUX]; EmCaloObj emcalo[NEMCALO_TMUX]; TkObj track[NTRACK_TMUX]; z0_t hwZPV;
    HadCaloObj calo_subem[NCALO_TMUX], calo_subem_ref[NCALO_TMUX]; 
    MuObj mu[NMU_TMUX];

    //printf("NTRACK = %i, NEMCALO = %i, NCALO = %i, NMU = %i, MP7_NCHANN = %i \n", NTRACK, NEMCALO, NCALO, NMU, MP7_NCHANN);

    std::cout<<mp7DataLength<<std::endl;
    const int listLength = NFRAMES_APX_GEN0*((NTEST*NREGIONS_TMUX*NETA_TMUX*NPHI_TMUX)+(TMUX_IN-TMUX_OUT));
    std::cout<<listLength<<std::endl;
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
        

        // initialize TP objects
        for (int i = 0; i < NTRACK_TMUX; ++i) {
            track[i].hwPt = 0; track[i].hwPtErr = 0; track[i].hwEta = 0; track[i].hwPhi = 0; track[i].hwZ0 = 0; 
        }
        for (int i = 0; i < NCALO_TMUX; ++i) {
            calo[i].hwPt = 0; calo[i].hwEmPt = 0; calo[i].hwEta = 0; calo[i].hwPhi = 0; calo[i].hwIsEM = 0; 
        }
        for (int i = 0; i < NEMCALO_TMUX; ++i) {
            emcalo[i].hwPt = 0; emcalo[i].hwPtErr = 0;  emcalo[i].hwEta = 0; emcalo[i].hwPhi = 0;
        }
        for (int i = 0; i < NMU_TMUX; ++i) {
            mu[i].hwPt = 0; mu[i].hwPtErr = 0; mu[i].hwEta = 0; mu[i].hwPhi = 0;
        }

        // get the inputs from the input object
        if (!inputs.nextRegion(calo, emcalo, track, mu, hwZPV)) break;

        /*for (int i = 0; i < NTRACK_TMUX; ++i) {
            std::cout<<track[i].hwPt<<"\t "<<track[i].hwEta<<"\t "<<track[i].hwPhi<<std::endl;
        }
        for (int i = 0; i < NCALO_TMUX; ++i) {
            std::cout<<calo[i].hwPt<<"\t "<<calo[i].hwEta<<"\t "<<calo[i].hwPhi<<std::endl;
        }
        for (int i = 0; i < NEMCALO_TMUX; ++i) {
            std::cout<<emcalo[i].hwPt<<"\t "<<emcalo[i].hwEta<<"\t "<<emcalo[i].hwPhi<<std::endl;
        }
        for (int i = 0; i < NMU_TMUX; ++i) {
            std::cout<<mu[i].hwPt<<"\t "<<mu[i].hwEta<<"\t "<<mu[i].hwPhi<<std::endl;
        }*/

        VtxObj curvtx;    
        simple_vtx_ref(track,&curvtx);
        //printf("Vertex Z   %i\n",(int)(curvtx.hwZ0));

        for (unsigned int ie = 0; ie < NETA_TMUX; ie++) {
            for (unsigned int ip = 0; ip < NPHI_TMUX; ip++) {
                //std::cout<<"ie"<<ie<<" ip"<<ip<<std::endl;
                HadCaloObj calo_temp[NREGIONS_TMUX][NCALO]; EmCaloObj emcalo_temp[NREGIONS_TMUX][NEMCALO]; TkObj track_temp[NREGIONS_TMUX][NTRACK]; MuObj mu_temp[NREGIONS_TMUX][NMU];
                for (int ir = 0; ir < NREGIONS_TMUX; ir++) {
                    // initialize temp objects
                    for (int i = 0; i < NTRACK; ++i) {
                        track_temp[ir][i].hwPt = 0; track_temp[ir][i].hwPtErr = 0; track_temp[ir][i].hwEta = 0; track_temp[ir][i].hwPhi = 0; track_temp[ir][i].hwZ0 = 0; 
                    }
                    for (int i = 0; i < NCALO; ++i) {
                        calo_temp[ir][i].hwPt = 0; calo_temp[ir][i].hwEmPt = 0; calo_temp[ir][i].hwEta = 0; calo_temp[ir][i].hwPhi = 0; calo_temp[ir][i].hwIsEM = 0; 
                    }
                    for (int i = 0; i < NEMCALO; ++i) {
                        emcalo_temp[ir][i].hwPt = 0; emcalo_temp[ir][i].hwPtErr = 0;  emcalo_temp[ir][i].hwEta = 0; emcalo_temp[ir][i].hwPhi = 0;
                    }
                    for (int i = 0; i < NMU; ++i) {
                        mu_temp[ir][i].hwPt = 0; mu_temp[ir][i].hwPtErr = 0; mu_temp[ir][i].hwEta = 0; mu_temp[ir][i].hwPhi = 0;
                    }
                }
                // fill temp containers
                int etalo = -MAXETA_INT+int(float(2*MAXETA_INT*ie)/float(NETA_TMUX));
                int etahi = -MAXETA_INT+int(float(2*MAXETA_INT*(ie+1))/float(NETA_TMUX));
                int philo = -MAXPHI_INT+int(float(2*MAXPHI_INT*ip)/float(NPHI_TMUX));
                int phihi = -MAXPHI_INT+int(float(2*MAXPHI_INT*(ip+1))/float(NPHI_TMUX));
                //std::cout<<etalo<<" "<<etahi<<" "<<philo<<" "<<phihi<<" "<<std::endl;

                int i_temp = 0;
                int ireg = 0;
                for (int i = 0; i < NTRACK_TMUX; ++i) {
                    if (int(track[i].hwEta) < etalo or int(track[i].hwEta) > etahi) continue;
                    if (int(track[i].hwPhi) < philo or int(track[i].hwPhi) > phihi) continue;
                    if (int(track[i].hwPt) == 0) continue;
                    track_temp[ireg][i_temp] = track[i];
                    ireg++;
                    if (ireg == NREGIONS_TMUX) {i_temp++; ireg=0;}
                    if (i_temp == NTRACK) {break;}
                }
                i_temp = 0;
                ireg = 0;
                for (int i = 0; i < NCALO_TMUX; ++i) {
                    if (int(calo[i].hwEta) < etalo or int(calo[i].hwEta) > etahi) continue;
                    if (int(calo[i].hwPhi) < philo or int(calo[i].hwPhi) > phihi) continue;
                    if (int(calo[i].hwPt) == 0) continue;
                    calo_temp[ireg][i_temp] = calo[i];
                    ireg++;
                    if (ireg == NREGIONS_TMUX) {i_temp++; ireg=0;}
                    if (i_temp == NCALO) {break;}
                }
                i_temp = 0;
                ireg = 0;
                for (int i = 0; i < NEMCALO_TMUX; ++i) {
                    if (int(emcalo[i].hwEta) < etalo or int(emcalo[i].hwEta) > etahi) continue;
                    if (int(emcalo[i].hwPhi) < philo or int(emcalo[i].hwPhi) > phihi) continue;
                    if (int(emcalo[i].hwPt) == 0) continue;
                    emcalo_temp[ireg][i_temp] = emcalo[i];
                    ireg++;
                    if (ireg == NREGIONS_TMUX) {i_temp++; ireg=0;}
                    if (i_temp == NEMCALO) {break;}
                }
                i_temp = 0;
                ireg = 0;
                for (int i = 0; i < NMU_TMUX; ++i) {
                    if (int(mu[i].hwEta) < etalo or int(mu[i].hwEta) > etahi) continue;
                    if (int(mu[i].hwPhi) < philo or int(mu[i].hwPhi) > phihi) continue;
                    if (int(mu[i].hwPt) == 0) continue;
                    mu_temp[ireg][i_temp] = mu[i];
                    ireg++;
                    if (ireg == NREGIONS_TMUX) {i_temp++; ireg=0;}
                    if (i_temp == NMU) {break;}
                }

                for (int ir = 0; ir < NREGIONS_TMUX; ir++) {

                    MP7DataWord data_in[MP7_NCHANN];
                    mp7wrapped_pack_in(emcalo_temp[ir], calo_temp[ir], track_temp[ir], mu_temp[ir], data_in);
                    //for (unsigned int in = 0; in < MP7_NCHANN; in++){
                    //    printf("data_in[%i] = %i \n", in, (int) data_in[in]);
                    //}
    
                    offset = NFRAMES_APX_GEN0*TMUX_OUT*((test*NETA_TMUX*NPHI_TMUX*NREGIONS_TMUX+ie*NPHI_TMUX*NREGIONS_TMUX+ip*NREGIONS_TMUX+ir)/TMUX_OUT);
                    //std::cout<<offset<<std::endl;
    
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
                        //std::cout<<"link_ctr = "<<link_ctr<<" link_off = "<<link_off<<" offset = "<<offset<<" index = "<<index<<std::endl;
                    }
                    
                    link_ctr++;
                    if (link_ctr >= TMUX_OUT) {
                        link_off += link_ctr;
                        link_ctr = 0;
                    }
                    if (link_off>= TMUX_IN) link_off = 0;
    
                    //std::cout<<"-----------"<<std::endl;
                }
            }
        }
    }

    for (int ib = 0; ib < listLength; ib++){
        // std::cout << ib << " ";
        std::cout << "0x" << std::setfill('0') << std::setw(4) << std::hex << ib << "  " <<std::dec;
        for (int ia = 0; ia < NLINKS_APX_GEN0; ia++){
            //datawords[ia][ib] = "0x0000000000000000";
            std::cout << datawords[ia][ib] << "   ";
        }
        std::cout << std::endl;
    }

    return 0;
}
