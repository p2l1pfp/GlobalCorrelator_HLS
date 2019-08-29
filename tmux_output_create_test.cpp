#include "tmux_create_test.h"

#define NETA_SMALL 2
#define NPHI_SMALL 9

int mp7DataLength = NTRACK+NCALO+NEMCALO+NMU;
int objDataLength[4] = {2*NTRACK, 2*(NEMCALO+NTRACK), 2*(NEMCALO+NCALO+NTRACK), 2*(NEMCALO+NCALO+NTRACK+NMU)};
int link_max[4] = {NLINKS_PER_TRACK, NLINKS_PER_TRACK+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO+NLINKS_PER_MU};
int link_min[4] = {0, NLINKS_PER_TRACK, NLINKS_PER_TRACK+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_EMCALO+NLINKS_PER_CALO};
unsigned int theEtaRegion = 0;
unsigned int thePhiRegion = 0;

//unsigned int outputOrder[TMUX_OUT] = {0,2,4,6,8,10,12,14,16,15,17,1,3,5,7,9,11,13};//NPHI x NETA
//unsigned int outputOrder[TMUX_OUT] = {0,11,1,12,2,13,3,14,4,15,5,16,6,17,7,9,8,10};//NPHI x NETA
unsigned int outputOrder[TMUX_OUT] = {0,2,1,3,11,4,12,5,13,6,14,7,15,8,16,9,17,10};//NPHI x NETA
//  mapping from Ryan:
//  eta, phi — 0,0 --- 1,8 
//  0,0 - 0,2 - 0,1 - 0,3 - 1,2 - 0,4 - 1,3 - 0,5 - 1,4 - 0,6 - 1,5 - 0,7 - 1,6 - 0,8 - 1,7 - 1,0 - 1,8 - 1,1
//  0     2     1     3     11    4     12    5     13    6     14    7     15    8     16    9     17    10


int main() {

    if (theEtaRegion>=NETA_TMUX) theEtaRegion = NETA_TMUX-1;
    if (thePhiRegion>=NPHI_TMUX) thePhiRegion = NPHI_TMUX-1;

    // input format: could be random or coming from simulation
    //RandomPFInputs inputs(37); // 37 is a good random number
    //DiscretePFInputs inputs("regions_TTbar_PU140.dump");
    //DiscretePFInputs inputs("barrel_sectors_1x1_TTbar_PU140.dump");
    //DiscretePFInputs inputs("barrel_sectors_1x1_TTbar_PU200.dump");
    DiscretePFInputs inputs("dummy.dump");
    
    // input TP objects
    HadCaloObj calo[NCALO_TMUX]; EmCaloObj emcalo[NEMCALO_TMUX]; TkObj track[NTRACK_TMUX]; z0_t hwZPV;
    HadCaloObj calo_subem[NCALO_TMUX], calo_subem_ref[NCALO_TMUX]; 
    MuObj mu[NMU_TMUX];

    //printf("NTRACK = %i, NEMCALO = %i, NCALO = %i, NMU = %i, MP7_NCHANN = %i \n", NTRACK, NEMCALO, NCALO, NMU, MP7_NCHANN);

    //std::cout<<mp7DataLength<<std::endl;
    const int listLength = NFRAMES_APX_GEN0*((NTEST*TMUX_OUT)+(TMUX_IN-TMUX_OUT));
    //std::cout<<listLength<<std::endl;
    std::string datawords[NTEST*TMUX_OUT][mp7DataLength+1];
    for (int ia = 0; ia < NTEST*TMUX_OUT; ia++){
        for (int ib = 0; ib < mp7DataLength+1; ib++){
            datawords[ia][ib] = "0000000000000000";
        }
    }

    // -----------------------------------------
    // run multiple tests

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
        if (!inputs.nextRegion_tmux(calo, emcalo, track, mu, hwZPV)) break;

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

        //VtxObj curvtx;    
        //simple_vtx_ref(track,&curvtx);
        //printf("Vertex Z   %i\n",(int)(curvtx.hwZ0));

        unsigned int ie = theEtaRegion;
        unsigned int ip = thePhiRegion;
        //std::cout<<"ie"<<ie<<" ip"<<ip<<std::endl;
        HadCaloObj calo_temp[TMUX_OUT][NCALO]; EmCaloObj emcalo_temp[TMUX_OUT][NEMCALO]; TkObj track_temp[TMUX_OUT][NTRACK]; MuObj mu_temp[TMUX_OUT][NMU];
        for (int ir = 0; ir < TMUX_OUT; ir++) {
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
        int etalo = -MAXETA_INT+int(float(2*MAXETA_INT*ie)/float(NETA_TMUX))-ETA_BUFFER;
        int etahi = -MAXETA_INT+int(float(2*MAXETA_INT*(ie+1))/float(NETA_TMUX))+ETA_BUFFER;
        int philo = MAXPHI_INT-NPHI_INT+int(float(NPHI_INT*ip)/float(NPHI_TMUX))-PHI_BUFFER;
        int phihi = MAXPHI_INT-NPHI_INT+int(float(NPHI_INT*(ip+1))/float(NPHI_TMUX))+PHI_BUFFER;
        //std::cout<<etalo<<" "<<etahi<<" "<<philo<<" "<<phihi<<" "<<std::endl;

	int etaremainder=(2*MAXETA_INT)%(NETA_TMUX*NETA_SMALL);
	int phiremainder=NPHI_INT%(NPHI_TMUX*NPHI_SMALL);
	int e1,e2,p1,p2;

        int i_temp[TMUX_OUT] = {0};
        int ireg = 0;
        int ntracks[TMUX_OUT] = {0};
        int ncalos[TMUX_OUT] = {0};
        int nemcalos[TMUX_OUT] = {0};
        int nmus[TMUX_OUT] = {0};
        for (int i = 0; i < NTRACK_TMUX; ++i) {
            if (int(track[i].hwEta) < etalo or int(track[i].hwEta) > etahi) continue;
            if (int(track[i].hwPhi) < philo or int(track[i].hwPhi) > phihi) continue;
            if (int(track[i].hwPt) == 0) continue;
            std::cout<<"\t"<<track[i].hwEta<<" "<<track[i].hwPhi<<std::endl;
            for (int ies = 0; ies < NETA_SMALL; ies++) {
                e1=etalo+int(float(2*MAXETA_INT)/float(NETA_TMUX*NETA_SMALL))*ies +std::min(ies,etaremainder);
	        e2=etalo+(2*ETA_BUFFER)+int(float(2*MAXETA_INT)/float(NETA_TMUX*NETA_SMALL))*(ies+1) + std::min(ies+1,etaremainder);
                std::cout<<"checking eta: "<<e1<<" "<<e2<<" (offsets "<<std::min(ies,etaremainder)<<" "<<std::min(ies+1,etaremainder)<<")"<<std::endl;
                if (int(track[i].hwEta) <= e2 and int(track[i].hwEta) > e1) {
                    for (int ips = 0; ips < NPHI_SMALL; ips++) {
                        p1=philo+int(float(2*MAXPHI_INT)/float(NPHI_TMUX*NPHI_SMALL))*ips + std::min(ips,phiremainder);
                        p2=philo+(2*PHI_BUFFER)+int(float(2*MAXPHI_INT)/float(NPHI_TMUX*NPHI_SMALL))*(ips+1) + std::min(ips+1,phiremainder);
                        std::cout<<"checking phi: "<<p1<<" "<<p2<<" (offsets "<<std::min(ips,phiremainder)<<" "<<std::min(ips+1,phiremainder)<<")"<<std::endl;
                        if ( isInPhiRegion(track[i].hwPhi, p1, p2) ) { // check "p1<=test<p2" accounting for phi wraparound
                            if (i_temp[ies*NPHI_SMALL+ips]==NTRACK) continue;
                            std::cout<<"\tX -- ("<<ies<<","<<ips<<")"<<std::endl;
                            track_temp[ies*NPHI_SMALL+ips][i_temp[ies*NPHI_SMALL+ips]] = track[i];
                            i_temp[ies*NPHI_SMALL+ips] += 1;
                            ntracks[ies*NPHI_SMALL+ips]++;
                        }
                    }
                }
            }
        }
        std::fill(i_temp, i_temp+TMUX_OUT, 0);
        for (int i = 0; i < NCALO_TMUX; ++i) {
            if (int(calo[i].hwEta) < etalo or int(calo[i].hwEta) > etahi) continue;
            if (int(calo[i].hwPhi) < philo or int(calo[i].hwPhi) > phihi) continue;
            if (int(calo[i].hwPt) == 0) continue;
            for (int ies = 0; ies < NETA_SMALL; ies++) {
                if (int(calo[i].hwEta) <= etalo+(2*ETA_BUFFER)+int(float(2*MAXETA_INT)/float(NETA_TMUX*NETA_SMALL))*(ies+1)
                and int(calo[i].hwEta) > etalo+int(float(2*MAXETA_INT)/float(NETA_TMUX*NETA_SMALL))*ies) {
                    for (int ips = 0; ips < NPHI_SMALL; ips++) {
                        if (int(calo[i].hwPhi) <= philo+(2*PHI_BUFFER)+int(float(2*MAXPHI_INT)/float(NPHI_TMUX*NPHI_SMALL))*(ips+1)
                        and int(calo[i].hwPhi) > philo+int(float(2*MAXPHI_INT)/float(NPHI_TMUX*NPHI_SMALL))*ips) {
                            if (i_temp[ies*NPHI_SMALL+ips]==NCALO) continue;
                            calo_temp[ies*NPHI_SMALL+ips][i_temp[ies*NPHI_SMALL+ips]] = calo[i];
                            i_temp[ies*NPHI_SMALL+ips] += 1;
                            ncalos[ies*NPHI_SMALL+ips]++;
                        }
                    }
                }
            }
        }
        std::fill(i_temp, i_temp+TMUX_OUT, 0);
        for (int i = 0; i < NEMCALO_TMUX; ++i) {
            if (int(emcalo[i].hwEta) < etalo or int(emcalo[i].hwEta) > etahi) continue;
            if (int(emcalo[i].hwPhi) < philo or int(emcalo[i].hwPhi) > phihi) continue;
            if (int(emcalo[i].hwPt) == 0) continue;
            for (int ies = 0; ies < NETA_SMALL; ies++) {
                if (int(emcalo[i].hwEta) <= etalo+(2*ETA_BUFFER)+int(float(2*MAXETA_INT)/float(NETA_TMUX*NETA_SMALL))*(ies+1)
                and int(emcalo[i].hwEta) > etalo+int(float(2*MAXETA_INT)/float(NETA_TMUX*NETA_SMALL))*ies) {
                    for (int ips = 0; ips < NPHI_SMALL; ips++) {
                        if (int(emcalo[i].hwPhi) <= philo+(2*PHI_BUFFER)+int(float(2*MAXPHI_INT)/float(NPHI_TMUX*NPHI_SMALL))*(ips+1)
                        and int(emcalo[i].hwPhi) > philo+int(float(2*MAXPHI_INT)/float(NPHI_TMUX*NPHI_SMALL))*ips) {
                            if (i_temp[ies*NPHI_SMALL+ips]==NEMCALO) continue;
                            emcalo_temp[ies*NPHI_SMALL+ips][i_temp[ies*NPHI_SMALL+ips]] = emcalo[i];
                            i_temp[ies*NPHI_SMALL+ips] += 1;
                            nemcalos[ies*NPHI_SMALL+ips]++;
                        }
                    }
                }
            }
        }
        std::fill(i_temp, i_temp+TMUX_OUT, 0);
        for (int i = 0; i < NMU_TMUX; ++i) {
            if (int(mu[i].hwEta) < etalo or int(mu[i].hwEta) > etahi) continue;
            if (int(mu[i].hwPhi) < philo or int(mu[i].hwPhi) > phihi) continue;
            if (int(mu[i].hwPt) == 0) continue;
            for (int ies = 0; ies < NETA_SMALL; ies++) {
                if (int(mu[i].hwEta) <= etalo+(2*ETA_BUFFER)+int(float(2*MAXETA_INT)/float(NETA_TMUX*NETA_SMALL))*(ies+1)
                and int(mu[i].hwEta) > etalo+int(float(2*MAXETA_INT)/float(NETA_TMUX*NETA_SMALL))*ies) {
                    for (int ips = 0; ips < NPHI_SMALL; ips++) {
                        if (int(mu[i].hwPhi) <= philo+(2*PHI_BUFFER)+int(float(2*MAXPHI_INT)/float(NPHI_TMUX*NPHI_SMALL))*(ips+1)
                        and int(mu[i].hwPhi) > philo+int(float(2*MAXPHI_INT)/float(NPHI_TMUX*NPHI_SMALL))*ips) {
                            if (i_temp[ies*NPHI_SMALL+ips]==NMU) continue;
                            mu_temp[ies*NPHI_SMALL+ips][i_temp[ies*NPHI_SMALL+ips]] = mu[i];
                            i_temp[ies*NPHI_SMALL+ips] += 1;
                            nmus[ies*NPHI_SMALL+ips]++;
                        }
                    }
                }
            }
        }

        for (int ir = 0; ir < TMUX_OUT; ir++) {

            /*std::cout<<"Totals: ("<<test<<", "<<ir<<")"<<std::endl;
            std::cout<<"\ttrack  = "<<ntracks[ir]<<std::endl;
            std::cout<<"\tcalo   = "<<ncalos[ir]<<std::endl;
            std::cout<<"\temcalo = "<<nemcalos[ir]<<std::endl;
            std::cout<<"\tmu     = "<<nmus[ir]<<std::endl;*/

            MP7DataWord data_in[MP7_NCHANN];
            mp7wrapped_pack_in_reorder(emcalo_temp[ir], calo_temp[ir], track_temp[ir], mu_temp[ir], data_in);
            //for (unsigned int in = 0; in < MP7_NCHANN; in++){
            //    printf("data_in[%i] = %i \n", in, (int) data_in[in]);
            //}
                

            for (int id = 0; id < mp7DataLength; id++) {
                std::stringstream stream1;
                stream1 << std::uppercase << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id*2+1]);
                stream1 << std::uppercase << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id*2+0]);
                datawords[test*TMUX_OUT+ir][id] = stream1.str();
            }
            std::stringstream stream1;
            stream1 << "00000000";
            stream1 << std::uppercase << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(hwZPV.range(9,0))) << 14) << "00";
            datawords[test*TMUX_OUT+ir][mp7DataLength] = stream1.str();
            
        }
    
    }


    int iclk = 0;
    for (int ia = 0; ia < NTEST; ia++){
        for (int io = 0; io < TMUX_OUT; io++){
            std::cout << "0x" << std::setfill('0') << std::setw(4) << std::hex << iclk << "   " <<std::dec;
            for (int ib = 0; ib < mp7DataLength+1; ib++){
                std::cout << datawords[ia*TMUX_OUT+outputOrder[io]][ib] << "    ";
            }
            std::cout << std::endl;
            iclk++;
        }
    }

    return 0;
}
