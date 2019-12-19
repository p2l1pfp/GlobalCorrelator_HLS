#include "tmux_create_test.h"

#define NETA_SMALL 2
#define NPHI_SMALL 9

int mp7DataLength = NTRACK+NCALO+NEMCALO+NMU;
int objDataLength[4] = {2*NTRACK, 2*(NEMCALO+NTRACK), 2*(NEMCALO+NCALO+NTRACK), 2*(NEMCALO+NCALO+NTRACK+NMU)};
int link_max[4] = {NLINKS_PER_TRACK, NLINKS_PER_TRACK+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO+NLINKS_PER_MU};
int link_min[4] = {0, NLINKS_PER_TRACK, NLINKS_PER_TRACK+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_EMCALO+NLINKS_PER_CALO};
unsigned int theEtaRegion = 0;
unsigned int thePhiRegion = 0;

//unsigned int outputOrder[TMUX_IN] = {0,2,4,6,8,10,12,14,16,15,17,1,3,5,7,9,11,13};//NPHI x NETA
//unsigned int outputOrder[TMUX_IN] = {0,11,1,12,2,13,3,14,4,15,5,16,6,17,7,9,8,10};//NPHI x NETA
unsigned int outputOrder[TMUX_IN] = {0,2,1,3,11,4,12,5,13,6,14,7,15,8,16,9,17,10};//NPHI x NETA
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
    DiscretePFInputs inputs("barrel_sectors_1x1_TTbar_PU200.dump");
    //DiscretePFInputs inputs("dummy.dump");
    
    // input TP objects
    HadCaloObj calo[NCALO_TMUX]; EmCaloObj emcalo[NEMCALO_TMUX]; TkObj track[NTRACK_TMUX]; z0_t hwZPV;
    HadCaloObj calo_subem[NCALO_TMUX], calo_subem_ref[NCALO_TMUX]; 
    MuObj mu[NMU_TMUX];

    //printf("NTRACK = %i, NEMCALO = %i, NCALO = %i, NMU = %i, MP7_NCHANN = %i \n", NTRACK, NEMCALO, NCALO, NMU, MP7_NCHANN);

    //std::cout<<mp7DataLength<<std::endl;
    // const int listLength = NFRAMES_APX_GEN0*((NTEST*TMUX_OUT)+(TMUX_IN-TMUX_OUT));
    //std::cout<<listLength<<std::endl;
    std::string datawords[NTEST*TMUX_IN][mp7DataLength+1];
    for (int ia = 0; ia < NTEST*TMUX_IN; ia++){
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

        //VtxObj curvtx;    
        //simple_vtx_ref(track,&curvtx);
        //printf("Vertex Z   %i\n",(int)(curvtx.hwZ0));

        unsigned int ie = theEtaRegion;
        unsigned int ip = thePhiRegion;
        //std::cout<<"ie"<<ie<<" ip"<<ip<<std::endl;
        HadCaloObj calo_temp[TMUX_IN][NCALO]; EmCaloObj emcalo_temp[TMUX_IN][NEMCALO]; TkObj track_temp[TMUX_IN][NTRACK]; MuObj mu_temp[TMUX_IN][NMU];
        for (int ir = 0; ir < TMUX_IN; ir++) {
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


        // Determine phi boundaries (include wraparound)
        std::vector<int> phi_bounds_lo{};
        std::vector<int> phi_bounds_hi{};
        const int phi_step = int(NPHI_INT)/int(NPHI_SMALL);
        const int phi_rmdr = int(NPHI_INT)%int(NPHI_SMALL);
        int p1,p2;
        for (int ip = 0; ip < NPHI_SMALL; ++ip) {
            p1 = (MAXPHI_INT-NPHI_INT) - PHI_BUFFER + phi_step*ip + std::min(ip,phi_rmdr);
            p2 = (MAXPHI_INT-NPHI_INT) + PHI_BUFFER + phi_step*(ip+1) + std::min(ip+1,phi_rmdr);
            phi_bounds_lo.push_back( p1 );
            phi_bounds_hi.push_back( p2 );
            std::cout << "TEST " << phi_bounds_lo[ip] << "  to  " << phi_bounds_hi[ip] << std::endl;
        }

        // Determine eta boundaries
        //   Want eta regions to be +/- symmetric. this implementation and these assumptions 
        //   only makes sense if ETA_TMUX==2 (all regions are doubled!)
        assert(NETA_TMUX==2);
        std::vector<int> pos_eta_bounds_lo{};
        std::vector<int> pos_eta_bounds_hi{};
        std::vector<int> eta_bounds_lo{};
        std::vector<int> eta_bounds_hi{};
        const int eta_step = int(MAXETA_INT-MINETA_INT)/int(NETA_SMALL);
        const int eta_rmdr = int(MAXETA_INT-MINETA_INT)%int(NETA_SMALL);
        for(int ie=0;ie<NETA_SMALL;++ie){
            pos_eta_bounds_lo.push_back(MINETA_INT-ETA_BUFFER+ie*eta_step+std::min(ie,eta_rmdr));
            pos_eta_bounds_hi.push_back(MINETA_INT+ETA_BUFFER+(ie+1)*eta_step+std::min(ie+1,eta_rmdr));
        }
        for(int ie=0;ie<NETA_SMALL;++ie){
            eta_bounds_lo.push_back( -pos_eta_bounds_hi[NETA_SMALL-1-ie] );
            eta_bounds_hi.push_back( -pos_eta_bounds_lo[NETA_SMALL-1-ie] );
        }
        for(int ie=0;ie<NETA_SMALL;++ie){
            eta_bounds_lo.push_back( pos_eta_bounds_lo[ie] );
            eta_bounds_hi.push_back( pos_eta_bounds_hi[ie] );
        }

        // todo - for barrel, manually implemening eta throwout
        if (eta_bounds_lo.front() < -243) eta_bounds_lo.front() = -243;
        if (eta_bounds_hi.back() > 243) eta_bounds_hi.back() = 243;

        // std::cout  << " HERE! ";
        // for (int i=0;i<eta_bounds_hi.size();i++) std::cout << eta_bounds_hi[i] << "  ";
        // std::cout << "  \n";
        eta_bounds_hi[1]=31; // clip from 32 for now

        int i_temp[TMUX_IN] = {0};
        int ireg = 0;
        int ntracks[TMUX_IN] = {0};
        int ncalos[TMUX_IN] = {0};
        int nemcalos[TMUX_IN] = {0};
        int nmus[TMUX_IN] = {0};

        int Ntracks=0;
        int Ncalos=0;
        int Nemcalos=0;
        int Nmus=0;

        for (int i = 0; i < NTRACK_TMUX; ++i) {
            if (int(track[i].hwPt) == 0) continue;
            // std::cout<<"\t"<<track[i].hwEta<<" "<<track[i].hwPhi<<std::endl;
            for (int ies = 0; ies < NETA_SMALL; ies++) {
                if (int(track[i].hwEta) >= eta_bounds_lo[ies] and int(track[i].hwEta) < eta_bounds_hi[ies]) {
                    for (int ips = 0; ips < NPHI_SMALL; ips++) {
                        if ( isInPhiRegion(track[i].hwPhi, phi_bounds_lo[ips], phi_bounds_hi[ips]) ) { 
                            // checks "p1<=test<p2" accounting for phi wraparound
                            if (i_temp[ies*NPHI_SMALL+ips]==NTRACK) continue;
                            // std::cout<<"\tX -- ("<<ies<<","<<ips<<")"<<std::endl;
                            track_temp[ies*NPHI_SMALL+ips][i_temp[ies*NPHI_SMALL+ips]] = track[i];
                            i_temp[ies*NPHI_SMALL+ips] += 1;
                            ntracks[ies*NPHI_SMALL+ips]++;
                            Ntracks++;
                        }
                    }
                }
            }
        }
        std::fill(i_temp, i_temp+TMUX_IN, 0);
        for (int i = 0; i < NCALO_TMUX; ++i) {
            //if(int(calo[i].hwPt)==49 && int(calo[i].hwEta)==32 && int(calo[i].hwPhi)==-131) std::cout << "FOUND IT AA\n"; //"000df42000070031"
            if (int(calo[i].hwPt) == 0) continue;
            //if(int(calo[i].hwPt)==49 && int(calo[i].hwEta)==32 && int(calo[i].hwPhi)==-131) std::cout << "FOUND IT BB\n"; //"000df42000070031"
            // std::cout<<"\t"<<calo[i].hwEta<<" "<<calo[i].hwPhi<<std::endl;
            for (int ies = 0; ies < NETA_SMALL; ies++) {
                if (int(calo[i].hwEta) >= eta_bounds_lo[ies] and int(calo[i].hwEta) < eta_bounds_hi[ies]) {
                    //if(int(calo[i].hwPt)==49 && int(calo[i].hwEta)==32 && int(calo[i].hwPhi)==-131) std::cout << "FOUND IT XX\n"; //"000df42000070031"
                    for (int ips = 0; ips < NPHI_SMALL; ips++) {
                        if ( isInPhiRegion(calo[i].hwPhi, phi_bounds_lo[ips], phi_bounds_hi[ips]) ) { 
                            //if(int(calo[i].hwPt)==49 && int(calo[i].hwEta)==32 && int(calo[i].hwPhi)==-131) std::cout << "FOUND IT YY\n"; //"000df42000070031"
                            if (i_temp[ies*NPHI_SMALL+ips]==NCALO) continue;
                            calo_temp[ies*NPHI_SMALL+ips][i_temp[ies*NPHI_SMALL+ips]] = calo[i];
                            i_temp[ies*NPHI_SMALL+ips] += 1;
                            ncalos[ies*NPHI_SMALL+ips]++;
                            Ncalos++;
                            //if(int(calo[i].hwPt)==49 && int(calo[i].hwEta)==32 && int(calo[i].hwPhi)==-131) std::cout << "FOUND IT ZZ\n"; //"000df42000070031"
                        }
                    }
                }
            }
        }
        std::fill(i_temp, i_temp+TMUX_IN, 0);
        for (int i = 0; i < NEMCALO_TMUX; ++i) {
            if (int(emcalo[i].hwPt) == 0) continue;
            // std::cout<<"\t"<<emcalo[i].hwEta<<" "<<emcalo[i].hwPhi<<std::endl;
            for (int ies = 0; ies < NETA_SMALL; ies++) {
                if (int(emcalo[i].hwEta) >= eta_bounds_lo[ies] and int(emcalo[i].hwEta) < eta_bounds_hi[ies]) {
                    for (int ips = 0; ips < NPHI_SMALL; ips++) {
                        if ( isInPhiRegion(emcalo[i].hwPhi, phi_bounds_lo[ips], phi_bounds_hi[ips]) ) { 
                            if (i_temp[ies*NPHI_SMALL+ips]==NEMCALO) continue;
                            emcalo_temp[ies*NPHI_SMALL+ips][i_temp[ies*NPHI_SMALL+ips]] = emcalo[i];
                            i_temp[ies*NPHI_SMALL+ips] += 1;
                            nemcalos[ies*NPHI_SMALL+ips]++;
                            Nemcalos++;
                        }
                    }
                }
            }
        }
        std::fill(i_temp, i_temp+TMUX_IN, 0);
        for (int i = 0; i < NMU_TMUX; ++i) {
            if (int(mu[i].hwPt) == 0) continue;
            for (int ies = 0; ies < NETA_SMALL; ies++) {
                if (int(mu[i].hwEta) >= eta_bounds_lo[ies] and int(mu[i].hwEta) < eta_bounds_hi[ies]) {
                    for (int ips = 0; ips < NPHI_SMALL; ips++) {
                        if ( isInPhiRegion(mu[i].hwPhi, phi_bounds_lo[ips], phi_bounds_hi[ips]) ) { 
                            if (i_temp[ies*NPHI_SMALL+ips]==NMU) continue;
                            mu_temp[ies*NPHI_SMALL+ips][i_temp[ies*NPHI_SMALL+ips]] = mu[i];
                            i_temp[ies*NPHI_SMALL+ips] += 1;
                            nmus[ies*NPHI_SMALL+ips]++;
                            Nmus++;
                        }
                    }
                }
            }
        }

        // std::cout<<"\ttrack  = "<<Ntracks<<std::endl;
        // std::cout<<"\tcalo   = "<<Ncalos<<std::endl;
        // std::cout<<"\temcalo = "<<Nemcalos<<std::endl;
        // std::cout<<"\tmu     = "<<Nmus<<std::endl;
 
        for (int ir = 0; ir < TMUX_IN; ir++) {

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
                datawords[test*TMUX_IN+ir][id] = stream1.str();
            }
            std::stringstream stream1;
            stream1 << "00000000";
            stream1 << std::uppercase << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(hwZPV.range(9,0))) << 14) << "00";
            datawords[test*TMUX_IN+ir][mp7DataLength] = stream1.str();
            
        }
    
    }


    int iclk = 0;
    for (int ia = 0; ia < NTEST; ia++){
        for (int io = 0; io < TMUX_IN; io++){
            std::cout << "0x" << std::setfill('0') << std::setw(4) << std::hex << iclk << "   " <<std::dec;
            for (int ib = 0; ib < mp7DataLength+1; ib++){
                std::cout << datawords[ia*TMUX_IN+outputOrder[io]][ib] << "    ";
            }
            std::cout << std::endl;
            iclk++;
        }
    }

    return 0;
}
