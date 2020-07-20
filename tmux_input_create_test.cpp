#include "tmux_create_test.h"

//int mp7DataLength = 2*(NTRACK+NCALO+NEMCALO+NMU);
int objDataLength[4] = {(NWORDS_TRACK*NTRACK), ((NWORDS_EMCALO*NEMCALO)+(NWORDS_TRACK*NTRACK)), ((NWORDS_EMCALO*NEMCALO)+(NWORDS_CALO*NCALO)+(NWORDS_TRACK*NTRACK)), ((NWORDS_EMCALO*NEMCALO)+(NWORDS_CALO*NCALO)+(NWORDS_TRACK*NTRACK)+(NWORDS_MU*NMU))};
int link_max[4] = {NLINKS_PER_TRACK, NLINKS_PER_TRACK+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO+NLINKS_PER_MU};
int link_min[4] = {0, NLINKS_PER_TRACK, NLINKS_PER_TRACK+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_EMCALO+NLINKS_PER_CALO};
unsigned int theEtaRegion = 0;
unsigned int thePhiRegion = 0;

void pick_link_em(int &link, l1tpf_int::CaloCluster in) {
    link++;
    if (link>=NLINKS_PER_EMCALO) link = 0;
}
void pick_link_had(int &link, l1tpf_int::CaloCluster in) {
    link++;
    if (link>=NLINKS_PER_CALO) link = 0;
}
void pick_link(int &link, l1tpf_int::PropagatedTrack in) {
    link = int((in.floatPhi()+M_PI)*float(TT_NPHI_SECTORS)/(2.*M_PI));
}
void pick_link(int &link, l1tpf_int::Muon in) {
    link++;
    if (link>=NLINKS_PER_MU) link = 0;
}

int main() {

    bool doSimple = false;
    bool debugWords = false;
    int n_alltracks = 0;
    int n_allcalos = 0;
    int n_allemcalos = 0;
    int n_allmus = 0;

    if (theEtaRegion>=NETA_TMUX) theEtaRegion = NETA_TMUX-1;
    if (thePhiRegion>=NPHI_TMUX) thePhiRegion = NPHI_TMUX-1;

    // input format: could be random or coming from simulation
    //RandomPFInputs inputs(37); // 37 is a good random number
    //DiscretePFInputs inputs("regions_TTbar_PU140.dump");
    //DiscretePFInputs inputs("barrel_sectors_1x1_TTbar_PU140.dump");
    DiscretePFInputs inputs("barrel_sectors_1x1_TTbar_PU200.dump");
    //DiscretePFInputs inputs("dummy.dump");
    
    // input TP objects
    //HadCaloObj calo_temp[TMUX_OUT][NCALO]; EmCaloObj emcalo_temp[TMUX_OUT][NEMCALO]; TkObj track_temp[TMUX_OUT][NTRACK]; 
    //MuObj mu_temp[TMUX_OUT][NMU];
    l1tpf_int::CaloCluster calo[NCALO_TMUX]; l1tpf_int::CaloCluster emcalo[NEMCALO_TMUX]; l1tpf_int::PropagatedTrack track[NTRACK_TMUX]; l1tpf_int::Muon mu[NMU_TMUX];
    z0_t hwZPV;

    //printf("NTRACK = %i, NEMCALO = %i, NCALO = %i, NMU = %i, MP7_NCHANN = %i \n", NTRACK, NEMCALO, NCALO, NMU, MP7_NCHANN);

    // NCLK_PER_BX is the number of frames per bx (320 mhz / 40mhz)
    // (320 mhz / 40mhz) * (Nevts * TM18 + (TM36-TM18)) assuming we alternate between two link groups
    //   Though the math is identical if instead we leave the links 1/3 empty and use all 3 groups
    const int listLength = NCLK_PER_BX*((NTEST*TMUX_OUT)+(TMUX_IN-TMUX_OUT));
    //std::cout<<"listLength = "<<listLength<<std::endl;
    std::string datawords [NLINKS_APX_GEN0][listLength];
    for (int ia = 0; ia < NLINKS_APX_GEN0; ia++){
        for (int ib = 0; ib < listLength; ib++){
            datawords[ia][ib] = "0x0000000000000000";
        }
    }

    // -----------------------------------------
    // run multiple tests

    int link_off = 0;
    int offset = 0;
    for (int test = 0; test < NTEST; ++test) {
        

        // initialize TP objects
        /*for (int i = 0; i < NTRACK_TMUX; ++i) {
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
        }*/
        for (int i = 0; i < NTRACK_TMUX; ++i) {
            track[i].hwPt = 0; track[i].hwEta = 0; track[i].hwPhi = 0; track[i].hwCaloPtErr = 0; track[i].hwZ0 = 0; track[i].hwStubs = 0; track[i].hwChi2 = 0;
        }
        for (int i = 0; i < NCALO_TMUX; ++i) {
            calo[i].hwPt = 0; calo[i].hwEta = 0; calo[i].hwPhi = 0; calo[i].isEM = 0; calo[i].hwEmPt = 0;
        }
        for (int i = 0; i < NEMCALO_TMUX; ++i) {
            emcalo[i].hwPt = 0; emcalo[i].hwEta = 0; emcalo[i].hwPhi = 0; emcalo[i].hwPtErr = 0;
        }
        for (int i = 0; i < NMU_TMUX; ++i) {
            mu[i].hwPt = 0; mu[i].hwEta = 0; mu[i].hwPhi = 0;
        }

        // insert some test data on the fly
        if(1){
            // get the inputs from the input object
            //if (!inputs.nextRegion_tmux(calo, emcalo, track, mu, hwZPV)) break;
            if (!inputs.nextRegion_tmux_raw(calo, emcalo, track, mu, hwZPV)) break;
        } else if(0) {
            // quickly set some fake data
            track[0].hwPt  = 1 + test * 16; // one of each object per event
            calo[0].hwPt   = 2 + test * 16;
            emcalo[0].hwPt = 3 + test * 16;
            mu[0].hwPt     = 4 + test * 16;
            hwZPV          = 7 + test * 16;
            //hwZPV          = 0;
        }

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

        std::vector<l1tpf_int::CaloCluster> calo_temp[NLINKS_PER_CALO]; std::vector<l1tpf_int::CaloCluster> emcalo_temp[NLINKS_PER_EMCALO]; std::vector<l1tpf_int::PropagatedTrack> track_temp[NLINKS_PER_TRACK]; std::vector<l1tpf_int::Muon> mu_temp[NLINKS_PER_MU];
        // initialize temp objects
        l1tpf_int::PropagatedTrack track_dummy; track_dummy.hwPt = 0; track_dummy.hwEta = 0; track_dummy.hwPhi = 0; track_dummy.hwCaloPtErr = 0; track_dummy.hwZ0 = 0; track_dummy.hwStubs = 0; track_dummy.hwChi2 = 0;
        l1tpf_int::CaloCluster calo_dummy; calo_dummy.hwPt = 0; calo_dummy.hwEta = 0; calo_dummy.hwPhi = 0; calo_dummy.isEM = 0; calo_dummy.hwEmPt = 0;
        l1tpf_int::CaloCluster emcalo_dummy; emcalo_dummy.hwPt = 0; emcalo_dummy.hwEta = 0; emcalo_dummy.hwPhi = 0; emcalo_dummy.hwPtErr = 0;
        l1tpf_int::Muon mu_dummy; mu_dummy.hwPt = 0; mu_dummy.hwEta = 0; mu_dummy.hwPhi = 0;
        // fill temp containers
        int etalo = -MAXETA_INT+int(float(2*MAXETA_INT*ie)/float(NETA_TMUX))-ETA_BUFFER;
        int etahi = -MAXETA_INT+int(float(2*MAXETA_INT*(ie+1))/float(NETA_TMUX))+ETA_BUFFER;
        int philo = -MAXPHI_INT+int(float(2*MAXPHI_INT*ip)/float(NPHI_TMUX))-PHI_BUFFER;
        int phihi = -MAXPHI_INT+int(float(2*MAXPHI_INT*(ip+1))/float(NPHI_TMUX))+PHI_BUFFER;

        int ilink = 0;
        int ntracks = 0;
        int ncalos = 0;
        int nemcalos = 0;
        int nmus = 0;
        if(1){
            ilink = NLINKS_PER_TRACK;
            for (int i = 0; i < NTRACK_TMUX; ++i) {
                if (int(track[i].hwEta) < etalo or int(track[i].hwEta) > etahi) continue;
                if (int(track[i].hwPhi) < philo or int(track[i].hwPhi) > phihi) continue;
                if (int(track[i].hwPt) == 0) continue;
                pick_link(ilink, track[i]);
                track_temp[ilink].push_back(track[i]);
                ntracks++;
            }
            ilink = NLINKS_PER_CALO;
            for (int i = 0; i < NCALO_TMUX; ++i) {
                if (int(calo[i].hwEta) < etalo or int(calo[i].hwEta) > etahi) continue;
                if (int(calo[i].hwPhi) < philo or int(calo[i].hwPhi) > phihi) continue;
                if (int(calo[i].hwPt) == 0) continue;
                pick_link_had(ilink, calo[i]);
                calo_temp[ilink].push_back(calo[i]);
                ncalos++;
            }
            ilink = NLINKS_PER_EMCALO;
            for (int i = 0; i < NEMCALO_TMUX; ++i) {
                if (int(emcalo[i].hwEta) < etalo or int(emcalo[i].hwEta) > etahi) continue;
                if (int(emcalo[i].hwPhi) < philo or int(emcalo[i].hwPhi) > phihi) continue;
                if (int(emcalo[i].hwPt) == 0) continue;
                pick_link_em(ilink, emcalo[i]);
                emcalo_temp[ilink].push_back(emcalo[i]);
                nemcalos++;
            }
            ilink = NLINKS_PER_MU;
            for (int i = 0; i < NMU_TMUX; ++i) {
                if (int(mu[i].hwEta) < etalo or int(mu[i].hwEta) > etahi) continue;
                if (int(mu[i].hwPhi) < philo or int(mu[i].hwPhi) > phihi) continue;
                if (int(mu[i].hwPt) == 0) continue;
                pick_link(ilink, mu[i]);
                mu_temp[ilink].push_back(mu[i]);
                nmus++;
            }

        } else {
            // FOR TESTING MULTIPLE INPUTS IN ONE TMUX REGION
            // ntracks and so on aren't used later, so don't bother to set them
            if(0 && test==0){
                for (int ir = 0; ir < NLINKS_PER_TRACK; ir++) {
                    track_temp [ir].push_back(track_dummy);
                    track_temp [ir][0].hwPt = 1 + ir * 16; // one of each object per TMUX_OUT
                }
                for (int ir = 0; ir < NLINKS_PER_CALO; ir++) {
                    calo_temp  [ir].push_back(calo_dummy);
                    calo_temp  [ir][0].hwPt = 2 + ir * 16;
                }
                for (int ir = 0; ir < NLINKS_PER_EMCALO; ir++) {
                    emcalo_temp[ir].push_back(emcalo_dummy);
                    emcalo_temp[ir][0].hwPt = 3 + ir * 16;
                }
                for (int ir = 0; ir < NLINKS_PER_MU; ir++) {
                    mu_temp    [ir].push_back(mu_dummy);
                    mu_temp    [ir][0].hwPt = 4 + ir * 16;
                }
            }
            if(0 && test==0){
                for(int it=0;it<NTRACK;it++) {
                    track_temp[0].push_back(track_dummy);
                    track_temp [0][it].hwPt = 1 + (it%16) * 16 + (it/16) * 16*16; // one of each object per TMUX_OUT
                }
            }
        }
        //resize to ensure number of objects can be sent in one link group
        for (int il = 0; il < NLINKS_PER_TRACK; il++) {
            track_temp[il].resize(((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_TRACK, track_dummy);
        }
        for (int il = 0; il < NLINKS_PER_CALO; il++) {
            calo_temp[il].resize(((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_CALO, calo_dummy);
        }
        for (int il = 0; il < NLINKS_PER_EMCALO; il++) {
            emcalo_temp[il].resize(((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_EMCALO, emcalo_dummy);
        }
        for (int il = 0; il < NLINKS_PER_MU; il++) {
            mu_temp[il].resize(((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_MU, mu_dummy);
        }

        /*std::cout<<"Totals:"<<std::endl;
        std::cout<<"\ttrack  = "<<ntracks<<std::endl;
        std::cout<<"\tcalo   = "<<ncalos<<std::endl;
        std::cout<<"\temcalo = "<<nemcalos<<std::endl;
        std::cout<<"\tmu     = "<<nmus<<std::endl;*/
        n_alltracks  += ntracks;
        n_allcalos   += ncalos;
        n_allemcalos += nemcalos;
        n_allmus     += nmus;    


        offset = NCLK_PER_BX*TMUX_OUT*test; // 8 * 18 * nevt
        //std::cout<<offset<<std::endl;

        unsigned int link_type = 0; //0=track, 1=emcalo, 2=calo, 3=mu
        unsigned int obj_link_no;

        for (int link_ctr = 0; link_ctr < NLINKS_PER_REG; link_ctr++) {

            if      (link_ctr < link_max[0]) link_type = 0;
            else if (link_ctr < link_max[1]) link_type = 1;
            else if (link_ctr < link_max[2]) link_type = 2;
            else if (link_ctr < link_max[3]) link_type = 3;
            obj_link_no = link_ctr-link_min[link_type];

            if (link_type == 0) {
                write_track_vector_to_link(track_temp[obj_link_no], datawords[link_off+link_ctr], offset, obj_link_no); 
            } else if (link_type == 1) {
                write_emcalo_vector_to_link(emcalo_temp[obj_link_no], datawords[link_off+link_ctr], offset); 
            } else if (link_type == 2) {
                write_calo_vector_to_link(calo_temp[obj_link_no], datawords[link_off+link_ctr], offset);
            } else if (link_type == 3) {
                write_mu_vector_to_link(mu_temp[obj_link_no], datawords[link_off+link_ctr], offset);
            }
        }
        
        std::stringstream stream2;
        stream2.str("");
        stream2 << "0x00000000" << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(hwZPV.range(9,0))) << 14) << "00"; 
        datawords[link_off+NLINKS_PER_REG][offset] = stream2.str();
        

        link_off += NLINKS_PER_REG+1;
        // 2 schemes should be equivalent in 18 / 6 setup
        // in this scheme, the first link is maximally saturated
        if (0 && link_off>= (NLINKS_PER_REG+1)*TMUX_IN/TMUX_OUT) link_off = 0;
        // in this scheme, inputs are spread across all links
        if (1 && link_off+(NLINKS_PER_REG+1) > NLINKS_APX_GEN0) link_off = 0;
    
    }
    std::ofstream outfile;
    outfile.open("../../../../inputs.txt");
    for (int ib = 0; ib < listLength; ib++){
        // std::cout << ib << " ";
        outfile << "0x" << std::setfill('0') << std::setw(4) << std::hex << ib << "   " <<std::dec;
        //for (int ia = 0; ia < NLINKS_APX_GEN0; ia++){
        for (int ia = NLINKS_APX_GEN0-1; ia >=0; ia--){
            //datawords[ia][ib] = "0x0000000000000000";
            outfile << datawords[ia][ib] << "    ";
        }
        outfile << std::endl;
    }
    outfile.close();

    // std::cout<<"For all events: "<<std::endl;
    // std::cout<<"\ttrack  = "<<n_alltracks<<std::endl;
    // std::cout<<"\tcalo   = "<<n_allcalos<<std::endl;
    // std::cout<<"\temcalo = "<<n_allemcalos<<std::endl;
    // std::cout<<"\tmu     = "<<n_allmus<<std::endl;

    return 0;
}
