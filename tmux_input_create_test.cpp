#include "tmux_create_test.h"

int mp7DataLength = 2*(NTRACK+NCALO+NEMCALO+NMU);
int objDataLength[4] = {2*NTRACK, 2*(NEMCALO+NTRACK), 2*(NEMCALO+NCALO+NTRACK), 2*(NEMCALO+NCALO+NTRACK+NMU)};
int link_max[4] = {NLINKS_PER_TRACK, NLINKS_PER_TRACK+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO+NLINKS_PER_MU};
int link_min[4] = {0, NLINKS_PER_TRACK, NLINKS_PER_TRACK+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_EMCALO+NLINKS_PER_CALO};
unsigned int theEtaRegion = 0;
unsigned int thePhiRegion = 0;

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
    HadCaloObj calo[NCALO_TMUX]; EmCaloObj emcalo[NEMCALO_TMUX]; TkObj track[NTRACK_TMUX]; z0_t hwZPV;
    HadCaloObj calo_subem[NCALO_TMUX], calo_subem_ref[NCALO_TMUX]; 
    MuObj mu[NMU_TMUX];

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


        // insert some test data on the fly
        if(1){
            // get the inputs from the input object
            if (!inputs.nextRegion_tmux(calo, emcalo, track, mu, hwZPV)) break;
            // if (!inputs.nextRegion_tmux(calo, emcalo, track, mu, hwZPV)) break;
            // if (!inputs.nextRegion_tmux(calo, emcalo, track, mu, hwZPV)) break;
            // if (!inputs.nextRegion_tmux(calo, emcalo, track, mu, hwZPV)) break;
            // if (!inputs.nextRegion_tmux(calo, emcalo, track, mu, hwZPV)) break;
            // if (!inputs.nextRegion_tmux(calo, emcalo, track, mu, hwZPV)) break;
        } else if(0) {
            track[0].hwPt  = 1 + test * 16; // one of each object per event
            calo[0].hwPt   = 2 + test * 16;
            emcalo[0].hwPt = 3 + test * 16;
            mu[0].hwPt     = 4 + test * 16;
            hwZPV          = 7 + test * 16;
            //hwZPV          = 0;
        }

        //hwZPV          = 7;


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
        int philo = -MAXPHI_INT+int(float(2*MAXPHI_INT*ip)/float(NPHI_TMUX))-PHI_BUFFER;
        int phihi = -MAXPHI_INT+int(float(2*MAXPHI_INT*(ip+1))/float(NPHI_TMUX))+PHI_BUFFER;
        //std::cout<<etalo<<" "<<etahi<<" "<<philo<<" "<<phihi<<" "<<std::endl;

        int i_temp = 0;
        int ireg = 0;
        int ntracks = 0;
        int ncalos = 0;
        int nemcalos = 0;
        int nmus = 0;
        if(1){
            for (int i = 0; i < NTRACK_TMUX; ++i) {
                if (int(track[i].hwEta) < etalo or int(track[i].hwEta) > etahi) continue;
                if (int(track[i].hwPhi) < philo or int(track[i].hwPhi) > phihi) continue;
                if (int(track[i].hwPt) == 0) continue;
                track_temp[ireg][i_temp] = track[i];
                ntracks++;
                ireg++;
                if (ireg == TMUX_OUT) {i_temp++; ireg=0;}
                if (i_temp == NTRACK) {break;}
            }
            i_temp = 0;
            ireg = 0;
            for (int i = 0; i < NCALO_TMUX; ++i) {
                if (int(calo[i].hwEta) < etalo or int(calo[i].hwEta) > etahi) continue;
                if (int(calo[i].hwPhi) < philo or int(calo[i].hwPhi) > phihi) continue;
                if (int(calo[i].hwPt) == 0) continue;
                calo_temp[ireg][i_temp] = calo[i];
                ncalos++;
                ireg++;
                if (ireg == TMUX_OUT) {i_temp++; ireg=0;}
                if (i_temp == NCALO) {break;}
            }
            i_temp = 0;
            ireg = 0;
            for (int i = 0; i < NEMCALO_TMUX; ++i) {
                if (int(emcalo[i].hwEta) < etalo or int(emcalo[i].hwEta) > etahi) continue;
                if (int(emcalo[i].hwPhi) < philo or int(emcalo[i].hwPhi) > phihi) continue;
                if (int(emcalo[i].hwPt) == 0) continue;
                emcalo_temp[ireg][i_temp] = emcalo[i];
                nemcalos++;
                ireg++;
                if (ireg == TMUX_OUT) {i_temp++; ireg=0;}
                if (i_temp == NEMCALO) {break;}
            }
            i_temp = 0;
            ireg = 0;
            for (int i = 0; i < NMU_TMUX; ++i) {
                if (int(mu[i].hwEta) < etalo or int(mu[i].hwEta) > etahi) continue;
                if (int(mu[i].hwPhi) < philo or int(mu[i].hwPhi) > phihi) continue;
                if (int(mu[i].hwPt) == 0) continue;
                mu_temp[ireg][i_temp] = mu[i];
                nmus++;
                ireg++;
                if (ireg == TMUX_OUT) {i_temp++; ireg=0;}
                if (i_temp == NMU) {break;}
            }

        } else {
            // FOR TESTING MULTIPLE INPUTS IN ONE TMUX REGION
            // ntracks and so on aren't used later, so don't bother to set them
            if(0 && test==0){
                for (int ir = 0; ir < TMUX_OUT; ir++) {
                    track_temp [ir][0].hwPt = 1 + ir * 16; // one of each object per TMUX_OUT
                    calo_temp  [ir][0].hwPt = 2 + ir * 16;
                    emcalo_temp[ir][0].hwPt = 3 + ir * 16;
                    mu_temp    [ir][0].hwPt = 4 + ir * 16;
                }
            }
            if(0 && test==0){
                for(int it=0;it<NTRACK;it++)
                    track_temp [0][it].hwPt = 1 + (it%16) * 16 + (it/16) * 16*16; // one of each object per TMUX_OUT
            }
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

        for (int ir = 0; ir < TMUX_OUT; ir++) {

            MP7DataWord data_in[MP7_NCHANN];
            mp7wrapped_pack_in_reorder(emcalo_temp[ir], calo_temp[ir], track_temp[ir], mu_temp[ir], data_in);
            //for (unsigned int in = 0; in < MP7_NCHANN; in++){
            //    printf("data_in[%i] = %i \n", in, (int) data_in[in]);
            //}

            offset = NCLK_PER_BX*TMUX_OUT*test; // 8 * 18 * nevt
            //std::cout<<offset<<std::endl;

            float tot_perc[4];
            int link_start[4];
            int add_off[4];

            tot_perc[0] = float(ir*NLINKS_PER_TRACK)/float(TMUX_OUT);
            tot_perc[1] = float(ir*NLINKS_PER_EMCALO)/float(TMUX_OUT);
            tot_perc[2] = float(ir*NLINKS_PER_CALO)/float(TMUX_OUT);
            tot_perc[3] = float(ir*NLINKS_PER_MU)/float(TMUX_OUT);
            link_start[0] = int(tot_perc[0]);
            link_start[1] = int(tot_perc[1]);
            link_start[2] = int(tot_perc[2]);
            link_start[3] = int(tot_perc[3]);
            add_off[0] = int(float(NCLK_PER_BX*TMUX_IN)*tot_perc[0])%(NCLK_PER_BX*TMUX_IN); // 8 * 36
            add_off[1] = int(float(NCLK_PER_BX*TMUX_IN)*tot_perc[1])%(NCLK_PER_BX*TMUX_IN);
            add_off[2] = int(float(NCLK_PER_BX*TMUX_IN)*tot_perc[2])%(NCLK_PER_BX*TMUX_IN);
            add_off[3] = int(float(NCLK_PER_BX*TMUX_IN)*tot_perc[3])%(NCLK_PER_BX*TMUX_IN);

            int id = 0;
            unsigned int link_type = 0; //0=track, 1=emcalo, 2=calo, 3=mu

            for (int link_ctr = 0; link_ctr < NLINKS_PER_REG; link_ctr++) {

                if      (link_ctr < link_max[0]) link_type = 0;
                else if (link_ctr < link_max[1]) link_type = 1;
                else if (link_ctr < link_max[2]) link_type = 2;
                else if (link_ctr < link_max[3]) link_type = 3;

                //std::cout<<"Region "<<ir<<" Link "<<link_ctr<<" : "<<tot_perc[link_type]<<" "<<link_start[link_type]<<" "<<add_off[link_type]<<std::endl;

                if (link_ctr < link_min[link_type]+link_start[link_type]) {continue;}
                
                std::stringstream stream1;
                std::stringstream tmp1,tmp2;
                int index = add_off[link_type]; // 0 for small region 0, else (ir/18*10)*(8*18) mod 8*18 or whatever
                if(doSimple){
                    // don't send 32b data in this word along with the vertex. (keeps things aligned later on)
                    stream1 << "0x";
                    stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (0);
                    stream1 << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(hwZPV.range(9,0))) << 14) << "00";
                    datawords[link_off+link_ctr][offset] = stream1.str();
                    index++;
                } else {
                    if (index==0) { // if first entry (clk of TMIN*CLKpBX) of any of any link
                        tmp1.str(""); tmp2.str("");
                        stream1 << "0x";
                        stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]); 
                        tmp1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]);
                        // zero first 8b for beginning of frame!
                        stream1 << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(hwZPV.range(9,0))) << 14) << "00"; 
                        tmp2 << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(hwZPV.range(9,0))) << 14) << "00";
                        datawords[link_off+link_ctr][offset] = stream1.str();
                        if(debugWords) printf("datawords[%d][%d] = %s (%s + %s) <-- index 0\n",link_off+link_ctr,offset,datawords[link_off+link_ctr][offset].c_str(),tmp1.str().c_str(),tmp2.str().c_str());
                        id++;
                        index++;
                    }
                    else {
                        // this seems like it does nothing since the same thing is written to the same data word as above...
                        tmp1.str(""); tmp2.str("");
                        stream1 << datawords[link_off+link_ctr][offset].substr(0,10);
                        tmp1 << datawords[link_off+link_ctr][offset].substr(0,10);
                        stream1 << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(hwZPV.range(9,0))) << 14) << "00";
                        tmp2 << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(hwZPV.range(9,0))) << 14) << "00";
                        datawords[link_off+link_ctr][offset] = stream1.str();
                        if(debugWords) printf("datawords[%d][%d] = %s (%s + %s) <-- copied?? 0\n",link_off+link_ctr,offset,datawords[link_off+link_ctr][offset].c_str(),tmp1.str().c_str(),tmp2.str().c_str());
                    }
                }

                //std::cout<<"index="<<index<<" id"<<id<<std::endl;

                // put the data on the link number = link_ctr;
                while (index < NCLK_PER_BX*TMUX_IN) { // 8*36
                    stream1.str("");
                    tmp1.str(""); tmp2.str("");
                    stream1 << "0x";
                    if(doSimple){
                        // do simple, aligned inputs :)
                        stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id+1]);
                        stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]);
                        id+=2;
                    } else {
                        if (index%NFRAMES_APX_GEN0==0) {
                            // 1 frame = 3 x 64b words. First 8b of a frame is reserved for header. 
                            // Zero the first 32b for simplicity.
                            stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]) << "00000000";
                            tmp1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]);
                            tmp2.str("00000000 A");
                            id+=1;
                        }
                        else if (id == objDataLength[link_type]-1) { // here theres only a half-word (32b) to write? (or trailer bits?)
                            stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]) << "00000000";
                            tmp1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]);
                            tmp2.str("00000000 B");
                            id+=1;
                        }
                        else {
                            stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id+1]);
                            stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]);
                            tmp1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id+1]);
                            tmp2 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]) << " C";
                            id+=2;
                        }
                    }
                    datawords[link_off+link_ctr][offset+index] = stream1.str();
                    if(debugWords) printf("datawords[%d][%d] = %s (%s + %s) id=%d\n",link_off+link_ctr,offset+index,datawords[link_off+link_ctr][offset+index].c_str(),tmp1.str().c_str(),tmp2.str().c_str(),id);
                    index++;
                    // std::cout << stream1.str() << std::endl;
                    // std::cout << datawords[link_ctr][offset+(id/2)] << std::endl;
                    //printf("test = %i, link ctr = %i, clock = %i, offset = %i \n", test, link_ctr, offset+((id+1)/2), offset);
                    //std::cout<<"link_ctr = "<<link_ctr<<" link_off = "<<link_off<<" offset = "<<offset<<" index = "<<index;
                    //std::cout<<"  "<<datawords[link_off+link_ctr][offset+index-1]<<std::endl;
                    if (id >= objDataLength[link_type]) {
                        break;
                    }
                }
                if (id >= objDataLength[link_type]) link_ctr = link_max[link_type]-1;
            }
            
            //std::cout<<"-----------"<<std::endl;
        }
        link_off += NLINKS_PER_REG;
        // in this scheme, the first link is maximally saturated
        if (0 && link_off>= NLINKS_PER_REG*TMUX_IN/TMUX_OUT) link_off = 0;
        // in this scheme, inputs are spread across all links
        if (1 && link_off+NLINKS_PER_REG > NLINKS_APX_GEN0) link_off = 0;
        //std::cout<<"\t"<<test<<std::endl;
    
    }
    //std::cout<<std::endl;
    for (int ib = 0; ib < listLength; ib++){
        // std::cout << ib << " ";
        std::cout << "0x" << std::setfill('0') << std::setw(4) << std::hex << ib << "   " <<std::dec;
        for (int ia = 0; ia < NLINKS_APX_GEN0; ia++){
            //datawords[ia][ib] = "0x0000000000000000";
            std::cout << datawords[ia][ib] << "    ";
        }
        std::cout << std::endl;
    }

    std::cout<<"For all events: "<<std::endl;
    std::cout<<"\ttrack  = "<<n_alltracks<<std::endl;
    std::cout<<"\tcalo   = "<<n_allcalos<<std::endl;
    std::cout<<"\temcalo = "<<n_allemcalos<<std::endl;
    std::cout<<"\tmu     = "<<n_allmus<<std::endl;

    return 0;
}
