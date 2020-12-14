#include "tmux_create_test.h"
#include <map>

#define NETA_SMALL 2
#define NPHI_SMALL 9
#define NREGIONS (NETA_SMALL*NPHI_SMALL)

int mp7DataLength = NTRACK+NCALO+NEMCALO+NMU;
int objDataLength[4] = {(NWORDS_TRACK*NTRACK), ((NWORDS_EMCALO*NEMCALO)+(NWORDS_TRACK*NTRACK)), ((NWORDS_EMCALO*NEMCALO)+(NWORDS_CALO*NCALO)+(NWORDS_TRACK*NTRACK)), ((NWORDS_EMCALO*NEMCALO)+(NWORDS_CALO*NCALO)+(NWORDS_TRACK*NTRACK)+(NWORDS_MU*NMU))};
int link_max[4] = {NLINKS_PER_TRACK, NLINKS_PER_TRACK+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO+NLINKS_PER_MU};
int link_min[4] = {0, NLINKS_PER_TRACK, NLINKS_PER_TRACK+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_EMCALO+NLINKS_PER_CALO};

//unsigned int outputOrder[TMUX_IN] = {0,2,4,6,8,10,12,14,16,15,17,1,3,5,7,9,11,13};//NPHI x NETA
//unsigned int outputOrder[TMUX_IN] = {0,11,1,12,2,13,3,14,4,15,5,16,6,17,7,9,8,10};//NPHI x NETA
unsigned int outputOrder[TMUX_IN] = {0,2,1,3,11,4,12,5,13,6,14,7,15,8,16,9,17,10};//NPHI x NETA
//  mapping from Ryan:
//  eta, phi — 0,0 --- 1,8 
//  0,0 - 0,2 - 0,1 - 0,3 - 1,2 - 0,4 - 1,3 - 0,5 - 1,4 - 0,6 - 1,5 - 0,7 - 1,6 - 0,8 - 1,7 - 1,0 - 1,8 - 1,1
//  0     2     1     3     11    4     12    5     13    6     14    7     15    8     16    9     17    10

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

    

    // make ordering for (link,clock)
    std::map<std::pair<int,int>, int> ilink_iclk_order;
    std::map<int, std::pair<int,int> > order_ilink_iclk;

    int CLKMAX = ((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_TRACK;
    int counter=0;
    for(int iclk_plus_ilink=0; iclk_plus_ilink < CLKMAX+NLINKS_PER_TRACK-1; iclk_plus_ilink++){ // clock
        for(int iclk=0;iclk<CLKMAX && iclk<=iclk_plus_ilink; iclk++){ // ilink
            //for(int il=0;il<NLINKS_PER_TRACK && il<=iclk_plus_ilink; il++){ // ilink
            int il = iclk_plus_ilink-iclk;
            if (il >= NLINKS_PER_TRACK) continue;
            // cout << il << " " << iclk << " -> " << counter << endl;

            // record correspondence
            auto p = std::make_pair(il,iclk);
            ilink_iclk_order[p]=counter;
            order_ilink_iclk[counter]=p;

            counter++;
        }
    }

    // input format: could be random or coming from simulation
    //RandomPFInputs inputs(37); // 37 is a good random number
    DiscretePFInputs inputs("barrel_sectors_1x1_TTbar_PU200.dump");
    //DiscretePFInputs inputs("dummy.dump");
    
    // input TP objects (structs for reading from dump files)
    l1tpf_int::CaloCluster calo[NCALO_TMUX]; 
    l1tpf_int::CaloCluster emcalo[NEMCALO_TMUX]; 
    l1tpf_int::PropagatedTrack track[NTRACK_TMUX]; 
    l1tpf_int::Muon mu[NMU_TMUX];
    z0_t hwZPV;

    // Data words for input to regionizer
    //   NCLK_PER_BX is the number of frames per bx (320 mhz / 40mhz)
    //   (320 mhz / 40mhz) * (Nevts * TM18 + (TM36-TM18)) assuming we alternate between two link groups
    //   Though the math is identical if instead we leave the links 1/3 empty and use all 3 groups
    const int input_listLength = NCLK_PER_BX*((NTEST*TMUX_OUT)+(TMUX_IN-TMUX_OUT));
    std::string input_datawords [NLINKS_APX_GEN0][input_listLength];
    std::string input_datawords_cvt [NLINKS_APX_GEN0][input_listLength];
    for (int ia = 0; ia < NLINKS_APX_GEN0; ia++){
        for (int ib = 0; ib < input_listLength; ib++){
            input_datawords[ia][ib] = "0x0000000000000000";
            input_datawords_cvt[ia][ib] = "0x0000000000000000";
        }
    }

    // Regionizer output data words
    std::string output_datawords[NTEST*TMUX_IN][mp7DataLength+1];
    for (int ia = 0; ia < NTEST*TMUX_IN; ia++){
        for (int ib = 0; ib < mp7DataLength+1; ib++){
            output_datawords[ia][ib] = "0000000000000000";
        }
    }

    // Layer 1 (PF+PUPPI) output data words
    const int layer1_listLength = NTEST*NREGIONS;
    std::string layer1_datawords [layer1_listLength][NLINKS_APX_GEN0];
    for (int ia = 0; ia < layer1_listLength; ia++){
        for (int ib = 0; ib < NLINKS_APX_GEN0; ib++){
            layer1_datawords[ia][ib] = "0000000000000000";
        }
    }

    // -----------------------------------------
    // Regionizer setup:
    // Determine eta and phi region boundaries

    // Phi region boundaries, including wraparound
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
        // std::cout << "TEST " << phi_bounds_lo[ip] << "  to  " << phi_bounds_hi[ip] << std::endl;
    }

    // Eta region boundaries
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
    int etalo = eta_bounds_lo.front();
    int etahi = eta_bounds_hi.back();

    // for now, we want to just consider the
    // [-243,-90) and [-154,32) regions
    etahi = 31;

    // Log region boundaries for later reference
    std::vector< std::string > str_lims;
    for (int ies = 0; ies < NETA_SMALL; ies++) {
        for (int ips = 0; ips < NPHI_SMALL; ips++) {
            char str[80];
            sprintf(str, "eta in [%d,%d) and phi in [%d,%d)", eta_bounds_lo[ies], eta_bounds_hi[ies],
                    phi_bounds_lo[ips], phi_bounds_hi[ips]);
            str_lims.push_back( str );
        }
    }
    std::ofstream outfile_srs;
    outfile_srs.open("../../../../SR_limits.txt");
    for(int i=0;i<TMUX_IN;i++)
        outfile_srs << "Small region " << i << " has " << str_lims[outputOrder[i]] << std::endl;
    outfile_srs.close();


    // -----------------------------------------
    // run multiple tests (events)

    int link_off = 0;
    int input_offset = 0;
    for (int test = 0; test < NTEST; ++test) {

        // initialize TP objects receiving input data
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
            // get the TP inputs from the 'inputs' object
            if (!inputs.nextRegion_tmux_raw(calo, emcalo, track, mu, hwZPV)) break;
        } else if(0) {
            // quickly set some fake data for testing
            track[0].hwPt  = 1 + test * 16; // one of each object per event
            calo[0].hwPt   = 2 + test * 16;
            emcalo[0].hwPt = 3 + test * 16;
            mu[0].hwPt     = 4 + test * 16;
            hwZPV          = 7 + test * 16;
        }
        
        // initialize TP object structures, distributed across input links
        std::vector<l1tpf_int::CaloCluster> calo_tp[NLINKS_PER_CALO]; 
        std::vector<l1tpf_int::CaloCluster> emcalo_tp[NLINKS_PER_EMCALO]; 
        std::vector<l1tpf_int::PropagatedTrack> track_tp[NLINKS_PER_TRACK]; 
        std::vector<l1tpf_int::Muon> mu_tp[NLINKS_PER_MU];
        l1tpf_int::PropagatedTrack track_dummy; track_dummy.hwPt = 0; track_dummy.hwEta = 0; track_dummy.hwPhi = 0; track_dummy.hwCaloPtErr = 0; track_dummy.hwZ0 = 0; track_dummy.hwStubs = 0; track_dummy.hwChi2 = 0;
        l1tpf_int::CaloCluster calo_dummy; calo_dummy.hwPt = 0; calo_dummy.hwEta = 0; calo_dummy.hwPhi = 0; calo_dummy.isEM = 0; calo_dummy.hwEmPt = 0;
        l1tpf_int::CaloCluster emcalo_dummy; emcalo_dummy.hwPt = 0; emcalo_dummy.hwEta = 0; emcalo_dummy.hwPhi = 0; emcalo_dummy.hwPtErr = 0;
        l1tpf_int::Muon mu_dummy; mu_dummy.hwPt = 0; mu_dummy.hwEta = 0; mu_dummy.hwPhi = 0;

        // analogous structures for converted input objects
        std::vector<HadCaloObj> calo_cvt[NLINKS_PER_CALO]; 
        std::vector<EmCaloObj> emcalo_cvt[NLINKS_PER_EMCALO]; 
        std::vector<TkObj> track_cvt[NLINKS_PER_TRACK]; 
        std::vector<MuObj> mu_cvt[NLINKS_PER_MU];
        HadCaloObj calo_cvt_dummy; calo_cvt_dummy.hwPt = 0; calo_cvt_dummy.hwEmPt = 0; calo_cvt_dummy.hwEta = 0; calo_cvt_dummy.hwPhi = 0; calo_cvt_dummy.hwIsEM = 0; 
        EmCaloObj emcalo_cvt_dummy; emcalo_cvt_dummy.hwPt = 0; emcalo_cvt_dummy.hwPtErr = 0;  emcalo_cvt_dummy.hwEta = 0; emcalo_cvt_dummy.hwPhi = 0;
        TkObj track_cvt_dummy; track_cvt_dummy.hwPt = 0; track_cvt_dummy.hwPtErr = 0; track_cvt_dummy.hwEta = 0; track_cvt_dummy.hwPhi = 0; track_cvt_dummy.hwZ0 = 0; 
        MuObj mu_cvt_dummy; mu_cvt_dummy.hwPt = 0; mu_cvt_dummy.hwPtErr = 0; mu_cvt_dummy.hwEta = 0; mu_cvt_dummy.hwPhi = 0;



        // -----------------------------------------
        // Fill input links, convert TP to PF objects

        // input link book-keeping
        int ilink = 0;
        int ntracks = 0;
        int ncalos = 0;
        int nemcalos = 0;
        int nmus = 0;

        ilink = NLINKS_PER_TRACK;
        for (int i = 0; i < NTRACK_TMUX; ++i) {
            if (int(track[i].hwEta) < etalo or int(track[i].hwEta) > etahi) continue;
            if (int(track[i].hwPt) == 0) continue;
            // pick input link and record
            pick_link(ilink, track[i]);
            track_tp[ilink].push_back(track[i]);
            // convert from L1Tk input format to PF format
            TkObj pf_track; 
            track_convert(track[i], pf_track, ilink);
            //cout << " - " << track[i].to_string(16) << "  " << pf_track.to_string(16) << endl;
            track_cvt[ilink].push_back(pf_track);
            ntracks++;
        }
        ilink = NLINKS_PER_CALO;
        for (int i = 0; i < NCALO_TMUX; ++i) {
            if (int(calo[i].hwEta) < etalo or int(calo[i].hwEta) > etahi) continue;
            if (int(calo[i].hwPt) == 0) continue;
            pick_link_had(ilink, calo[i]);
            calo_tp[ilink].push_back(calo[i]);
            // trivial conversion for now
            HadCaloObj pf_calo; 
            dpf2fw::convert(calo[i], pf_calo);
            calo_cvt[ilink].push_back(pf_calo);
            ncalos++;
        }
        ilink = NLINKS_PER_EMCALO;
        for (int i = 0; i < NEMCALO_TMUX; ++i) {
            if (int(emcalo[i].hwEta) < etalo or int(emcalo[i].hwEta) > etahi) continue;
            if (int(emcalo[i].hwPt) == 0) continue;
            pick_link_em(ilink, emcalo[i]);
            emcalo_tp[ilink].push_back(emcalo[i]);
            // trivial conversion for now
            EmCaloObj pf_emcalo; 
            dpf2fw::convert(emcalo[i], pf_emcalo);
            emcalo_cvt[ilink].push_back(pf_emcalo);
            nemcalos++;
        }
        ilink = NLINKS_PER_MU;
        for (int i = 0; i < NMU_TMUX; ++i) {
            if (int(mu[i].hwEta) < etalo or int(mu[i].hwEta) > etahi) continue;
            if (int(mu[i].hwPt) == 0) continue;
            pick_link(ilink, mu[i]);
            mu_tp[ilink].push_back(mu[i]);
            // trivial conversion for now
            MuObj pf_mu; 
            dpf2fw::convert(mu[i], pf_mu);
            mu_cvt[ilink].push_back(pf_mu);
            nmus++;
        }

        // resize to ensure number of objects can be sent in one link group
        for (int il = 0; il < NLINKS_PER_TRACK; il++) {
            if (track_tp[il].size() > ((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_TRACK)
                cout << "Event " << test << ": truncating input tracks on link " << il << endl;
            track_tp [il].resize(((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_TRACK, track_dummy);
            track_cvt[il].resize(((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_TRACK, track_cvt_dummy);
        }
        for (int il = 0; il < NLINKS_PER_CALO; il++) {
            if (calo_tp[il].size() > ((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_CALO)
                cout << "Event " << test << ": truncating input hadcalos on link " << il << endl;
            calo_tp [il].resize(((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_CALO, calo_dummy);
            calo_cvt[il].resize(((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_CALO, calo_cvt_dummy);
        }
        for (int il = 0; il < NLINKS_PER_EMCALO; il++) {
            if (emcalo_tp[il].size() > ((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_EMCALO)
                cout << "Event " << test << ": truncating input emcalos on link " << il << endl;
            emcalo_tp [il].resize(((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_EMCALO, emcalo_dummy);
            emcalo_cvt[il].resize(((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_EMCALO, emcalo_cvt_dummy);
        }
        for (int il = 0; il < NLINKS_PER_MU; il++) {
            if (mu_tp[il].size() > ((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_MU)
                cout << "Event " << test << ": truncating input muons on link " << il << endl;
            mu_tp [il].resize(((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_MU, mu_dummy);
            mu_cvt[il].resize(((NCLK_PER_BX*TMUX_IN-1)*2)/NWORDS_MU, mu_cvt_dummy);
        }

        // -----------------------------------------
        // Write input TP objects to the input links
        unsigned int link_type = 0; //0=track, 1=emcalo, 2=calo, 3=mu
        unsigned int obj_link_no;
        input_offset = NCLK_PER_BX*TMUX_OUT*test; // 8 * 18 * nevt
        for (int link_ctr = 0; link_ctr < NLINKS_PER_REG; link_ctr++) {

            if      (link_ctr < link_max[0]) link_type = 0;
            else if (link_ctr < link_max[1]) link_type = 1;
            else if (link_ctr < link_max[2]) link_type = 2;
            else if (link_ctr < link_max[3]) link_type = 3;
            obj_link_no = link_ctr-link_min[link_type];

            
            if (link_type == 0) {                
                write_track_vector_to_link(track_tp[obj_link_no], input_datawords[link_off+link_ctr], input_offset);
                write_track_vector_to_link(track_cvt[obj_link_no], input_datawords_cvt[link_off+link_ctr], input_offset);
            } else if (link_type == 1) {
                write_emcalo_vector_to_link(emcalo_tp[obj_link_no], input_datawords[link_off+link_ctr], input_offset); 
                write_emcalo_vector_to_link(emcalo_tp[obj_link_no], input_datawords_cvt[link_off+link_ctr], input_offset);                
            } else if (link_type == 2) {
                write_calo_vector_to_link(calo_tp[obj_link_no], input_datawords[link_off+link_ctr], input_offset);
                write_calo_vector_to_link(calo_tp[obj_link_no], input_datawords_cvt[link_off+link_ctr], input_offset);
            } else if (link_type == 3) {
                write_mu_vector_to_link(mu_tp[obj_link_no], input_datawords[link_off+link_ctr], input_offset);
                write_mu_vector_to_link(mu_tp[obj_link_no], input_datawords_cvt[link_off+link_ctr], input_offset);
            }
        }

        // turn off vtx for now
        if (false){
            std::stringstream stream2;
            stream2.str("");
            stream2 << "0x00000000" << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(hwZPV.range(9,0))) << 14) << "00"; 
            input_datawords[link_off+NLINKS_PER_REG][input_offset] = stream2.str();
            input_datawords_cvt[link_off+NLINKS_PER_REG][input_offset] = stream2.str();
        }
        
        link_off += NLINKS_PER_REG+1;
        // 2 schemes should be equivalent in 18 / 6 setup
        // in this scheme, the first link is maximally saturated
        if (0 && link_off>= (NLINKS_PER_REG+1)*TMUX_IN/TMUX_OUT) link_off = 0;
        // in this scheme, inputs are spread across all links
        if (1 && link_off+(NLINKS_PER_REG+1) > NLINKS_APX_GEN0) link_off = 0;


        // -----------------------------------------
        // Begin regionizer: distribute to small regions

        // regional sorting book-keeping
        int i_temp[TMUX_IN] = {0};
        int ireg = 0;
        int ntracks_smallreg[TMUX_IN] = {0};
        int ncalos_smallreg[TMUX_IN] = {0};
        int nemcalos_smallreg[TMUX_IN] = {0};
        int nmus_smallreg[TMUX_IN] = {0};
        int ntracks_reg = 0;
        int ncalos_reg = 0;
        int nemcalos_reg = 0;
        int nmus_reg = 0;

        // initialize converted (PF input) objects, sorted into regions
        HadCaloObj calo_pf_in[TMUX_IN][NCALO]; EmCaloObj emcalo_pf_in[TMUX_IN][NEMCALO]; TkObj track_pf_in[TMUX_IN][NTRACK]; MuObj mu_pf_in[TMUX_IN][NMU];
        for (int ir = 0; ir < TMUX_IN; ir++) {
            for (int i = 0; i < NTRACK; ++i) {
                track_pf_in[ir][i].hwPt = 0; track_pf_in[ir][i].hwPtErr = 0; track_pf_in[ir][i].hwEta = 0; track_pf_in[ir][i].hwPhi = 0; track_pf_in[ir][i].hwZ0 = 0; 
            }
            for (int i = 0; i < NCALO; ++i) {
                calo_pf_in[ir][i].hwPt = 0; calo_pf_in[ir][i].hwEmPt = 0; calo_pf_in[ir][i].hwEta = 0; calo_pf_in[ir][i].hwPhi = 0; calo_pf_in[ir][i].hwIsEM = 0; 
            }
            for (int i = 0; i < NEMCALO; ++i) {
                emcalo_pf_in[ir][i].hwPt = 0; emcalo_pf_in[ir][i].hwPtErr = 0;  emcalo_pf_in[ir][i].hwEta = 0; emcalo_pf_in[ir][i].hwPhi = 0;
            }
            for (int i = 0; i < NMU; ++i) {
                mu_pf_in[ir][i].hwPt = 0; mu_pf_in[ir][i].hwPtErr = 0; mu_pf_in[ir][i].hwEta = 0; mu_pf_in[ir][i].hwPhi = 0;
            }
        }

        // fill regions, starting from beginning of each link
        std::fill(i_temp, i_temp+TMUX_IN, 0);
        for(int i=0; i<CLKMAX*NLINKS_PER_TRACK;i++){
            auto p = order_ilink_iclk[i];
            int il = p.first;
            int ind = p.second;

            auto pf_track = track_cvt[il].at(ind);
            if (!(pf_track.hwPt>0)) continue;
            // find region
            for (int ies = 0; ies < NETA_SMALL; ies++) {
                if (int(pf_track.hwEta) >= eta_bounds_lo[ies] and int(pf_track.hwEta) < eta_bounds_hi[ies]) {
                    for (int ips = 0; ips < NPHI_SMALL; ips++) {
                        if ( isInPhiRegion(pf_track.hwPhi, phi_bounds_lo[ips], phi_bounds_hi[ips]) ) { 
                            if (i_temp[ies*NPHI_SMALL+ips]==NTRACK){
                                continue;
                            }
                            track_pf_in[ies*NPHI_SMALL+ips][i_temp[ies*NPHI_SMALL+ips]] = pf_track;
                            i_temp[ies*NPHI_SMALL+ips] += 1;
                            ntracks_smallreg[ies*NPHI_SMALL+ips]++;
                            ntracks_reg++;
                        }
                    }
                }
            }
        }
        std::fill(i_temp, i_temp+TMUX_IN, 0);
        for (int ind = 0; ind < calo_tp[0].size(); ind++) { // works since same size for all links
            for (int il = 0; il < NLINKS_PER_CALO; il++) {
                auto pf_calo = calo_cvt[il].at(ind);
                // find region
                for (int ies = 0; ies < NETA_SMALL; ies++) {
                    if (int(pf_calo.hwEta) >= eta_bounds_lo[ies] and int(pf_calo.hwEta) < eta_bounds_hi[ies]) {
                        for (int ips = 0; ips < NPHI_SMALL; ips++) {
                            if ( isInPhiRegion(pf_calo.hwPhi, phi_bounds_lo[ips], phi_bounds_hi[ips]) ) { 
                                if (i_temp[ies*NPHI_SMALL+ips]==NCALO) continue;
                                calo_pf_in[ies*NPHI_SMALL+ips][i_temp[ies*NPHI_SMALL+ips]] = pf_calo;
                                i_temp[ies*NPHI_SMALL+ips] += 1;
                                ncalos_smallreg[ies*NPHI_SMALL+ips]++;
                                ncalos_reg++;
                            }
                        }
                    }
                }
            }
        }
        std::fill(i_temp, i_temp+TMUX_IN, 0);
        for (int ind = 0; ind < emcalo_tp[0].size(); ind++) { // works since same size for all links
            for (int il = 0; il < NLINKS_PER_EMCALO; il++) {
                auto pf_emcalo = emcalo_cvt[il].at(ind);
                // find region
                for (int ies = 0; ies < NETA_SMALL; ies++) {
                    if (int(pf_emcalo.hwEta) >= eta_bounds_lo[ies] and int(pf_emcalo.hwEta) < eta_bounds_hi[ies]) {
                        for (int ips = 0; ips < NPHI_SMALL; ips++) {
                            if ( isInPhiRegion(pf_emcalo.hwPhi, phi_bounds_lo[ips], phi_bounds_hi[ips]) ) { 
                                if (i_temp[ies*NPHI_SMALL+ips]==NEMCALO) continue;
                                emcalo_pf_in[ies*NPHI_SMALL+ips][i_temp[ies*NPHI_SMALL+ips]] = pf_emcalo;
                                i_temp[ies*NPHI_SMALL+ips] += 1;
                                nemcalos_smallreg[ies*NPHI_SMALL+ips]++;
                                nemcalos_reg++;
                            }
                        }
                    }
                }
            }
        }
        std::fill(i_temp, i_temp+TMUX_IN, 0);
        for (int ind = 0; ind < mu_tp[0].size(); ind++) { // works since same size for all links
            for (int il = 0; il < NLINKS_PER_MU; il++) {
                auto pf_mu = mu_cvt[il].at(ind);
                // find region
                for (int ies = 0; ies < NETA_SMALL; ies++) {
                    if (int(pf_mu.hwEta) >= eta_bounds_lo[ies] and int(pf_mu.hwEta) < eta_bounds_hi[ies]) {
                        for (int ips = 0; ips < NPHI_SMALL; ips++) {
                            if ( isInPhiRegion(pf_mu.hwPhi, phi_bounds_lo[ips], phi_bounds_hi[ips]) ) { 
                                if (i_temp[ies*NPHI_SMALL+ips]==NMU) continue;
                                mu_pf_in[ies*NPHI_SMALL+ips][i_temp[ies*NPHI_SMALL+ips]] = pf_mu;
                                i_temp[ies*NPHI_SMALL+ips] += 1;
                                nmus_smallreg[ies*NPHI_SMALL+ips]++;
                                nmus_reg++;
                            }
                        }
                    }
                }
            }
        }


        // -----------------------------------------
        // Run PF+PUPPI in small regions

        for (int ir = 0; ir < TMUX_IN; ir++) {
            MP7DataWord data_in_reg[MP7_NCHANN];
            MP7DataWord data_in[MP7_NCHANN];
            MP7DataWord data_out[MP7_NCHANN];
            for (unsigned int idm = 0; idm < MP7_NCHANN; idm++) {
                data_in_reg[idm] = 0;
                data_in[idm] = 0;
                data_out[idm] = 0;
            }
            
            // write the output of the regionizer
            mp7wrapped_pack_in_reorder(emcalo_pf_in[ir], calo_pf_in[ir], track_pf_in[ir], mu_pf_in[ir], data_in_reg);
            for (int id = 0; id < mp7DataLength; id++) {
                std::stringstream stream1;
                stream1 << std::uppercase << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in_reg[id*2+1]);
                stream1 << std::uppercase << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in_reg[id*2+0]);
                output_datawords[test*TMUX_IN+ir][id] = stream1.str();
            }
            std::stringstream stream1;
            stream1 << "00000000";
            stream1 << std::uppercase << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(hwZPV.range(9,0))) << 14) << "00";
            output_datawords[test*TMUX_IN+ir][mp7DataLength] = stream1.str();

            // pack input to PF+PUPPI block (only diff from above due to reordering)
            mp7wrapped_pack_in(emcalo_pf_in[ir], calo_pf_in[ir], track_pf_in[ir], mu_pf_in[ir], data_in);
            mp7wrapped_pfalgo3_full(data_in, data_out, hwZPV);

            for (int id = 0; id < mp7DataLength; id++) {
                std::stringstream stream1;
                stream1 << std::uppercase << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_out[id*2+1]);
                stream1 << std::uppercase << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_out[id*2+0]);
                layer1_datawords[test*NREGIONS+ir][id] = stream1.str();
            }
        }


    }

    // -----------------------------------------
    // write the emulation results

    std::ofstream outfile_inputs;
    outfile_inputs.open("../../../../inputs.txt");
    std::ofstream outfile_inputs_cvt;
    outfile_inputs_cvt.open("../../../../inputs_converted.txt");
    for (int ib = 0; ib < input_listLength; ib++){
        outfile_inputs << "0x" << std::setfill('0') << std::setw(4) << std::hex << ib << "   " <<std::dec;
        outfile_inputs_cvt << "0x" << std::setfill('0') << std::setw(4) << std::hex << ib << "   " <<std::dec;
        //for (int ia = 0; ia < NLINKS_APX_GEN0; ia++){
        for (int ia = NLINKS_APX_GEN0-1; ia >=0; ia--){ // write backwards
            outfile_inputs << input_datawords[ia][ib] << "    ";
            outfile_inputs_cvt << input_datawords_cvt[ia][ib] << "    ";
        }
        outfile_inputs << std::endl;
        outfile_inputs_cvt << std::endl;
    }
    outfile_inputs.close();
    outfile_inputs_cvt.close();

    std::ofstream outfile_outputs;
    outfile_outputs.open("../../../../output.txt");
    int iclk = 0;
    for (int ia = 0; ia < NTEST; ia++){
        for (int io = 0; io < TMUX_IN; io++){
            outfile_outputs << "0x" << std::setfill('0') << std::setw(4) << std::hex << iclk << "   " <<std::dec;
            for (int ib = 0; ib < mp7DataLength+1; ib++){
                outfile_outputs << output_datawords[ia*TMUX_IN+outputOrder[io]][ib] << "    ";
            }
            outfile_outputs << std::endl;
            iclk++;
        }
    }
    outfile_outputs.close();

    std::ofstream outfile_layer1;
    outfile_layer1.open("../../../../layer1.txt");
    iclk = 0;
    for (int ia = 0; ia < NTEST; ia++){
        for (int io = 0; io < NREGIONS; io++){
            outfile_layer1 << "0x" << std::setfill('0') << std::setw(4) << std::hex << iclk << "   " <<std::dec;
            //for (int ib = 0; ib < mp7DataLength+1; ib++){
            for (int ib = NLINKS_APX_GEN0-1; ib >= 0; ib--){
                //std::cout << ia*NREGIONS+outputOrder[io] << "-" << ib << "    ";
                outfile_layer1 << layer1_datawords[ia*NREGIONS+outputOrder[io]][ib] << "    ";
            }
            outfile_layer1 << std::endl;
            iclk++;
        }
    }
    outfile_layer1.close();

    return 0;
}
