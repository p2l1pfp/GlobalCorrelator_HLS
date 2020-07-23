#include "test.h"

int main(){
    // set random seed
    srand(42); 

    // run piecewise tests
    test_tanlambda_to_eta();
    test_resolution();
    test_prop_tanlambda();
    test_prop_phi();
    test_pt_inversion();
    // test_pack_l1tk(); // these print to cout, disable for now...
    // test_pack_pf();

    // run full test
    test_input_converter();
    
    return 0;
}

void get_random_l1track(unsigned int ntest,
                        float &rinv       ,
                        float &tkphi      ,
                        float &tanlam     ,
                        float &tkz0       ,
                        float &tkd0       ,
                        float &chi2rphi   ,
                        float &chi2rz     ,
                        float &bendChi2   ,
                        int   &hit        ,
                        int   &trackMVA   ,
                        int   &extraMVA   ,
                        int   &valid){
    // (not very random so far...)
    static int itest=0;
    rinv     = 0.5;
    tkphi    = 0.2;
    tanlam   = 1.;
    tkz0     = -2.25;
    tkd0     = 1.5;
    chi2rphi = 1;
    chi2rz   = 1;
    bendChi2 = 1;
    hit      = 0;
    trackMVA = 1;
    extraMVA = 1;
    valid    = 1;
    rinv = 1./(itest+2); // pt = 2,3,4,...
    tanlam = 7./(itest+1) * (itest%2 ? 1:-1); // 7, -3.5, 2.3, ...
    itest++;
}

void test_input_converter(){
    unsigned int ntrials = 100;

    // reference values
    float rinv           ;
    float tkphi          ;
    float tanlam         ;
    float tkz0           ;
    float tkd0           ;
    float chi2rphi       ;
    float chi2rz         ;
    float bendChi2       ;
    int   hit            ;
    int   trackMVA       ;
    int   extraMVA       ;
    int   valid          ;
    float pf_pt          ;
    float pf_pterr       ;
    float pf_eta         ;
    float pf_phi         ;
    float pf_z0          ; 
    bool  pf_TightQuality;
    // hw tks
    rinv_t     hw_rinv           ;
    tkphi_t    hw_tkphi          ;
    tanlam_t   hw_tanlam         ;
    tkz0_t     hw_tkz0           ;
    tkd0_t     hw_tkd0           ;
    chi2rphi_t hw_chi2rphi       ;
    chi2rz_t   hw_chi2rz         ;
    bendChi2_t hw_bendChi2       ;
    hit_t      hw_hit            ;
    trackMVA_t hw_trackMVA       ;
    extraMVA_t hw_extraMVA       ;
    valid_t    hw_valid          ;
    pt_t       hw_pf_pt          ;
    pt_t       hw_pf_pterr       ;
    etaphi_t   hw_pf_eta         ;
    etaphi_t   hw_pf_phi         ;
    z0_t       hw_pf_z0          ;
    bool       hw_pf_TightQuality;

    std::ofstream outfile;
    outfile.open("../../../../tests/results/test.txt");
    outfile << "rinv tkphi tanlam tkz0 tkd0 chi2rphi chi2rz bendChi2 hit trackMVA extraMVA valid pf_pt pf_pterr pf_eta pf_phi pf_z0 pf_TightQuality hw_rinv hw_tkphi hw_tanlam hw_tkz0 hw_tkd0 hw_chi2rphi hw_chi2rz hw_bendChi2 hw_hit hw_trackMVA hw_extraMVA hw_valid hw_pf_pt hw_pf_pterr hw_pf_eta hw_pf_phi hw_pf_z0 hw_pf_TightQuality" << endl;

    for(unsigned int itrial = 0; itrial< ntrials; itrial++){
        // get a random track
        get_random_l1track(ntrials, rinv, tkphi, tanlam, 
                           tkz0, tkd0, chi2rphi, chi2rz, bendChi2, 
                           hit, trackMVA, extraMVA, valid);
        // reference conversion
        pf_input_track_conv_ref(rinv, tkphi, tanlam, tkz0, tkd0, 
                                chi2rphi, chi2rz, bendChi2, hit, 
                                trackMVA, extraMVA, valid, pf_pt, 
                                pf_pterr, pf_eta, pf_phi, 
                                pf_z0, pf_TightQuality);

        // get hw inputs from random track and pack
        hw_rinv     = rinv    ; 
        hw_tkphi    = tkphi   ; 
        hw_tanlam   = tanlam  ; 
        hw_tkz0     = tkz0    ; 
        hw_tkd0     = tkd0    ; 
        hw_chi2rphi = chi2rphi; 
        hw_chi2rz   = chi2rz  ; 
        hw_bendChi2 = bendChi2; 
        hw_hit      = hit     ; 
        hw_trackMVA = trackMVA; 
        hw_extraMVA = extraMVA; 
        hw_valid    = valid   ; 
    
        l1tk_word_t in_tk(0);
        pf_tk_word_t out_tk(0);
        pack_L1T_track(in_tk, hw_rinv, hw_tkphi, hw_tanlam, hw_tkz0, 
                       hw_tkd0, hw_chi2rphi, hw_chi2rz, hw_bendChi2, 
                       hw_hit, hw_trackMVA, hw_extraMVA, hw_valid);
        pf_input_track_conv_hw(in_tk, out_tk, 0/*linkNo*/);
        unpack_pf_track(out_tk, hw_pf_pt, hw_pf_pterr, hw_pf_eta, hw_pf_phi, 
                        hw_pf_z0, hw_pf_TightQuality);

        outfile << rinv               << " "
                << tkphi              << " "
                << tanlam             << " "
                << tkz0               << " "
                << tkd0               << " "
                << chi2rphi           << " "
                << chi2rz             << " "
                << bendChi2           << " "
                << hit                << " "
                << trackMVA           << " "
                << extraMVA           << " "
                << valid              << " "
                << pf_pt              << " "
                << pf_pterr           << " "
                << pf_eta             << " "
                << pf_phi             << " "
                << pf_z0              << " "
                << pf_TightQuality    << " "
                << hw_rinv            << " "
                << hw_tkphi           << " "
                << hw_tanlam          << " "
                << hw_tkz0            << " "
                << hw_tkd0            << " "
                << hw_chi2rphi        << " "
                << hw_chi2rz          << " "
                << hw_bendChi2        << " "
                << hw_hit             << " "
                << hw_trackMVA        << " "
                << hw_extraMVA        << " "
                << hw_valid           << " "
                << hw_pf_pt           << " "
                << hw_pf_pterr        << " "
                << hw_pf_eta          << " "
                << hw_pf_phi          << " "
                << hw_pf_z0           << " "
                << hw_pf_TightQuality << "\n";
    }
    outfile.close();
}
