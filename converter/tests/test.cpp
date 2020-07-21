#include "test.h"

int main(){
    // set random seed
    srand(42); 

    // run piecewise tests
    // test_tanlambda_to_eta();
    // test_resolution();
    // test_prop_tanlambda();
    // test_prop_phi();
    // test_pt_inversion();
    // test_pack_l1tk();
    // test_pack_pf();

    // run full test
    test_input_converter();
    return 0;
}


l1tk_word_t get_random_l1track(unsigned int ntest){
    static int itest=0;
    // L1T tk parts for testing
    rinv_t     rinv     = 0.5;
    tkphi_t    tkphi    = 0.2;
    tanlam_t   tanlam   = 1.;
    tkz0_t     tkz0     = -2.25;
    tkd0_t     tkd0     = 1.5;
    chi2rphi_t chi2rphi = 1;
    chi2rz_t   chi2rz   = 1;
    bendChi2_t bendChi2 = 1;
    hit_t      hit      = 0;
    trackMVA_t trackMVA = 1;
    extraMVA_t extraMVA = 1;
    valid_t    valid    = 1;

    l1tk_word_t in_tk;
    rinv = 1./(itest+2); // pt = 2,3,4,...
    tanlam = 7./(itest+1) * (itest%2 ? 1:-1); // 7, -3.5, 2.3, ...

    pack_L1T_track(in_tk, rinv, tkphi, tanlam, tkz0, tkd0, 
                   chi2rphi, chi2rz, bendChi2, hit, 
                   trackMVA, extraMVA, valid);
    itest++;
    return in_tk;
}

void test_input_converter(){

    unsigned int ntest=5;

    // L1T tk parts for testing
    // rinv_t     rinv     = 0.5;
    // tkphi_t    tkphi    = 0.2;
    // tanlam_t   tanlam   = 1.;
    // tkz0_t     tkz0     = -2.25;
    // tkd0_t     tkd0     = 1.5;
    // chi2rphi_t chi2rphi = 1;
    // chi2rz_t   chi2rz   = 1;
    // bendChi2_t bendChi2 = 1;
    // hit_t      hit      = 0;
    // trackMVA_t trackMVA = 1;
    // extraMVA_t extraMVA = 1;
    // valid_t    valid    = 1;

    // PF tk parts for testing
    pt_t pf_pt = 20;
    pt_t pf_pterr = 22;
    etaphi_t pf_eta = -12;
    etaphi_t pf_phi = 14;
    z0_t pf_z0 = -3;
    bool pf_TightQuality = true;

    for(unsigned long itest=0; itest<ntest; itest++){
        l1tk_word_t in_tk = get_random_l1track(ntest);
        pf_tk_word_t out_tk(0);

        //pack_L1T_track(in_tk, rinv, tkphi, tanlam, tkz0, tkd0, chi2rphi, chi2rz, bendChi2, hit, trackMVA, extraMVA, valid);
        pf_input_track_conv_hw(in_tk, out_tk, 0);
        unpack_pf_track(out_tk, pf_pt, pf_pterr, pf_eta, pf_phi, pf_z0, pf_TightQuality);

        // std::cout << "pT: TB (" << rinv.to_double() << ") versus HW (" << pf_pt.to_double << ")" << std::endl;
        // std::cout << "pT: TB (" << 1./rinv.to_double() << ") versus HW (" << pf_pt.to_double()/PF_PT_SCALE << ")" << std::endl;
        // std::cout << "eta: TB (" << tanlam.to_double() << ") versus HW (" << pf_eta.to_double() << ")" << std::endl;
        // std::cout << "eta: TB (" << tanlam_to_eta(tanlam.to_double()) << ") versus HW (" << pf_eta.to_double()/PF_ETAPHI_SCALE << ")" << std::endl;
    }
    //
}

