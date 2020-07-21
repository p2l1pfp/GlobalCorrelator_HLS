#include "test.h"

int main(){
    srand(42); // random seed
    //test_tanlambda_to_eta();
    test_resolution();

    //test_input_converter();
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

void test_converter(){

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


void test_tk_pack(){
    unsigned int ntest=3;

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

    for(unsigned long itest=0; itest<ntest; itest++){
        ap_uint<96> tk;
        rinv = 1./(itest+2); // pt = 2,3,4,...

        std::cout << "Before packing:" << std::endl;
        std::cout<<"rinv     "; for(int i=rinv_t    ::width-1; i>=0;i--){ std::cout << int( rinv    [i]);}std::cout << "  " << rinv    .to_double() << std::endl;
        std::cout<<"tkphi    "; for(int i=tkphi_t   ::width-1; i>=0;i--){ std::cout << int( tkphi   [i]);}std::cout << "  " << tkphi   .to_double() << std::endl;
        std::cout<<"tanlam   "; for(int i=tanlam_t  ::width-1; i>=0;i--){ std::cout << int( tanlam  [i]);}std::cout << "  " << tanlam  .to_double() << std::endl;
        std::cout<<"tkz0     "; for(int i=tkz0_t    ::width-1; i>=0;i--){ std::cout << int( tkz0    [i]);}std::cout << "  " << tkz0    .to_double() << std::endl;
        std::cout<<"tkd0     "; for(int i=tkd0_t    ::width-1; i>=0;i--){ std::cout << int( tkd0    [i]);}std::cout << "  " << tkd0    .to_double() << std::endl;
        std::cout<<"chi2rphi "; for(int i=chi2rphi_t::width-1; i>=0;i--){ std::cout << int( chi2rphi[i]);}std::cout << "  " << chi2rphi.to_double() << std::endl;
        std::cout<<"chi2rz   "; for(int i=chi2rz_t  ::width-1; i>=0;i--){ std::cout << int( chi2rz  [i]);}std::cout << "  " << chi2rz  .to_double() << std::endl;
        std::cout<<"bendChi2 "; for(int i=bendChi2_t::width-1; i>=0;i--){ std::cout << int( bendChi2[i]);}std::cout << "  " << bendChi2.to_double() << std::endl;
        std::cout<<"hit      "; for(int i=hit_t     ::width-1; i>=0;i--){ std::cout << int( hit     [i]);}std::cout << "  " << hit     .to_double() << std::endl;
        std::cout<<"trackMVA "; for(int i=trackMVA_t::width-1; i>=0;i--){ std::cout << int( trackMVA[i]);}std::cout << "  " << trackMVA.to_double() << std::endl;
        std::cout<<"extraMVA "; for(int i=extraMVA_t::width-1; i>=0;i--){ std::cout << int( extraMVA[i]);}std::cout << "  " << extraMVA.to_double() << std::endl;
        std::cout<<"valid    "; for(int i=valid_t   ::width-1; i>=0;i--){ std::cout << int( valid   [i]);}std::cout << "  " << valid   .to_double() << std::endl;

        pack_L1T_track(tk, rinv, tkphi, tanlam, tkz0, tkd0, chi2rphi, chi2rz, bendChi2, hit, trackMVA, extraMVA, valid);
        unpack_L1T_track(tk, rinv, tkphi, tanlam, tkz0, tkd0, chi2rphi, chi2rz, bendChi2, hit, trackMVA, extraMVA, valid);

        std::cout << "After packing:" << std::endl;
        std::cout<<"rinv     "; for(int i=rinv_t    ::width-1; i>=0;i--){ std::cout << int( rinv    [i]);}std::cout << "  " << rinv    .to_double() << std::endl;
        std::cout<<"tkphi    "; for(int i=tkphi_t   ::width-1; i>=0;i--){ std::cout << int( tkphi   [i]);}std::cout << "  " << tkphi   .to_double() << std::endl;
        std::cout<<"tanlam   "; for(int i=tanlam_t  ::width-1; i>=0;i--){ std::cout << int( tanlam  [i]);}std::cout << "  " << tanlam  .to_double() << std::endl;
        std::cout<<"tkz0     "; for(int i=tkz0_t    ::width-1; i>=0;i--){ std::cout << int( tkz0    [i]);}std::cout << "  " << tkz0    .to_double() << std::endl;
        std::cout<<"tkd0     "; for(int i=tkd0_t    ::width-1; i>=0;i--){ std::cout << int( tkd0    [i]);}std::cout << "  " << tkd0    .to_double() << std::endl;
        std::cout<<"chi2rphi "; for(int i=chi2rphi_t::width-1; i>=0;i--){ std::cout << int( chi2rphi[i]);}std::cout << "  " << chi2rphi.to_double() << std::endl;
        std::cout<<"chi2rz   "; for(int i=chi2rz_t  ::width-1; i>=0;i--){ std::cout << int( chi2rz  [i]);}std::cout << "  " << chi2rz  .to_double() << std::endl;
        std::cout<<"bendChi2 "; for(int i=bendChi2_t::width-1; i>=0;i--){ std::cout << int( bendChi2[i]);}std::cout << "  " << bendChi2.to_double() << std::endl;
        std::cout<<"hit      "; for(int i=hit_t     ::width-1; i>=0;i--){ std::cout << int( hit     [i]);}std::cout << "  " << hit     .to_double() << std::endl;
        std::cout<<"trackMVA "; for(int i=trackMVA_t::width-1; i>=0;i--){ std::cout << int( trackMVA[i]);}std::cout << "  " << trackMVA.to_double() << std::endl;
        std::cout<<"extraMVA "; for(int i=extraMVA_t::width-1; i>=0;i--){ std::cout << int( extraMVA[i]);}std::cout << "  " << extraMVA.to_double() << std::endl;
        std::cout<<"valid    "; for(int i=valid_t   ::width-1; i>=0;i--){ std::cout << int( valid   [i]);}std::cout << "  " << valid   .to_double() << std::endl;
        std::cout<<std::endl;

    }
}
void test_pf_pack(){
    unsigned int ntest=3;

    // PF parts for testing
    pt_t pf_pt = 20;
    pt_t pf_pterr = 22;
    etaphi_t pf_eta = -12;
    etaphi_t pf_phi = 14;
    z0_t pf_z0 = -3;
    bool pf_TightQuality = true;

    for(unsigned long itest=0; itest<ntest; itest++){
        pf_pt = 2+itest;
        ap_uint<64> tk;

        std::cout << "Before packing:" << std::endl;
        std::cout<<"pf_pt    "; for(int i=pt_t    ::width-1; i>=0;i--){ std::cout << int(pf_pt   [i]);}std::cout << "  " << pf_pt   .to_double() << std::endl;
        std::cout<<"pf_pterr "; for(int i=pt_t    ::width-1; i>=0;i--){ std::cout << int(pf_pterr[i]);}std::cout << "  " << pf_pterr.to_double() << std::endl;
        std::cout<<"pf_eta   "; for(int i=etaphi_t::width-1; i>=0;i--){ std::cout << int(pf_eta  [i]);}std::cout << "  " << pf_eta  .to_double() << std::endl;
        std::cout<<"pf_phi   "; for(int i=etaphi_t::width-1; i>=0;i--){ std::cout << int(pf_phi  [i]);}std::cout << "  " << pf_phi  .to_double() << std::endl;
        std::cout<<"pf_z0    "; for(int i=z0_t    ::width-1; i>=0;i--){ std::cout << int(pf_z0   [i]);}std::cout << "  " << pf_z0   .to_double() << std::endl; 
        std::cout<<"pf_TightQuality    " << pf_TightQuality << std::endl;
        std::cout<<"tk       "; for(int i=ap_uint<64>::width-1; i>=0;i--){ std::cout << int(tk   [i]);}std::cout << std::endl;

        pack_pf_track(tk, pf_pt, pf_pterr, pf_eta, pf_phi, pf_z0, pf_TightQuality);
        unpack_pf_track(tk, pf_pt, pf_pterr, pf_eta, pf_phi, pf_z0, pf_TightQuality);

        std::cout << "After packing:" << std::endl;
        std::cout<<"pf_pt    "; for(int i=pt_t    ::width-1; i>=0;i--){ std::cout << int(pf_pt   [i]);}std::cout << "  " << pf_pt   .to_double() << std::endl;
        std::cout<<"pf_pterr "; for(int i=pt_t    ::width-1; i>=0;i--){ std::cout << int(pf_pterr[i]);}std::cout << "  " << pf_pterr.to_double() << std::endl;
        std::cout<<"pf_eta   "; for(int i=etaphi_t::width-1; i>=0;i--){ std::cout << int(pf_eta  [i]);}std::cout << "  " << pf_eta  .to_double() << std::endl;
        std::cout<<"pf_phi   "; for(int i=etaphi_t::width-1; i>=0;i--){ std::cout << int(pf_phi  [i]);}std::cout << "  " << pf_phi  .to_double() << std::endl;
        std::cout<<"pf_z0    "; for(int i=z0_t    ::width-1; i>=0;i--){ std::cout << int(pf_z0   [i]);}std::cout << "  " << pf_z0   .to_double() << std::endl; 
        std::cout<<"pf_TightQuality    " << pf_TightQuality << std::endl;
        std::cout<<"tk       "; for(int i=ap_uint<64>::width-1; i>=0;i--){ std::cout << int(tk   [i]);}std::cout << std::endl;
        std::cout<<std::endl;
    }
}
