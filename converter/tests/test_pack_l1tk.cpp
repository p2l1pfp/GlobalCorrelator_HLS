//
#include "test.h"

void test_pack_l1tk(){
    unsigned int ntrials = 10;

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

    for(unsigned int itrial = 0; itrial< ntrials; itrial++){
        ap_uint<96> tk;
        rinv = 1./(itrial+2); // pt = 2,3,4,...

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
