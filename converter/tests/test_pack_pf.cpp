//
#include "test.h"

void test_pack_pf(){
    unsigned int ntrials = 10;

    // PF parts for testing
    pt_t pf_pt = 20;
    pt_t pf_pterr = 22;
    etaphi_t pf_eta = -12;
    etaphi_t pf_phi = 14;
    z0_t pf_z0 = -3;
    bool pf_TightQuality = true;

    for(unsigned int itrial = 0; itrial< ntrials; itrial++){
        pf_pt = 2+itrial;
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
