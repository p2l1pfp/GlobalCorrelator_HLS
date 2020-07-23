//
#include "test.h"

void test_prop_tanlambda(){
    unsigned int ntrials = 10000;

    float tl, tl_calo, z0;
    tanlam_t tl_hw, tl_calo_hw;
    tkz0_t z0_hw;

    std::ofstream outfile;
    outfile.open("../../../../tests/results/test_prop_tanlambda.txt");
    outfile << "# tl tl_calo z0 tl_hw tl_calo_hw z0_hw" << endl;
    for(unsigned int itrial = 0; itrial< ntrials; itrial++){
        tl = urand(-8,8);
        z0 = urand(-30,30);
        tl_calo = propagate_tanlam_ref(z0, tl);

        tl_hw = tl;
        z0_hw = z0;
        propagate_tanlam(z0_hw, tl_hw, tl_calo_hw);
        
        outfile << tl    << " " << tl_calo    << " " << z0    << " "
                << tl_hw << " " << tl_calo_hw << " " << z0_hw << "\n";
    }
    outfile.close();
}

