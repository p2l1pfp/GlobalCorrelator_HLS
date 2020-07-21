//
#include "test.h"

void test_prop_phi(){
    unsigned int ntrials = 100;

    float pt_inv, dphi;
    rinv_t pt_inv_hw;
    etaphi_t dphi_hw;

    std::ofstream outfile;
    outfile.open("../../../../tests/results/test_prop_phi.txt");
    outfile << "pt dphi pt_inv dphi" << endl;
    for(unsigned int itrial = 0; itrial< ntrials; itrial++){
        pt_inv = urand(0.01,0.5);
        dphi = convert_dphi_ref( 1./pt_inv );

        pt_inv_hw = pt_inv;
        convert_dphi(pt_inv_hw, dphi_hw);
        
        outfile << pt_inv    << " " << dphi    << " "
                << pt_inv_hw << " " << dphi_hw << "\n";
    }
    outfile.close();
}


