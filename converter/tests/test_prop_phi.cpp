//
#include "test.h"

void test_prop_phi(){
    unsigned int ntrials = 10000;

    float pt_inv, dphi;
    rinv_t pt_inv_hw;
    etaphi_t dphi_hw;

    std::ofstream outfile;
    outfile.open("../../../../tests/results/test_prop_phi.txt",
                 std::ofstream::app); // append, for comparisons
    if(outfile.tellp()==0) // write header if new file
        outfile << "# pt_inv dphi pt_inv_hw pt_inv_hwf dphi_hw dphi_hwf tab_size" << endl;
    for(unsigned int itrial = 0; itrial< ntrials; itrial++){
        pt_inv = urand(0.01,0.5);
        dphi = convert_dphi_ref( 1./pt_inv );

        pt_inv_hw = pt_inv;
        convert_dphi(pt_inv_hw, dphi_hw);
        
        outfile << pt_inv    << " " << dphi    << " "
                << pt_inv_hw << " " << float(pt_inv_hw)/PF_PT_SCALE << " "
                << dphi_hw   << " " << float(dphi_hw)/PF_ETAPHI_SCALE << " "
                << DPHI_TAB_SIZE << "\n";

    }

    outfile.close();
}


