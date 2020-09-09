//
#include "test.h"

void test_pt_inversion(){
    unsigned int ntrials = 10000;

    float pt;
    rinv_t pt_inv_hw;
    pt_t pt_hw;

    std::ofstream outfile;
    outfile.open("../../../../tests/results/test_pt_inversion.txt",
                 std::ofstream::app); // append, for comparisons
    if(outfile.tellp()==0) // write header if new file
        outfile << "# pt pt_inv pt_hw pt_hwf pt_inv_hw tab_size max_bits" << endl;
    for(unsigned int itrial = 0; itrial< ntrials; itrial++){
        pt = 1. / urand(1/140.,0.5);

        pt_inv_hw = 1./pt;
        convert_pt(pt_inv_hw, pt_hw);
        
        outfile << pt    << " " << 1./pt    << " "
                << pt_hw << " " << float(pt_hw)/(PF_PT_SCALE) << " "
                << pt_inv_hw << " "
                << PT_INV_TAB_SIZE << " "
                << PT_INV_MAX_BITS << "\n";
    }
    outfile.close();
}


