//
#include "test.h"

void test_tanlambda_to_eta(){
    unsigned int ntrials = 10000;

    float tanlambda, eta;
    tanlam_t tanlambda_hw;
    etaphi_t eta_hw;

    std::ofstream outfile;
    outfile.open("../../../../tests/results/test_tanlambda_to_eta.txt",
                 std::ofstream::app); // append, for comparisons
    if(outfile.tellp()==0) // write header if new file
        outfile << "# tanlambda eta tanlambda_hw eta_hw eta_hwf tab_size" << endl;
    for(unsigned int itrial = 0; itrial< ntrials; itrial++){
        tanlambda = urand(-6,6); // eta=2.5
        tanlambda_hw = tanlambda;
        eta = tanlam_to_eta_ref(tanlambda);
        convert_eta(tanlambda_hw, eta_hw);
        outfile << tanlambda    << " " << eta    << " "
                << tanlambda_hw << " " << eta_hw << " "
                << eta_hw/PF_ETAPHI_SCALE << " "
                << ETA_TAB_SIZE << "\n";
    }
    outfile.close();
}

