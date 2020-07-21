//
#include "test.h"

void test_tanlambda_to_eta(){
    unsigned int ntrials = 100;

    float tanlambda, eta;
    tanlam_t tanlambda_hw;
    etaphi_t eta_hw;

    std::ofstream outfile;
    outfile.open("../../../../tests/results/test_tanlambda_to_eta.txt");
    outfile << "tanlambda eta tanlambda_hw eta_hw" << endl;
    for(unsigned int itrial = 0; itrial< ntrials; itrial++){
        tanlambda = urand(-8,8);
        tanlambda_hw = tanlambda;
        eta = tanlam_to_eta_ref(tanlambda);
        convert_eta(tanlambda_hw, eta_hw);
        outfile << tanlambda    << " " << eta    << " "
                << tanlambda_hw << " " << eta_hw << "\n";
    }
    outfile.close();
}

