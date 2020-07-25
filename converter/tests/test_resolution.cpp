//
#include "test.h"

void test_resolution(){
    unsigned int ntrials = 10000;

    float pt, pt_err, eta;
    pt_t pt_hw, pt_err_hw;
    etaphi_t eta_hw;

    std::ofstream outfile;
    outfile.open("../../../../tests/results/test_resolution.txt");
    outfile << "# pt pt_err eta pt_hw pt_hwf pt_err_hw pt_err_hwf eta_hw eta_hwf" << endl;
    for(unsigned int itrial = 0; itrial< ntrials; itrial++){
        eta = urand(-2.5,2.5);
        pt = 1./urand(0.01,0.5); // flat curv 2 to 100 GeV
        pt_err = reso_calo_ref(pt, eta);

        pt_hw = pt * PF_PT_SCALE;
        eta_hw = eta * PF_ETAPHI_SCALE;
        reso_calo(pt_hw, eta_hw, pt_err_hw);
        
        outfile << pt    << " " << pt_err    << " " << eta    << " "
                << pt_hw     << " " << float(pt_hw)/PF_PT_SCALE << " "
                << pt_err_hw << " " << float(pt_err_hw)/PF_PT_SCALE << " "
                << eta_hw    << " " << float(eta_hw)/PF_ETAPHI_SCALE << "\n";
    }
    outfile.close();
}

