#ifndef PUPPI_PuppiChecker_h
#define PUPPI_PuppiChecker_h
#include <cstdio>
#include <cmath>

class PuppiChecker {
    public:
        PuppiChecker() : 
            npt2_(0), nok_(0), n1bit_(0), nalmostok_(0), nbad_(0),
            sumDiff_(0), sumAbsDiff_(0) {}

        template<typename T, unsigned int N>
        void checkIntVsFloat(const T input[N], const PFNeutralObj puppi[N], const PFNeutralObj puppi_flt[N], bool verbose) ;

        template<unsigned int N>
        bool check(const PFNeutralObj puppi[N], const PFNeutralObj puppi_ref[N], const PFNeutralObj puppi_flt[N]) ;

        void printIntVsFloatReport() {
            int nmiss = n1bit_ + nalmostok_ + nbad_, nall = nok_ + nmiss;
            float fall = 1.0/std::max(nall,1), fmiss = 1.0/std::max(nmiss,1);
            printf("  - all        : %6d\n", nall);
            printf("  -   filled   : %6d  (%6.2f%% )   [ pT >= 2 GeV ]\n", npt2_,  npt2_ * 100.0 * fall);
            printf("  - exact match: %6d  (%6.2f%% )\n", nok_,   nok_ * 100.0 * fall);
            printf("  - mismatch   : %6d  (%6.2f%% )\n", nmiss, nmiss * 100.0 * fall);
            printf("  -   by 1*LSB : %6d  (%6.2f%% )   [ 1 unit, %.2f GeV ]\n", n1bit_, n1bit_ * 100.0 * fall, LINPUPPI_ptLSB);
            printf("  -      small : %6d  (%6.2f%% )   [ %.2f < delta(pt) <= 1 GeV + 1% ]\n", nalmostok_, nalmostok_ * 100.0 * fall, LINPUPPI_ptLSB);
            printf("  -      big   : %6d  (%6.2f%% )   [ delta(pt) > 1 GeV + 1% ]\n", (nbad_), (nbad_) * 100.0 * fall);
            printf("  - average pT  diff   %+8.4f  (on all)    %+8.4f  (on mismatch)\n", sumDiff_*fall, sumDiff_*fmiss);
            printf("  - average pT |diff|  % 8.4f  (on all)    % 8.4f  (on mismatch)\n", sumAbsDiff_*fall, sumAbsDiff_*fmiss);
        }

    private:
        int npt2_, nok_, n1bit_, nalmostok_, nbad_;
        float sumDiff_, sumAbsDiff_;
};

template<typename T, unsigned int N>
void PuppiChecker::checkIntVsFloat(const T input[N], const PFNeutralObj puppi[N], const PFNeutralObj puppi_flt[N], bool verbose) {
    for (int i = 0; i < N; ++i){
        if (input[i].hwPt > 0) {
            if (puppi_flt[i].hwPtPuppi*LINPUPPI_ptLSB >= 2) npt2_++;

            int hwPtDiff = (puppi_flt[i].hwPtPuppi - puppi[i].hwPtPuppi);
            float ptDiff = hwPtDiff * LINPUPPI_ptLSB;

            if (hwPtDiff == 0) {
                nok_++; 
            } else if (std::abs(hwPtDiff) == 1) {
                n1bit_++;
            } else if (std::abs(ptDiff)< 1 + 0.01 * puppi_flt[i].hwPtPuppi*LINPUPPI_ptLSB) {
                nalmostok_++;
            } else {
                nbad_++;
            }
            sumDiff_ += ptDiff;
            sumAbsDiff_ += std::abs(ptDiff);
            if (verbose)  printf("particle %02d pT %7.2f :  puppiPt_int %7.2f   puppiPt_flt %7.2f    diff %+7.2f\n",
                    i, input[i].hwPt * LINPUPPI_ptLSB, 
                    puppi[i].hwPtPuppi * LINPUPPI_ptLSB, puppi_flt[i].hwPtPuppi * LINPUPPI_ptLSB, ptDiff);
        }
    }
}


template<unsigned int N>
bool PuppiChecker::check(const PFNeutralObj puppi[N], const PFNeutralObj puppi_ref[N], const PFNeutralObj puppi_flt[N]) {
    bool ret = true;
    for (int i = 0; i < N; ++i){
        if (!puppi_equals(puppi_ref[i], puppi[i], "Puppi", i)) {
            ret = false;
        }
    }
    if (!ret) {
        for (int i = 0; i < N; ++i){
            printf("particle %02d:  puppiPt_hw %7.2f eta %+7d phi %+7d    puppiPt_ref %7.2f eta %+7d phi %+7d   puppiPt_flt %7.2f eta %+7d phi %+7d\n", i,
                    puppi[i].hwPtPuppi     * LINPUPPI_ptLSB, int(puppi[i].hwEta), int(puppi[i].hwPhi),
                    puppi_ref[i].hwPtPuppi * LINPUPPI_ptLSB, int(puppi[i].hwEta), int(puppi[i].hwPhi), 
                    puppi_flt[i].hwPtPuppi * LINPUPPI_ptLSB, int(puppi[i].hwEta), int(puppi[i].hwPhi));
        }
    }
    return ret;
}


#endif
