#include "firmware/linpuppi.h"
#include <cmath>
#include <algorithm>


template<typename T, int NIn, int NOut>
void puppisort_and_crop_ref(const T in[NIn], int minHwPtPuppi, T out[NOut]) {
    T tmp[NOut];

    for (int iout = 0; iout < NOut; ++iout) {
        clear(tmp[iout]);
    }

    for (int it = 0; it < NIn; ++it) {
        if (in[it].hwPtPuppi < minHwPtPuppi) continue;
        for (int iout = NOut-1; iout >= 0; --iout) {
            if (tmp[iout].hwPtPuppi <= in[it].hwPtPuppi) {
                if (iout == 0 || tmp[iout-1].hwPtPuppi > in[it].hwPtPuppi) {
                    tmp[iout] = in[it];
                } else {
                    tmp[iout] = tmp[iout-1];
                }
            }
        }
    }

    for (int iout = 0; iout < NOut; ++iout) {
        out[iout] = tmp[iout];
    }

}


void linpuppi_chs_ref(z0_t pvZ0, const PFChargedObj pfch[NTRACK], PFChargedObj outallch[NTRACK]) {
    for (unsigned int i = 0; i < NTRACK; ++i) {
        int z0diff = pfch[i].hwZ0 - pvZ0;
        if (std::abs(z0diff) <= LINPUPPI_dzCut || pfch[i].hwId == PID_Muon) {
            outallch[i] = pfch[i];
        } else {
            clear(outallch[i]);
        }
    }
}

pt_t linpuppi_ref_sum2puppiPt(uint64_t sum, pt_t pt, bool isEM, int icand, bool debug) {
    const int sum_bitShift = 15;
    const int x2_bits = 6;    // decimal bits the discriminator values
    const int alpha_bits = 4; // decimal bits of the alpha values
    const int alphaSlope_bits = 4; // decimal bits of the alphaSlope values
    const int ptSlope_bits = 6;    // decimal bits of the ptSlope values 
    const int weight_bits = 8;

    const int ptSlopeNe = LINPUPPI_ptSlopeNe * (1 << ptSlope_bits);
    const int ptSlopePh = LINPUPPI_ptSlopePh * (1 << ptSlope_bits);
    const int ptZeroNe = LINPUPPI_ptZeroNe / LINPUPPI_ptLSB; // in pt scale
    const int ptZeroPh = LINPUPPI_ptZeroPh / LINPUPPI_ptLSB; // in pt scale
    const int alphaCrop = LINPUPPI_alphaCrop * (1 << x2_bits);
    const int alphaSlopeNe = LINPUPPI_alphaSlope * std::log(2.) * (1 << alphaSlope_bits); // we put a log(2) here since we compute alpha as log2(sum) instead of ln(sum)
    const int alphaSlopePh = LINPUPPI_alphaSlope * std::log(2.) * (1 << alphaSlope_bits);
    const int alphaZeroNe = LINPUPPI_alphaZero / std::log(2.) * (1 << alpha_bits);
    const int alphaZeroPh = LINPUPPI_alphaZero / std::log(2.) * (1 << alpha_bits);
    const int priorNe = LINPUPPI_priorNe * (1 << x2_bits);
    const int priorPh = LINPUPPI_priorPh * (1 << x2_bits);

    // -- simplest version
    //int alpha = sum > 0 ? int(std::log2(float(sum) * LINPUPPI_pt2DR2_scale / (1<<sum_bitShift)) * (1 << alpha_bits) + 0.5) :  0;
    // -- re-written bringing terms out of the log
    //int alpha = sum > 0 ? int(std::log2(float(sum))*(1 << alpha_bits) + (std::log2(LINPUPPI_pt2DR2_scale) - sum_bitShift)*(1 << alpha_bits) + 0.5 ) :  0;
    // -- re-written for a LUT implementation of the log2
    const int log2lut_bits = 10;
    int alpha = 0; uint64_t logarg = sum;
    if (logarg > 0) {
        alpha = int((std::log2(LINPUPPI_pt2DR2_scale) - sum_bitShift)*(1 << alpha_bits) + 0.5);
        while (logarg >= (1 << log2lut_bits)) { logarg = logarg >> 1; alpha += (1 << alpha_bits); }
        alpha += int(std::log2(float(logarg))*(1 << alpha_bits)); // the maximum value of this term is log2lut_bits * (1 << alpha_bits) ~ 10*16 = 160 => fits in ap_uint<4+alpha_bits>
    }
    int alphaZero  = (isEM ? alphaZeroPh : alphaZeroNe);
    int alphaSlope = (isEM ? alphaSlopePh : alphaSlopeNe);
    int x2a = std::min(std::max( alphaSlope * (alpha - alphaZero)  >> (alphaSlope_bits + alpha_bits - x2_bits), -alphaCrop), alphaCrop);

    // -- re-written to fit in a single LUT
    int x2a_lut = - alphaSlope * alphaZero; logarg = sum;
    if (logarg > 0) {
        x2a_lut += alphaSlope * int((std::log2(LINPUPPI_pt2DR2_scale) - sum_bitShift)*(1 << alpha_bits) + 0.5);
        while (logarg >= (1 << log2lut_bits)) { 
            logarg = logarg >> 1; x2a_lut += alphaSlope * (1 << alpha_bits); 
        }
        x2a_lut += alphaSlope * int(std::log2(float(logarg))*(1 << alpha_bits)); 
        /*if (in <= 3) printf("ref [%d]:  x2a(sum = %9lu): logarg = %9lu, sumterm = %9d, table[logarg] = %9d, ret pre-crop = %9d\n", 
          in, sum, logarg, 
          alphaSlope * int((std::log2(LINPUPPI_pt2DR2_scale) - sum_bitShift)*(1 << alpha_bits) + 0.5) - alphaSlope * alphaZero,
          alphaSlope * int(std::log2(float(logarg))*(1 << alpha_bits)), 
          x2a_lut); */
    } else {
        //if (in <= 3) printf("ref [%d]:  x2a(sum = %9lu): logarg = %9lu, ret pre-crop = %9d\n", 
        //        in, sum, logarg, x2a_lut); 
    }
    x2a_lut = std::min(std::max( x2a_lut >> (alphaSlope_bits + alpha_bits - x2_bits), -alphaCrop), alphaCrop);
    assert( x2a_lut == x2a );

    int ptZero  = (isEM ? ptZeroPh : ptZeroNe);
    int ptSlope = (isEM ? ptSlopePh : ptSlopeNe);
    int x2pt    = ptSlope * (pt - ptZero) >> (ptSlope_bits + 2 - x2_bits);

    int prior  = (isEM ? priorPh : priorNe);

    int x2 = x2a + x2pt - prior;

    // -- simplest version
    // int weight = 1.0/(1.0 + std::exp(- float(x2)/(1<<x2_bits))) * ( 1 << weight_bits ) + 0.5;
    // -- make explicit the max range of |x| values that yield a weight different from 0 or 1, to allow implementing as LUT
    // 1/(1+exp(-x*2^-A)) = 2^-B -->  |x| < 2^A*log(2^B) = 2^A * B * log(2)   with [ A = x2_bits, B = weight_bits ] 
    assert(weight_bits * std::log(2.0f) < 8);
    const int x2_abs_max = 1 << ( x2_bits + 3 );

    int weight = 0;
    if (x2 >= x2_abs_max) { 
        weight = (1 << weight_bits); // note: in practice here we should just skip the multiplication
    } else if (x2 >= -x2_abs_max) {
        weight = int(1.0/(1.0 + std::exp(- float(x2)/(1<<x2_bits))) * ( 1 << weight_bits ) + 0.5);
        assert(weight >= 0 && weight <= (1 << weight_bits));
    }

    pt_t ptPuppi = ( pt * weight ) >> weight_bits;

    if (debug) printf("ref candidate %02d pt %7.2f  em %1d: alpha %+7.2f   x2a %+5d = %+7.3f  x2pt %+5d = %+7.3f   x2 %+5d = %+7.3f  --> weight %4d = %.4f  puppi pt %7.2f\n",
               icand, pt*LINPUPPI_ptLSB, int(isEM), 
               std::max<float>(alpha/float(1<<alpha_bits)*std::log(2.),-99.99f), 
               x2a, x2a/float(1<<x2_bits), x2pt, x2pt/float(1<<x2_bits), x2, x2/float(1<<x2_bits), 
               weight, weight/float( 1 << weight_bits ), 
               ptPuppi*LINPUPPI_ptLSB);

    return ptPuppi;
}


void fwdlinpuppi_ref(const HadCaloObj caloin[NCALO], PFNeutralObj pfallne[NCALO], PFNeutralObj pfselne[NNEUTRALS], bool debug) {
    const int DR2MAX = LINPUPPI_DR2MAX; 
    const int DR2MIN = LINPUPPI_DR2MIN; 
    const int PTMAX2 = (LINPUPPI_ptMax*LINPUPPI_ptMax/LINPUPPI_ptLSB/LINPUPPI_ptLSB);

    for (int i = 0; i < NCALO; ++i) {
        pfallne[i].hwPt      = caloin[i].hwPt;
        pfallne[i].hwEta     = caloin[i].hwEta;
        pfallne[i].hwPhi     = caloin[i].hwPhi;
        pfallne[i].hwId      = caloin[i].hwIsEM ? PID_Photon : PID_Neutral;
        pfallne[i].hwPtPuppi = 0;
    }

    const int sum_bitShift = 15;
    const int ptCut = LINPUPPI_ptCut / LINPUPPI_ptLSB; 

    for (int in = 0; in < NCALO; ++in) {
        if (caloin[in].hwPt == 0) continue;
        uint64_t sum = 0; // 2 ^ sum_bitShift times (int pt^2)/(int dr2)
        for (int it = 0; it < NCALO; ++it) {
            if (it == in || caloin[it].hwPt == 0) continue;
            int dr2 = dr2_int(caloin[it].hwEta, caloin[it].hwPhi, caloin[in].hwEta, caloin[in].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
                ap_uint<9> dr2short = (dr2 >= DR2MIN ? dr2 : DR2MIN) >> 5; // reduce precision to make divide LUT cheaper
                uint64_t pt2 = caloin[it].hwPt*caloin[it].hwPt;
                uint64_t term = std::min<uint64_t>(pt2 >> 5, PTMAX2 >> 5)*((1 << sum_bitShift)/int(dr2short));
                //      dr2short >= (DR2MIN >> 5) = 2
                //      num <= (PTMAX2 >> 5) << sum_bitShift = (2^11) << 15 = 2^26
                //      ==> term <= 2^25
                //printf("ref term [%2d,%2d]: dr = %8d  pt2_shift = %8lu  term = %12lu\n", in, it, dr2, std::min<uint64_t>(pt2 >> 5, PTMAX2 >> 5), term);
                assert(uint64_t(PTMAX2 << (sum_bitShift-5))/(DR2MIN >> 5) <= (1 << 25));
                assert(term < (1 << 25));
                sum += term;
                //printf("    pT cand %5.1f    pT item %5.1f    dR = %.3f   term = %.1f [dbl] = %lu [int]\n",
                //            LINPUPPI_ptLSB*caloin[in].hwPt, LINPUPPI_ptLSB*caloin[it].hwPt, std::sqrt(dr2*LINPUPPI_DR2LSB),
                //            double(std::min<uint64_t>(pt2 >> 5, 131071)<<15)/double(std::max<int>(dr2,DR2MIN) >> 5),
                //            term);
            }
        }

        pfallne[in].hwPtPuppi = linpuppi_ref_sum2puppiPt(sum, caloin[in].hwPt, caloin[in].hwIsEM, in, debug);

    }

    puppisort_and_crop_ref<PFNeutralObj,NCALO,NNEUTRALS>(pfallne, ptCut, pfselne);
}

void linpuppi_ref(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], PFNeutralObj outallne[NALLNEUTRALS], PFNeutralObj outselne[NNEUTRALS], bool debug) {
    const int DR2MAX = LINPUPPI_DR2MAX; 
    const int DR2MIN = LINPUPPI_DR2MIN; 
    const int PTMAX2 = (LINPUPPI_ptMax*LINPUPPI_ptMax/LINPUPPI_ptLSB/LINPUPPI_ptLSB);

    const int sum_bitShift = 15;
    const int ptCut = LINPUPPI_ptCut / LINPUPPI_ptLSB; 

    for (int in = 0; in < NALLNEUTRALS; ++in) {
        outselne[in] = pfallne[in];
        if (pfallne[in].hwPt == 0) continue;
        uint64_t sum = 0; // 2 ^ sum_bitShift times (int pt^2)/(int dr2)
        for (int it = 0; it < NTRACK; ++it) {
            if (it == in || track[it].hwPt == 0) continue;
            if (std::abs<int>(track[it].hwZ0 - pvZ0) > LINPUPPI_dzCut) continue;
            int dr2 = dr2_int(pfallne[in].hwEta, pfallne[in].hwPhi, track[it].hwEta, track[it].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
                ap_uint<9> dr2short = (dr2 >= DR2MIN ? dr2 : DR2MIN) >> 5; // reduce precision to make divide LUT cheaper
                uint64_t pt2 = track[it].hwPt*track[it].hwPt;
                uint64_t term = std::min<uint64_t>(pt2 >> 5, PTMAX2 >> 5)*((1 << sum_bitShift)/int(dr2short));
                //      dr2short >= (DR2MIN >> 5) = 2
                //      num <= (PTMAX2 >> 5) << sum_bitShift = (2^11) << 15 = 2^26
                //      ==> term <= 2^25
                //printf("ref term [%2d,%2d]: dr = %8d  pt2_shift = %8lu  term = %12lu\n", in, it, dr2, std::min<uint64_t>(pt2 >> 5, PTMAX2 >> 5), term);
                assert(uint64_t(PTMAX2 << (sum_bitShift-5))/(DR2MIN >> 5) <= (1 << 25));
                assert(term < (1 << 25));
                sum += term;
                //printf("    pT cand %5.1f    pT item %5.1f    dR = %.3f   term = %.1f [dbl] = %lu [int]\n",
                //            LINPUPPI_ptLSB*pfallne[in].hwPt, LINPUPPI_ptLSB*track[it].hwPt, std::sqrt(dr2*LINPUPPI_DR2LSB),
                //            double(std::min<uint64_t>(pt2 >> 5, 131071)<<15)/double(std::max<int>(dr2,DR2MIN) >> 5),
                //            term);
            }
        }

        bool isEM = (pfallne[in].hwId == PID_Photon);
        outallne[in].hwPtPuppi = linpuppi_ref_sum2puppiPt(sum, pfallne[in].hwPt, isEM, in, debug);

     }

    puppisort_and_crop_ref<PFNeutralObj,NCALO,NNEUTRALS>(outallne, ptCut, outselne);

}

float linpuppi_flt_sum2puppiPt(float sum, float pt, bool isEM, int icand, bool debug) {
    float alphaZero  = LINPUPPI_alphaZero, alphaSlope = LINPUPPI_alphaSlope, alphaCrop = LINPUPPI_alphaCrop;
    float alpha = sum > 0 ? std::log(sum) : -9e9;
    float x2a = std::min(std::max( alphaSlope * (alpha - alphaZero), -alphaCrop), alphaCrop);

    float ptZero  = (isEM ? LINPUPPI_ptZeroPh  : LINPUPPI_ptZeroNe);
    float ptSlope = (isEM ? LINPUPPI_ptSlopePh : LINPUPPI_ptSlopeNe);
    float x2pt    = ptSlope * (pt - ptZero);

    float prior  = (isEM ? LINPUPPI_priorPh : LINPUPPI_priorNe);

    float x2 = x2a + x2pt - prior;

    float weight = 1.0/(1.0 + std::exp(-x2));

    float puppiPt = pt *weight;
    if (debug) printf("flt candidate %02d pt %7.2f  em %1d: alpha %+7.2f   x2a         %+7.3f  x2pt         %+7.3f   x2         %+7.3f  --> weight        %.4f  puppi pt %7.2f\n",
                   icand, pt, int(isEM), std::max(alpha,-99.99f), x2a, x2pt, x2, weight, puppiPt);

    return puppiPt;
}


void fwdlinpuppi_flt(const HadCaloObj caloin[NCALO], PFNeutralObj pfallne[NCALO], PFNeutralObj pfselne[NNEUTRALS], bool debug) {
    const int DR2MAX = LINPUPPI_DR2MAX; 
    const int DR2MIN = LINPUPPI_DR2MIN; 
    const float f_ptMax = LINPUPPI_ptMax;
    const float f_ptCut = LINPUPPI_ptCut;

    for (int i = 0; i < NCALO; ++i) {
        pfallne[i].hwPt      = caloin[i].hwPt;
        pfallne[i].hwEta     = caloin[i].hwEta;
        pfallne[i].hwPhi     = caloin[i].hwPhi;
        pfallne[i].hwId      = caloin[i].hwIsEM ? PID_Photon : PID_Neutral;
        pfallne[i].hwPtPuppi = 0;
    }

    for (int in = 0; in < NCALO; ++in) {
        if (caloin[in].hwPt == 0) continue;
        float sum = 0;
        for (int it = 0; it < NCALO; ++it) {
            if (it == in || caloin[it].hwPt == 0) continue;
            int dr2 = dr2_int(caloin[it].hwEta, caloin[it].hwPhi, caloin[in].hwEta, caloin[in].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
                sum += std::pow(std::min<float>(caloin[it].hwPt*LINPUPPI_ptLSB,f_ptMax),2) / (std::max<int>(dr2,DR2MIN) * LINPUPPI_DR2LSB);
            }
        }

        float ptPuppi = linpuppi_flt_sum2puppiPt(sum, caloin[in].hwPt*LINPUPPI_ptLSB, caloin[in].hwIsEM, in, debug);
        pfallne[in].hwPtPuppi = int(ptPuppi / LINPUPPI_ptLSB);
    }

    puppisort_and_crop_ref<PFNeutralObj,NCALO,NNEUTRALS>(pfallne, int(f_ptCut/LINPUPPI_ptLSB), pfselne);
}

void linpuppi_flt(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], PFNeutralObj outallne[NALLNEUTRALS], PFNeutralObj outselne[NNEUTRALS], bool debug) {
    const int DR2MAX = LINPUPPI_DR2MAX; 
    const int DR2MIN = LINPUPPI_DR2MIN; 
    const float f_ptMax = LINPUPPI_ptMax;
    const float f_ptCut = LINPUPPI_ptCut;

    for (int in = 0; in < NALLNEUTRALS; ++in) {
        outselne[in] = pfallne[in];
        if (pfallne[in].hwPt == 0) continue;
        float sum = 0;
        for (int it = 0; it < NTRACK; ++it) {
            if (it == in || track[it].hwPt == 0) continue;
            if (std::abs<int>(track[it].hwZ0 - pvZ0) > LINPUPPI_dzCut) continue;
            int dr2 = dr2_int(pfallne[in].hwEta, pfallne[in].hwPhi, track[it].hwEta, track[it].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
                sum += std::pow(std::min<float>(track[it].hwPt*LINPUPPI_ptLSB,f_ptMax),2) / (std::max<int>(dr2,DR2MIN) * LINPUPPI_DR2LSB);
            }
        }
        bool isEM = (pfallne[in].hwId == PID_Photon);
        float ptPuppi = linpuppi_flt_sum2puppiPt(sum, pfallne[in].hwPt*LINPUPPI_ptLSB, isEM, in, debug);
        outallne[in].hwPtPuppi = int( ptPuppi / LINPUPPI_ptLSB);
    }
    puppisort_and_crop_ref<PFNeutralObj,NCALO,NNEUTRALS>(outallne, int(f_ptCut/LINPUPPI_ptLSB), outselne);
}



