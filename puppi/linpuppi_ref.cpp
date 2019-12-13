#include "linpuppi_ref.h"
#include <cmath>
#include <algorithm>


template<typename T>
void puppisort_and_crop_ref(unsigned int nIn, unsigned int nOut, const T in[/*nIn*/], T out[/*nOut*/]) {
    std::vector<T> tmp(nOut);

    for (int iout = 0; iout < nOut; ++iout) {
        clear(tmp[iout]);
    }

    for (int it = 0; it < nIn; ++it) {
        for (int iout = nOut-1; iout >= 0; --iout) {
            if (tmp[iout].hwPtPuppi <= in[it].hwPtPuppi) {
                if (iout == 0 || tmp[iout-1].hwPtPuppi > in[it].hwPtPuppi) {
                    tmp[iout] = in[it];
                } else {
                    tmp[iout] = tmp[iout-1];
                }
            }
        }
    }

    for (int iout = 0; iout < nOut; ++iout) {
        out[iout] = tmp[iout];
    }

}


void linpuppi_chs_ref(const linpuppi_config &cfg, z0_t pvZ0, const PFChargedObj pfch[/*cfg.nTrack*/], PFChargedObj outallch[/*cfg.nTrack*/]) {
    for (unsigned int i = 0; i < cfg.nTrack; ++i) {
        int z0diff = pfch[i].hwZ0 - pvZ0;
        if (std::abs(z0diff) <= cfg.dzCut || pfch[i].hwId == PID_Muon) {
            outallch[i] = pfch[i];
        } else {
            clear(outallch[i]);
        }
    }
}

unsigned int linpuppi_ieta_ref(const linpuppi_config &cfg, etaphi_t eta) {
    assert(cfg.absEtaBins.empty());
    if (cfg.absEtaBins.empty()) return 0;
}

pt_t linpuppi_ref_sum2puppiPt(const linpuppi_config &cfg, uint64_t sum, pt_t pt, unsigned int ieta, bool isEM, int icand, bool debug) {
    const int sum_bitShift = 15;
    const int x2_bits = 6;    // decimal bits the discriminator values
    const int alpha_bits = 4; // decimal bits of the alpha values
    const int alphaSlope_bits = 4; // decimal bits of the alphaSlope values
    const int ptSlope_bits = 6;    // decimal bits of the ptSlope values 
    const int weight_bits = 8;

    const int ptSlopeNe = cfg.ptSlopeNe[ieta] * (1 << ptSlope_bits);
    const int ptSlopePh = cfg.ptSlopePh[ieta] * (1 << ptSlope_bits);
    const int ptZeroNe = cfg.ptZeroNe[ieta] / LINPUPPI_ptLSB; // in pt scale
    const int ptZeroPh = cfg.ptZeroPh[ieta] / LINPUPPI_ptLSB; // in pt scale
    const int alphaCrop = cfg.alphaCrop[ieta] * (1 << x2_bits);
    const int alphaSlopeNe = cfg.alphaSlope[ieta] * std::log(2.) * (1 << alphaSlope_bits); // we put a log(2) here since we compute alpha as log2(sum) instead of ln(sum)
    const int alphaSlopePh = cfg.alphaSlope[ieta] * std::log(2.) * (1 << alphaSlope_bits);
    const int alphaZeroNe = cfg.alphaZero[ieta] / std::log(2.) * (1 << alpha_bits);
    const int alphaZeroPh = cfg.alphaZero[ieta] / std::log(2.) * (1 << alpha_bits);
    const int priorNe = cfg.priorNe[ieta] * (1 << x2_bits);
    const int priorPh = cfg.priorPh[ieta] * (1 << x2_bits);

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


void fwdlinpuppi_ref(const linpuppi_config &cfg, const HadCaloObj caloin[/*cfg.nIn*/], PFNeutralObj outallne_nocut[/*cfg.nIn*/], PFNeutralObj outallne[/*cfg.nIn*/], PFNeutralObj outselne[/*cfg.nOut*/], bool debug) {
    const int DR2MAX = LINPUPPI_DR2MAX; 
    const int DR2MIN = LINPUPPI_DR2MIN; 
    const int PTMAX2 = (LINPUPPI_ptMax*LINPUPPI_ptMax);

    for (int i = 0; i < cfg.nIn; ++i) {
        outallne_nocut[i].hwPt      = caloin[i].hwPt;
        outallne_nocut[i].hwEta     = caloin[i].hwEta;
        outallne_nocut[i].hwPhi     = caloin[i].hwPhi;
        outallne_nocut[i].hwId      = caloin[i].hwIsEM ? PID_Photon : PID_Neutral;
        outallne_nocut[i].hwPtPuppi = 0;
        clear(outallne[i]);
    }

    const int sum_bitShift = 15;

    for (int in = 0; in < cfg.nIn; ++in) {
        if (caloin[in].hwPt == 0) continue;
        uint64_t sum = 0; // 2 ^ sum_bitShift times (int pt^2)/(int dr2)
        for (int it = 0; it < cfg.nIn; ++it) {
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
        unsigned int ieta = linpuppi_ieta_ref(cfg, caloin[in].hwEta);
        outallne_nocut[in].hwPtPuppi = linpuppi_ref_sum2puppiPt(cfg, sum, caloin[in].hwPt, ieta, caloin[in].hwIsEM, in, debug);
        if (outallne_nocut[in].hwPtPuppi >= cfg.ptCut[ieta]) {
            outallne[in] = outallne_nocut[in];
        }
    }
    puppisort_and_crop_ref(cfg.nIn, cfg.nOut, outallne, outselne);
}

void linpuppi_ref(const linpuppi_config &cfg, const TkObj track[/*cfg.nTrack*/], z0_t pvZ0, const PFNeutralObj pfallne[/*cfg.nIn*/], PFNeutralObj outallne_nocut[/*cfg.nIn*/], PFNeutralObj outallne[/*cfg.nIn*/], PFNeutralObj outselne[/*cfg.nOut*/], bool debug) {
    const int DR2MAX = LINPUPPI_DR2MAX; 
    const int DR2MIN = LINPUPPI_DR2MIN; 
    const int PTMAX2 = (LINPUPPI_ptMax*LINPUPPI_ptMax);

    const int sum_bitShift = 15;

    for (int in = 0; in < cfg.nIn; ++in) {
        outallne_nocut[in] = pfallne[in]; clear(outallne[in]);
        if (pfallne[in].hwPt == 0) continue;
        uint64_t sum = 0; // 2 ^ sum_bitShift times (int pt^2)/(int dr2)
        for (int it = 0; it < cfg.nTrack; ++it) {
            if (it == in || track[it].hwPt == 0) continue;
            if (std::abs<int>(track[it].hwZ0 - pvZ0) > cfg.dzCut) continue;
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

        unsigned int ieta = linpuppi_ieta_ref(cfg, pfallne[in].hwEta);
        bool isEM = (pfallne[in].hwId == PID_Photon);
        outallne_nocut[in].hwPtPuppi = linpuppi_ref_sum2puppiPt(cfg, sum, pfallne[in].hwPt, ieta, isEM, in, debug);
        if (outallne_nocut[in].hwPtPuppi >= cfg.ptCut[ieta]) {
            outallne[in] = outallne_nocut[in];
        }
     }
    puppisort_and_crop_ref(cfg.nIn, cfg.nOut, outallne, outselne);

}

float linpuppi_flt_sum2puppiPt(const linpuppi_config &cfg, float sum, float pt, unsigned int ieta, bool isEM, int icand, bool debug) {
    float alphaZero  = cfg.alphaZero[ieta], alphaSlope = cfg.alphaSlope[ieta], alphaCrop = cfg.alphaCrop[ieta];
    float alpha = sum > 0 ? std::log(sum) : -9e9;
    float x2a = std::min(std::max( alphaSlope * (alpha - alphaZero), -alphaCrop), alphaCrop);

    float ptZero  = (isEM ? cfg.ptZeroPh[ieta]  : cfg.ptZeroNe[ieta]);
    float ptSlope = (isEM ? cfg.ptSlopePh[ieta] : cfg.ptSlopeNe[ieta]);
    float x2pt    = ptSlope * (pt - ptZero);

    float prior  = (isEM ? cfg.priorPh[ieta] : cfg.priorNe[ieta]);

    float x2 = x2a + x2pt - prior;

    float weight = 1.0/(1.0 + std::exp(-x2));

    float puppiPt = pt *weight;
    if (debug) printf("flt candidate %02d pt %7.2f  em %1d: alpha %+7.2f   x2a         %+7.3f  x2pt         %+7.3f   x2         %+7.3f  --> weight        %.4f  puppi pt %7.2f\n",
                   icand, pt, int(isEM), std::max(alpha,-99.99f), x2a, x2pt, x2, weight, puppiPt);

    return puppiPt;
}


void fwdlinpuppi_flt(const linpuppi_config &cfg, const HadCaloObj caloin[/*cfg.nIn*/], PFNeutralObj outallne_nocut[/*cfg.nIn*/], PFNeutralObj outallne[/*cfg.nIn*/], PFNeutralObj outselne[/*cfg.nOut*/], bool debug) {
    const int DR2MAX = cfg.dR2Max; 
    const int DR2MIN = cfg.dR2Min; 
    const float f_ptMax = LINPUPPI_ptMax * LINPUPPI_ptLSB;

    for (int i = 0; i < cfg.nIn; ++i) {
        outallne_nocut[i].hwPt      = caloin[i].hwPt;
        outallne_nocut[i].hwEta     = caloin[i].hwEta;
        outallne_nocut[i].hwPhi     = caloin[i].hwPhi;
        outallne_nocut[i].hwId      = caloin[i].hwIsEM ? PID_Photon : PID_Neutral;
        outallne_nocut[i].hwPtPuppi = 0;
        clear(outallne[i]);
    }

    for (int in = 0; in < cfg.nIn; ++in) {
        if (caloin[in].hwPt == 0) continue;
        float sum = 0;
        for (int it = 0; it < cfg.nIn; ++it) {
            if (it == in || caloin[it].hwPt == 0) continue;
            int dr2 = dr2_int(caloin[it].hwEta, caloin[it].hwPhi, caloin[in].hwEta, caloin[in].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
                sum += std::pow(std::min<float>(caloin[it].hwPt*LINPUPPI_ptLSB,f_ptMax),2) / (std::max<int>(dr2,DR2MIN) * LINPUPPI_DR2LSB);
            }
        }

        unsigned int ieta = linpuppi_ieta_ref(cfg, caloin[in].hwEta);
        float ptPuppi = linpuppi_flt_sum2puppiPt(cfg, sum, caloin[in].hwPt*LINPUPPI_ptLSB, ieta, caloin[in].hwIsEM, in, debug);
        outallne_nocut[in].hwPtPuppi = int(ptPuppi / LINPUPPI_ptLSB);
        if (outallne_nocut[in].hwPtPuppi >= cfg.ptCut[ieta]) {
            outallne[in] = outallne_nocut[in];
        }
    }

    puppisort_and_crop_ref(cfg.nIn, cfg.nOut, outallne, outselne);
}

void linpuppi_flt(const linpuppi_config &cfg, const TkObj track[/*cfg.nTrack*/], z0_t pvZ0, const PFNeutralObj pfallne[/*cfg.nIn*/], PFNeutralObj outallne_nocut[/*cfg.nIn*/], PFNeutralObj outallne[/*cfg.nIn*/], PFNeutralObj outselne[/*cfg.nOut*/], bool debug) {
    const int DR2MAX = cfg.dR2Max; 
    const int DR2MIN = cfg.dR2Min; 
    const float f_ptMax = LINPUPPI_ptMax * LINPUPPI_ptLSB;

    for (int in = 0; in < cfg.nIn; ++in) {
        outallne_nocut[in] = pfallne[in]; clear(outallne[in]);
        if (pfallne[in].hwPt == 0) continue;
        float sum = 0;
        for (int it = 0; it < cfg.nTrack; ++it) {
            if (it == in || track[it].hwPt == 0) continue;
            if (std::abs<int>(track[it].hwZ0 - pvZ0) > cfg.dzCut) continue;
            int dr2 = dr2_int(pfallne[in].hwEta, pfallne[in].hwPhi, track[it].hwEta, track[it].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
                sum += std::pow(std::min<float>(track[it].hwPt*LINPUPPI_ptLSB,f_ptMax),2) / (std::max<int>(dr2,DR2MIN) * LINPUPPI_DR2LSB);
            }
        }
        unsigned int ieta = linpuppi_ieta_ref(cfg, pfallne[in].hwEta);
        bool isEM = (pfallne[in].hwId == PID_Photon);
        float ptPuppi = linpuppi_flt_sum2puppiPt(cfg, sum, pfallne[in].hwPt*LINPUPPI_ptLSB, ieta, isEM, in, debug);
        outallne_nocut[in].hwPtPuppi = int(ptPuppi / LINPUPPI_ptLSB);
        if (outallne_nocut[in].hwPtPuppi >= cfg.ptCut[ieta]) {
            outallne[in] = outallne_nocut[in];
        }
    }
    puppisort_and_crop_ref(cfg.nIn, cfg.nOut, outallne, outselne);
}



