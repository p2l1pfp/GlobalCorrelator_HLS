#include "linpuppi.h"
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}


void _lut_shift15_invert_init(ap_uint<16> _table[512]) { // returns 2^15 / x
    _table[0] = 32768; // this is 2^15
    for (int i = 1; i <= 511; ++i) {
        _table[i] = (32768 / i);
    }
}
int _lut_shift15_divide(ap_uint<17> num, ap_uint<9> den) { // returns (num * 2^15) / den
    assert(int(den) >= 0 && int(den) <= 511);
    ap_uint<16> _table[512];
    _lut_shift15_invert_init(_table);
    return (num * _table[den]);
}

void fwdlinpuppi_init_x2a_short(ap_int<14> table[1024]) {
    for (int i = 0; i < 1024; ++i) {
        // NOTE: HLS doesn't wants this constants to be inside of the loop in order to properly infer a ROM :-/
        const int sum_bitShift = 15;
        const int alpha_bits = 4; // decimal bits of the alpha values
        const int alphaSlope_bits = 4; // decimal bits of the alphaSlope values
        const int alphaSlope = 2.2 * std::log(2) * (1 << alphaSlope_bits); // we put a log(2) here since we compute alpha as log2(sum) instead of ln(sum)
        const int alphaZero = 9.0 / std::log(2) * (1 << alpha_bits);
        const int C0 = - alphaSlope * alphaZero;
        const int C1 =   alphaSlope * int((std::log(0.25*0.25 / 1.9e-5)/std::log(2.) - sum_bitShift)*(1 << alpha_bits) + 0.5); 
                                    // note: ^^^^^^^^^^^^^^^^^^^^^^^^^^  HLS does not know std::log2(x)
        table[i] = C0 + (i >  0 ? alphaSlope * int(std::log(float(i))/std::log(2.)*(1 << alpha_bits)) + C1 : 0);
                               // note: ^^^^^^^^^^^^^^^^^^^^^^^^^^  HLS does not know std::log2(x) --> std::log/std::log(2) 
    }
}

int fwdlinpuppi_calc_x2a(ap_uint<32> sum) {
    static ap_int<14> table[1024];
    fwdlinpuppi_init_x2a_short(table);

    const int log2lut_bits = 10;
    const int x2_bits = 6;    // decimal bits the discriminator values
    const int alpha_bits = 4; // decimal bits of the alpha values
    const int alphaSlope_bits = 4; // decimal bits of the alphaSlope values
    const int alphaSlope = 2.2 * std::log(2) * (1 << alphaSlope_bits); // we put a log(2) here since we compute alpha as log2(sum) instead of ln(sum)
    const int alphaCrop = 4.0 * (1 << x2_bits);

    assert(sum >= 0);    
    int sumterm = 0, logarg = sum;
    for (int b = 31-log2lut_bits; b >=0; --b) {
        if (sum[b+log2lut_bits]) {
            logarg  = logarg >> (b + 1); 
            sumterm = (b + 1) * alphaSlope * (1 << alpha_bits); 
            break;
        }
    }

    assert(logarg >= 0 && logarg <= 1023);
    int ret = (table[logarg] + sumterm) >> (alphaSlope_bits + alpha_bits - x2_bits);
    //printf("hw  x2a(sum = %9d): logarg = %9d, sumterm = %9d, table[logarg] = %9d, ret pre-crop = %9d\n", 
    //            sum, logarg, sumterm, int(table[logarg]), ret);
    if (ret < -alphaCrop) {
        return -alphaCrop;
    } else if (ret > alphaCrop) {
        return +alphaCrop;
    } else {
        return ret; 
    }
}

void fwdlinpuppi_init_w(ap_uint<8> table[1024]) {
    const int xavg = 512;
    for (int i = 0; i <= 1023; ++i) {
        int x2 = i - xavg;
        int val = 1.0/(1.0 + std::exp(- float(x2)/(1<<6))) * ( 1 << 8 ) + 0.5;
        table[i] = val <= 255 ? val : 255; // HLS doesn't know std::min<int>(a,b)  :-(
    }
}

pt_t fwdlinpuppi_calc_wpt(pt_t pt, int x2) {
    static ap_uint<8> table[1024];
    fwdlinpuppi_init_w(table);

    const int xavg = 512;
    const int x2_max = 400;  
    // the smallest number at which the table overflows because the number would round to 256
    // for (int x2 = 0; x2 < 1024; ++x2) { if (int(1/(1+std::exp(-float(x2)/(1<<6))) * 256 + 0.5) == 256) { printf( "x2_max = %d\n", x2); break; } }
    if (x2 >= x2_max) {
        return pt;
    } else if (x2 <= -xavg) {
        return pt_t(0);
    } else {
        return pt_t( int(pt * table[x2 + xavg]) >> 8 );
    }
}

void fwdlinpuppiSum_hw(const HadCaloObj caloin[NCALO], ap_uint<32> sums[NCALO]) {
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=sums complete
    #pragma HLS inline

    const int DR2MAX = 4727; // 0.3 cone
    const int DR2MIN =  84; // 0.04 cone
    const int DR2MIN_SHIFT =  DR2MIN >> 5; // 0.04 cone
    const int PTMAX = (50*4);
    const int PTMAX2_SHIFT = (PTMAX)*(PTMAX) >> 5;

    ap_uint<17> pt2_shift[NCALO];
    #pragma HLS ARRAY_PARTITION variable=pt2_shift complete
    for (int it = 0; it < NCALO; ++it) {
        int mypt2 = (caloin[it].hwPt*caloin[it].hwPt) >> 5; // reduce precision to make multiplication smaller later 
        pt2_shift[it] = (mypt2 < PTMAX2_SHIFT? mypt2 : PTMAX2_SHIFT);
    }

    for (int in = 0; in < NCALO; ++in) {
        ap_uint<32> sum = 0;
        for (int it = 0; it < NCALO; ++it) {
            if (it == in) continue;
            int dr2 = dr2_int(caloin[it].hwEta, caloin[it].hwPhi, caloin[in].hwEta, caloin[in].hwPhi); 
            if (dr2 <= DR2MAX) { // if dr is inside puppi cone
                ap_uint<9> dr2short = dr2 >> 5; // reduce precision to make divide LUT cheaper
                if (dr2short < DR2MIN_SHIFT) dr2short = DR2MIN_SHIFT;
                int term = _lut_shift15_divide(pt2_shift[it], dr2short);
                //printf("hw  term [%2d,%2d]: dr = %8d  pt2_shift = %8d  term = %12d\n", in, it, dr2, int(pt2_shift[it]), term);
                sum += term;
            }
        }
        sums[in] = sum;
    }
}

void fwdlinpuppiPt_hw(const HadCaloObj caloin[NCALO], pt_t puppiPts[NCALO]) {
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=puppiPts complete
    #pragma HLS pipeline II=2

    ap_uint<32> sums[NCALO];
    #pragma HLS ARRAY_PARTITION variable=sums complete
    fwdlinpuppiSum_hw(caloin, sums);

    fwdlinpuppiSum2Pt_hw(caloin, sums, puppiPts);
}

void fwdlinpuppiSum2Pt_hw(const HadCaloObj caloin[NCALO], const ap_uint<32> sums[NCALO], pt_t puppiPts[NCALO]) {
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=sums complete
    #pragma HLS ARRAY_PARTITION variable=puppiPts complete
    #pragma HLS inline

    const int x2_bits = 6;    // decimal bits the discriminator values
    const int ptSlope_bits = 6;    // decimal bits of the ptSlope values 
    const int weight_bits = 8;

    const int ptSlopeNe = 0.3 * (1 << ptSlope_bits);
    const int ptSlopePh = 0.4 * (1 << ptSlope_bits);
    const int ptZeroNe = 9.0 * 4; // in pt scale
    const int ptZeroPh = 5.0 * 4; // in pt scale
    const int priorNe = 7.0 * (1 << x2_bits);
    const int priorPh = 5.0 * (1 << x2_bits);
    const int ptCut = 4.0 * 4; // 4 GeV

    ap_int<12>  x2a[NCALO], x2ptp[NCALO];
    #pragma HLS ARRAY_PARTITION variable=x2a complete    
    #pragma HLS ARRAY_PARTITION variable=x2ptp complete    

    for (int in = 0; in < NCALO; ++in) {
        x2a[in] = fwdlinpuppi_calc_x2a(sums[in]);
    }

    for (int in = 0; in < NCALO; ++in) {
        if (caloin[in].hwIsEM) {
            int val = (ptSlopePh*caloin[in].hwPt - ptSlopePh*ptZeroPh) >> (ptSlope_bits + 2 - x2_bits);
            x2ptp[in] =  val < 2047 ? val - priorPh : 2047; // saturate
        } else {
            int val = (ptSlopeNe*caloin[in].hwPt - ptSlopeNe*ptZeroNe) >> (ptSlope_bits + 2 - x2_bits);
            x2ptp[in] =  val < 2047 ? val - priorNe : 2047; // saturate
        }
    }

    for (int in = 0; in < NCALO; ++in) {
        int x2 = x2a[in]+x2ptp[in];
        puppiPts[in] = fwdlinpuppi_calc_wpt(caloin[in].hwPt, x2);
#ifndef __SYNTHESIS__
        if (caloin[in].hwPt == 0) continue;
        printf("hw  candidate %02d pt %7.2f  em %1d: alpha %+7.2f   x2a %+5d = %+7.3f  x2pt %+5d = %+7.3f   x2 %+5d = %+7.3f  -->                       puppi pt %7.2f\n",
                   in, caloin[in].hwPt* 0.25, int(caloin[in].hwIsEM), 
                   sums[in] > 0 ? std::log2(float(sums[in]) * 0.25*0.25 / 1.9e-5 / (1<<15))*std::log(2.) : 0., 
                   int(x2a[in]), x2a[in]/float(1<<x2_bits), 
                   (int(x2ptp[in]) + (caloin[in].hwIsEM ? priorPh : priorNe) ), 
                   (int(x2ptp[in]) + (caloin[in].hwIsEM ? priorPh : priorNe) )/float(1<<x2_bits), 
                   x2, x2/float(1<<x2_bits), 
                   puppiPts[in]*0.25);
#endif
    }
}

void fwdlinpuppiNoCrop_hw(const HadCaloObj caloin[NCALO], PFNeutralObj pfallne[NCALO]) {
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=pfselne complete
    #pragma HLS pipeline II=2

    pt_t puppiPts[NCALO];
    #pragma HLS ARRAY_PARTITION variable=puppiPts complete    

    fwdlinpuppiPt_hw(caloin, puppiPts);

    const int ptCut = 4.0 * 4; // 4 GeV
    for (int in = 0; in < NCALO; ++in) {
        if (puppiPts[in] >= ptCut) {
            pfallne[in].hwPt      = caloin[in].hwPt;
            pfallne[in].hwEta     = caloin[in].hwEta;
            pfallne[in].hwPhi     = caloin[in].hwPhi;
            pfallne[in].hwId      = caloin[in].hwIsEM ? PID_Photon : PID_Neutral;
            pfallne[in].hwPtPuppi = puppiPts[in];
        } else {
            pfallne[in].hwPt      = 0;
            pfallne[in].hwEta     = 0;
            pfallne[in].hwPhi     = 0;
            pfallne[in].hwId      = 0;
            pfallne[in].hwPtPuppi = 0;
        }
    }
}

void fwdlinpuppi_hw(const HadCaloObj caloin[NCALO], PFNeutralObj pfselne[NNEUTRALS]) {
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=pfselne complete
    #pragma HLS pipeline II=2

    pt_t puppiPts[NCALO];
    #pragma HLS ARRAY_PARTITION variable=puppiPts complete    

    fwdlinpuppiPt_hw(caloin, puppiPts);

    PFNeutralObj work[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=work complete    

    for (int out = 0; out < NNEUTRALS; ++out) {
        work[out].hwPt = 0;
        work[out].hwEta = 0;
        work[out].hwPhi = 0;
        work[out].hwId = 0;
        work[out].hwPtPuppi = 0;
    }

    const int ptCut = 4.0 * 4; // 4 GeV
    for (int in = 0; in < NCALO; ++in) {
        if (puppiPts[in] < ptCut) continue;
        for (int iout = NNEUTRALS-1; iout >= 0; --iout) {
            if (work[iout].hwPtPuppi <= puppiPts[in]) {
                if (iout == 0 || work[iout-1].hwPtPuppi > puppiPts[in]) {
                    work[iout].hwPt      = caloin[in].hwPt;
                    work[iout].hwEta     = caloin[in].hwEta;
                    work[iout].hwPhi     = caloin[in].hwPhi;
                    work[iout].hwId      = caloin[in].hwIsEM ? PID_Photon : PID_Neutral;
                    work[iout].hwPtPuppi = puppiPts[in];
                } else {
                    work[iout] = work[iout-1];
                }
            }
        }
    }

    for (int iout = 0; iout < NNEUTRALS; ++iout) {
        pfselne[iout] = work[iout];
    }
}

