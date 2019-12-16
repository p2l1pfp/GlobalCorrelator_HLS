#include "linpuppi.h"
#include <algorithm>
#include <cassert>

#ifndef __SYNTHESIS__
#include <cstdio>
int gdebug_;
void linpuppi_set_debug(bool debug) { gdebug_ = debug; }
#else
void linpuppi_set_debug(bool debug) {}
#endif

void fwdlinpuppiSum(const HadCaloObj caloin[NCALO], ap_uint<32> sums[NCALO]);
void fwdlinpuppiSum2Pt(const HadCaloObj caloin[NCALO], const ap_uint<32> sums[NCALO], pt_t puppiPts[NCALO]);
void fwdlinpuppiPt(const HadCaloObj caloin[NCALO], pt_t puppiPts[NCALO]);


int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    ap_int<etaphi_t::width+1> deta = (eta1-eta2);
    ap_int<etaphi_t::width+1> dphi = (phi1-phi2);
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
    static ap_uint<16> _table[512];
    _lut_shift15_invert_init(_table);
    return (num * _table[den]);
}

void fwdlinpuppi_init_x2a_short(ap_int<14> table[1024]) {
    for (int i = 0; i < 1024; ++i) {
        // NOTE: HLS doesn't wants this constants to be inside of the loop in order to properly infer a ROM :-/
        const int sum_bitShift = 15;
        const int alpha_bits = 4; // decimal bits of the alpha values
        const int alphaSlope_bits = 4; // decimal bits of the alphaSlope values
        const int alphaSlope = LINPUPPI_alphaSlope * std::log(2) * (1 << alphaSlope_bits); // we put a log(2) here since we compute alpha as log2(sum) instead of ln(sum)
        const int alphaZero = LINPUPPI_alphaZero / std::log(2) * (1 << alpha_bits);
        const int C0 = - alphaSlope * alphaZero;
        const int C1 =   alphaSlope * int((std::log2(LINPUPPI_pt2DR2_scale) - sum_bitShift)*(1 << alpha_bits) + 0.5); 
        table[i] = C0 + (i >  0 ? alphaSlope * int(std::log2(float(i))*(1 << alpha_bits)) + C1 : 0);
    }
}

int fwdlinpuppi_calc_x2a(ap_uint<32> sum) {
    static ap_int<14> table[1024];
#ifdef __SYNTHESIS__
    fwdlinpuppi_init_x2a_short(table);
#else // initialize the table only once, otherwise this is really slow
    static bool is_init = false;
    if (!is_init) { fwdlinpuppi_init_x2a_short(table); is_init = true; }
#endif

    const int log2lut_bits = 10;
    const int x2_bits = 6;    // decimal bits the discriminator values
    const int alpha_bits = 4; // decimal bits of the alpha values
    const int alphaSlope_bits = 4; // decimal bits of the alphaSlope values
    const int alphaSlope = LINPUPPI_alphaSlope * std::log(2) * (1 << alphaSlope_bits); // we put a log(2) here since we compute alpha as log2(sum) instead of ln(sum)
    const int alphaCrop = LINPUPPI_alphaCrop * (1 << x2_bits);

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

void fwdlinpuppi_init_w(ap_uint<9> table[1024]) {
    const int xavg = 512;
    for (int i = 0; i <= 1023; ++i) {
        int x2 = i - xavg;
        int val = 1.0/(1.0 + std::exp(- float(x2)/(1<<6))) * ( 1 << 8 ) + 0.5;
        table[i] = val;
    }
}

pt_t fwdlinpuppi_calc_wpt(pt_t pt, int x2) {
    static ap_uint<9> table[1024];
#ifdef __SYNTHESIS__
    fwdlinpuppi_init_w(table);
#else // initialize the table only once, otherwise this is really slow
    static bool is_init = false;
    if (!is_init) { fwdlinpuppi_init_w(table); is_init = true; }
#endif

    const int xavg = 512;
    int index;
    if (x2 < -xavg) index = 0;
    else if (x2 > xavg) index = 1023;
    else index = x2 + xavg;
    return pt_t( int(pt * table[index]) >> 8 );
}

void fwdlinpuppiSum(const HadCaloObj caloin[NCALO], ap_uint<32> sums[NCALO]) {
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=sums complete
    #pragma HLS inline

    const int DR2MAX = LINPUPPI_DR2MAX; 
    const int DR2MIN = LINPUPPI_DR2MIN; 
    const int DR2MIN_SHIFT =  DR2MIN >> 5; 
    const int PTMAX2_SHIFT = (LINPUPPI_ptMax)*(LINPUPPI_ptMax) >> 5;

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

void fwdlinpuppiPt(const HadCaloObj caloin[NCALO], pt_t puppiPts[NCALO]) {
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=puppiPts complete
    #pragma HLS pipeline II=2

    ap_uint<32> sums[NCALO];
    #pragma HLS ARRAY_PARTITION variable=sums complete
    fwdlinpuppiSum(caloin, sums);

    fwdlinpuppiSum2Pt(caloin, sums, puppiPts);
}

void fwdlinpuppiSum2Pt(const HadCaloObj caloin[NCALO], const ap_uint<32> sums[NCALO], pt_t puppiPts[NCALO]) {
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=sums complete
    #pragma HLS ARRAY_PARTITION variable=puppiPts complete
    #pragma HLS inline

    const int x2_bits = 6;    // decimal bits the discriminator values
    const int ptSlope_bits = 6;    // decimal bits of the ptSlope values 
    const int weight_bits = 8;

    const int ptSlopeNe = LINPUPPI_ptSlopeNe * (1 << ptSlope_bits);
    const int ptSlopePh = LINPUPPI_ptSlopePh * (1 << ptSlope_bits);
    const int ptZeroNe = LINPUPPI_ptZeroNe / LINPUPPI_ptLSB; // in pt scale
    const int ptZeroPh = LINPUPPI_ptZeroPh / LINPUPPI_ptLSB; // in pt scale
    const int priorNe = LINPUPPI_priorNe * (1 << x2_bits);
    const int priorPh = LINPUPPI_priorPh * (1 << x2_bits);
    const int ptCut = LINPUPPI_ptCut; 

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
        if (gdebug_) printf("hw  candidate %02d pt %7.2f  em %1d: alpha %+7.2f   x2a %+5d = %+7.3f  x2pt %+5d = %+7.3f   x2 %+5d = %+7.3f  -->                       puppi pt %7.2f\n",
                   in, caloin[in].hwPt*LINPUPPI_ptLSB, int(caloin[in].hwIsEM), 
                   sums[in] > 0 ? std::log2(float(sums[in]) * LINPUPPI_pt2DR2_scale / (1<<15))*std::log(2.) : 0., 
                   int(x2a[in]), x2a[in]/float(1<<x2_bits), 
                   (int(x2ptp[in]) + (caloin[in].hwIsEM ? priorPh : priorNe) ), 
                   (int(x2ptp[in]) + (caloin[in].hwIsEM ? priorPh : priorNe) )/float(1<<x2_bits), 
                   x2, x2/float(1<<x2_bits), 
                   puppiPts[in]*LINPUPPI_ptLSB);
#endif
    }
}

void fwdlinpuppiNoCrop(const HadCaloObj caloin[NCALO], PFNeutralObj pfallne[NCALO]) {
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=pfselne complete
    #pragma HLS pipeline II=2

    pt_t puppiPts[NCALO];
    #pragma HLS ARRAY_PARTITION variable=puppiPts complete    

    fwdlinpuppiPt(caloin, puppiPts);

    const int ptCut = LINPUPPI_ptCut;
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

void fwdlinpuppi(const HadCaloObj caloin[NCALO], PFNeutralObj pfselne[NNEUTRALS]) {
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=pfselne complete
    #pragma HLS pipeline II=2

    pt_t puppiPts[NCALO];
    #pragma HLS ARRAY_PARTITION variable=puppiPts complete    

    fwdlinpuppiPt(caloin, puppiPts);

    PFNeutralObj work[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=work complete    

    for (int out = 0; out < NNEUTRALS; ++out) {
        work[out].hwPt = 0;
        work[out].hwEta = 0;
        work[out].hwPhi = 0;
        work[out].hwId = 0;
        work[out].hwPtPuppi = 0;
    }

    const int ptCut = LINPUPPI_ptCut;
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

template<typename T>
inline bool linpuppi_fromPV(const T & obj, z0_t pvZ0) {
       int z0diff = obj.hwZ0 - pvZ0;
       if (z0diff < 0) z0diff = -z0diff; 
       return (z0diff <= LINPUPPI_dzCut);
}

void linpuppi_chs(z0_t pvZ0, const PFChargedObj pfch[NTRACK], PFChargedObj outallch[NTRACK]) {
    for (unsigned int i = 0; i < NTRACK; ++i) {
        if (linpuppi_fromPV(pfch[i], pvZ0) || pfch[i].hwId == PID_Muon) {
            outallch[i] = pfch[i];
        } else {
            clear(outallch[i]);
        }
    }
}

void linpuppiSum(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj caloin[NALLNEUTRALS], ap_uint<32> sums[NALLNEUTRALS]) {
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=sums complete
    #pragma HLS inline

    const int DR2MAX = LINPUPPI_DR2MAX; 
    const int DR2MIN = LINPUPPI_DR2MIN; 
    const int DR2MIN_SHIFT =  DR2MIN >> 5; 
    const int PTMAX2_SHIFT = (LINPUPPI_ptMax)*(LINPUPPI_ptMax) >> 5;

    bool fromPV[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=fromPV complete
    for (unsigned int i = 0; i < NTRACK; ++i) {
        fromPV[i] = linpuppi_fromPV(track[i], pvZ0);
    }

    ap_uint<17> pt2_shift[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=pt2_shift complete
    for (int it = 0; it < NTRACK; ++it) {
        int mypt2 = (track[it].hwPt*track[it].hwPt) >> 5; // reduce precision to make multiplication smaller later 
        pt2_shift[it] = (mypt2 < PTMAX2_SHIFT? mypt2 : PTMAX2_SHIFT);
    }

    for (int in = 0; in < NALLNEUTRALS; ++in) {
        ap_uint<32> sum = 0;
        for (int it = 0; it < NTRACK; ++it) {
            if (it == in) continue;
            int dr2 = dr2_int(track[it].hwEta, track[it].hwPhi, caloin[in].hwEta, caloin[in].hwPhi); 
            if (dr2 <= DR2MAX && fromPV[it]) { // if dr is inside puppi cone
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

void linpuppiSum2All(const PFNeutralObj caloin[NALLNEUTRALS], const ap_uint<32> sums[NALLNEUTRALS], PFNeutralObj out[NALLNEUTRALS]) {
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=sums complete
    #pragma HLS ARRAY_PARTITION variable=out complete

    const int x2_bits = 6;    // decimal bits the discriminator values
    const int ptSlope_bits = 6;    // decimal bits of the ptSlope values 
    const int weight_bits = 8;

    const int ptSlopeNe = LINPUPPI_ptSlopeNe * (1 << ptSlope_bits);
    const int ptSlopePh = LINPUPPI_ptSlopePh * (1 << ptSlope_bits);
    const int ptZeroNe = LINPUPPI_ptZeroNe / LINPUPPI_ptLSB; // in pt scale
    const int ptZeroPh = LINPUPPI_ptZeroPh / LINPUPPI_ptLSB; // in pt scale
    const int priorNe = LINPUPPI_priorNe * (1 << x2_bits);
    const int priorPh = LINPUPPI_priorPh * (1 << x2_bits);
    const int ptCut = LINPUPPI_ptCut; 

    ap_int<12>  x2a[NALLNEUTRALS], x2ptp[NALLNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=x2a complete    
    #pragma HLS ARRAY_PARTITION variable=x2ptp complete    

    for (int in = 0; in < NALLNEUTRALS; ++in) {
        x2a[in] = fwdlinpuppi_calc_x2a(sums[in]);
    }

    for (int in = 0; in < NALLNEUTRALS; ++in) {
        if (caloin[in].hwId == PID_Photon) {
            int val = (ptSlopePh*caloin[in].hwPt - ptSlopePh*ptZeroPh) >> (ptSlope_bits + 2 - x2_bits);
            x2ptp[in] =  val < 2047 ? val - priorPh : 2047; // saturate
        } else {
            int val = (ptSlopeNe*caloin[in].hwPt - ptSlopeNe*ptZeroNe) >> (ptSlope_bits + 2 - x2_bits);
            x2ptp[in] =  val < 2047 ? val - priorNe : 2047; // saturate
        }
    }

    for (int in = 0; in < NALLNEUTRALS; ++in) {
        int x2 = x2a[in]+x2ptp[in];
        pt_t puppiPt = fwdlinpuppi_calc_wpt(caloin[in].hwPt, x2);
        if (puppiPt >= LINPUPPI_ptCut) {
            out[in] = caloin[in];
            out[in].hwPtPuppi = puppiPt;
        } else {
            clear(out[in]);
        }
#ifndef __SYNTHESIS__
        if (caloin[in].hwPt == 0) continue;
        if (gdebug_) printf("hw  candidate %02d pt %7.2f  em %1d: alpha %+7.2f   x2a %+5d = %+7.3f  x2pt %+5d = %+7.3f   x2 %+5d = %+7.3f  -->                       puppi pt %7.2f\n",
                   in, caloin[in].hwPt*LINPUPPI_ptLSB, int(caloin[in].hwId == PID_Photon), 
                   sums[in] > 0 ? std::log2(float(sums[in]) * LINPUPPI_pt2DR2_scale / (1<<15))*std::log(2.) : 0., 
                   int(x2a[in]), x2a[in]/float(1<<x2_bits), 
                   (int(x2ptp[in]) + (caloin[in].hwId == PID_Photon ? priorPh : priorNe) ), 
                   (int(x2ptp[in]) + (caloin[in].hwId == PID_Photon ? priorPh : priorNe) )/float(1<<x2_bits), 
                   x2, x2/float(1<<x2_bits), 
                   puppiPt*LINPUPPI_ptLSB);
#endif
    }
}


void linpuppiNoCrop(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], PFNeutralObj outallne[NALLNEUTRALS]) {
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS ARRAY_PARTITION variable=outallne complete
    #pragma HLS pipeline II=2

    ap_uint<32> sums[NALLNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=sums complete
    linpuppiSum(track, pvZ0, pfallne, sums);

    linpuppiSum2All(pfallne, sums, outallne);
}

void linpuppi(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], PFNeutralObj outselne[NNEUTRALS]) {
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS ARRAY_PARTITION variable=outselne complete
    #pragma HLS pipeline II=2

    PFNeutralObj allne[NALLNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=allne complete

    linpuppiNoCrop(track, pvZ0, pfallne, allne);

    PFNeutralObj work[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=work complete

    for (int iout = 0; iout < NNEUTRALS; ++iout) {
        clear(work[iout]);
    }

    for (int in = 0; in < NALLNEUTRALS; ++in) {
        for (int iout = NNEUTRALS-1; iout >= 0; --iout) {
            if (work[iout].hwPtPuppi <= allne[in].hwPtPuppi) {
                if (iout == 0 || work[iout-1].hwPtPuppi > allne[in].hwPtPuppi) {
                    work[iout] = allne[in];
                } else {
                    work[iout] = work[iout-1];
                }
            }
        }
    }

    for (int iout = 0; iout < NNEUTRALS; ++iout) {
        outselne[iout] = work[iout];
    }
}

