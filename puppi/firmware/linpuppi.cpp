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


inline int dr2_int(eta_t eta1, phi_t phi1, eta_t eta2, phi_t phi2) {
    ap_int<eta_t::width+1> deta = (eta1-eta2);
    ap_int<phi_t::width+1> dphi = (phi1-phi2);
    //ap_int<phi_t::width> dphi = (phi1-phi2); // intentional wrap-around
#ifdef LINPUPPI_DR2_LATENCY4
    int deta2 = deta*deta;
    int dphi2 = dphi*dphi;
    #pragma HLS resource variable=deta2 latency=4
    #pragma HLS resource variable=dphi2 latency=4
    int ret = deta2 + dphi2;
    return ret;
#else
    return deta*deta + dphi*dphi;
#endif
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

#define fwdlinpuppi_init_x2a_table_size 1024
#define x2a_t ap_int<16>

// GP: I can't template it since slope and zero are #defined as floats, so I use a macro
#define MAKE_X2A_LUT(FUNC_NAME, ALPHA_SLOPE, ALPHA_ZERO)                                                                    \
void FUNC_NAME(x2a_t table[fwdlinpuppi_init_x2a_table_size]) {                                                              \
    for (int i = 0; i < fwdlinpuppi_init_x2a_table_size; ++i) {                                                             \
        /* NOTE: HLS doesn't wants this constants to be inside of the loop in order to properly infer a ROM :-(  */         \
        const int sum_bitShift = LINPUPPI_sum_bitShift;                                                                     \
        const int alpha_bits = LINPUPPI_alpha_bits; /* decimal bits of the alpha values */                                  \
        const int alphaSlope_bits = LINPUPPI_alphaSlope_bits; /* decimal bits of the alphaSlope values */                   \
        const int alphaSlope = ALPHA_SLOPE * std::log(2) * (1 << alphaSlope_bits);                                          \
                         /* we put a log(2) here since we compute alpha as log2(sum) instead of ln(sum) */                  \
        const int alphaZero = ALPHA_ZERO / std::log(2) * (1 << alpha_bits);                                                 \
        const int C0 = - alphaSlope * alphaZero;                                                                            \
        const int C1 =   alphaSlope * int((std::log2(LINPUPPI_pt2DR2_scale) - sum_bitShift)*(1 << alpha_bits) + 0.5);       \
        int val = C0 + (i >  0 ? alphaSlope * int(std::log2(float(i))*(1 << alpha_bits)) + C1 : 0);                         \
        if (!(val >= -(1<<(x2a_t::width-1)) && val < (1<<(x2a_t::width-1)))) {                                              \
            printf("ERROR: overflow in x2a table[%d] with ap_int<%d> at index %d, val = %d, maxval = %d\n",                 \
                    fwdlinpuppi_init_x2a_table_size, x2a_t::width, i, val, 1<<(x2a_t::width-1));                            \
        }                                                                                                                   \
        assert(val >= -(1<<(x2a_t::width-1)) && val < (1<<(x2a_t::width-1)));                                               \
        table[i] = val;                                                                                                     \
    }                                                                                                                       \
}


MAKE_X2A_LUT(fwdlinpuppi_init_x2a_short, LINPUPPI_alphaSlope, LINPUPPI_alphaZero)  
#if defined(LINPUPPI_etaBins) && LINPUPPI_etaBins == 2
MAKE_X2A_LUT(fwdlinpuppi_init_x2a_short_1, LINPUPPI_alphaSlope_1, LINPUPPI_alphaZero_1)  
#endif


int fwdlinpuppi_calc_x2a(ap_uint<32> sum, int alphaSlope, int alphaCrop, const x2a_t & table) {
    const int log2lut_bits = 10;
    const int x2_bits = LINPUPPI_x2_bits;    // decimal bits the discriminator values
    const int alpha_bits = LINPUPPI_alpha_bits; // decimal bits of the alpha values
    const int alphaSlope_bits = LINPUPPI_alphaSlope_bits; // decimal bits of the alphaSlope values

    assert(sum >= 0);    
    int sumterm = 0, logarg = sum;
    for (int b = 31-log2lut_bits; b >=0; --b) {
        if (sum[b+log2lut_bits]) {
            logarg  = logarg >> (b + 1); 
            sumterm = (b + 1) * alphaSlope * (1 << alpha_bits); 
            break;
        }
    }

#ifndef __SYNTHESIS__
    if (logarg < 0 || logarg >= fwdlinpuppi_init_x2a_table_size) {
        printf("hw  x2a(sum = %9d): sumterm = %9d, logarg = %9d, ERROR\n", int(sum), sumterm, logarg);
    }
#endif
    assert(logarg >= 0 && logarg < fwdlinpuppi_init_x2a_table_size);
    int ret = (table[logarg] + sumterm) >> (alphaSlope_bits + alpha_bits - x2_bits);
#ifndef __SYNTHESIS__
    //printf("hw  x2a(sum = %9d): logarg = %9d, sumterm = %9d, table[logarg] = %9d, ret pre-crop = %9d\n", 
    //            int(sum), logarg, sumterm, int(table[logarg]), ret);
#endif
    if (ret < -alphaCrop) {
        return -alphaCrop;
    } else if (ret > alphaCrop) {
        return +alphaCrop;
    } else {
        return ret; 
    }
}


int fwdlinpuppi_calc_x2a_step2(ap_uint<32> sum, int alphaSlope, int alphaCrop, const x2a_t table[fwdlinpuppi_init_x2a_table_size]) {
    #pragma HLS inline

    const int log2lut_bits = 10;
    const int x2_bits = LINPUPPI_x2_bits;    // decimal bits the discriminator values
    const int alpha_bits = LINPUPPI_alpha_bits; // decimal bits of the alpha values
    const int alphaSlope_bits = LINPUPPI_alphaSlope_bits; // decimal bits of the alphaSlope values

    assert(sum >= 0);    
    int sumterm = 0, logarg = sum;
    for (int b = 31-log2lut_bits; b >=0; --b) {
        if (sum[b+log2lut_bits]) {
            logarg  = logarg >> (b + 1); 
            sumterm = (b + 1) * alphaSlope * (1 << alpha_bits); 
            break;
        }
    }

#ifndef __SYNTHESIS__
    if (logarg < 0 || logarg >= fwdlinpuppi_init_x2a_table_size) {
        printf("hw  x2a(sum = %9d): sumterm = %9d, logarg = %9d, ERROR\n", int(sum), sumterm, logarg);
    }
#endif
    assert(logarg >= 0 && logarg < fwdlinpuppi_init_x2a_table_size);
    int ret = (table[logarg] + sumterm) >> (alphaSlope_bits + alpha_bits - x2_bits);
#ifndef __SYNTHESIS__
    //printf("hw  x2a(sum = %9d): logarg = %9d, sumterm = %9d, table[logarg] = %9d, ret pre-crop = %9d\n", 
    //            int(sum), logarg, sumterm, int(table[logarg]), ret);
#endif
    if (ret < -alphaCrop) {
        return -alphaCrop;
    } else if (ret > alphaCrop) {
        return +alphaCrop;
    } else {
        return ret; 
    }
}


int fwdlinpuppi_calc_x2a(ap_uint<32> sum) {
    static x2a_t table[fwdlinpuppi_init_x2a_table_size];
#ifdef __SYNTHESIS__
    fwdlinpuppi_init_x2a_short(table); 
#else // initialize the table only once, otherwise this is really slow
    static bool is_init = false;
    if (!is_init) { fwdlinpuppi_init_x2a_short(table); is_init = true; }
#endif

    const int alphaSlope = LINPUPPI_alphaSlope * std::log(2) * (1 << LINPUPPI_alphaSlope_bits); // we put a log(2) here since we compute alpha as log2(sum) instead of ln(sum)
    const int alphaCrop = LINPUPPI_alphaCrop * (1 << LINPUPPI_x2_bits);

    return fwdlinpuppi_calc_x2a_step2(sum, alphaSlope, alphaCrop, table);
}

#if defined(LINPUPPI_etaBins) && LINPUPPI_etaBins == 2
int fwdlinpuppi_calc_x2a(ap_uint<32> sum, bool ietaBin) {
    static x2a_t table0[fwdlinpuppi_init_x2a_table_size], table1[fwdlinpuppi_init_x2a_table_size];
#ifdef __SYNTHESIS__
    fwdlinpuppi_init_x2a_short(table0); fwdlinpuppi_init_x2a_short_1(table1); 
#else // initialize the table only once, otherwise this is really slow
    static bool is_init = false;
    if (!is_init) { fwdlinpuppi_init_x2a_short(table0); fwdlinpuppi_init_x2a_short_1(table1); is_init = true; }
#endif

    const int alphaSlope_0 = LINPUPPI_alphaSlope * std::log(2) * (1 << LINPUPPI_alphaSlope_bits); // we put a log(2) here since we compute alpha as log2(sum) instead of ln(sum)
    const int alphaCrop_0 = LINPUPPI_alphaCrop * (1 << LINPUPPI_x2_bits);
    const int alphaSlope_1 = LINPUPPI_alphaSlope_1 * std::log(2) * (1 << LINPUPPI_alphaSlope_bits); // we put a log(2) here since we compute alpha as log2(sum) instead of ln(sum)
    const int alphaCrop_1 = LINPUPPI_alphaCrop_1 * (1 << LINPUPPI_x2_bits);
    int alphaSlope = ietaBin ? alphaSlope_1 : alphaSlope_0;
    int alphaCrop = ietaBin ? alphaCrop_1 : alphaCrop_0;

    return fwdlinpuppi_calc_x2a_step2(sum, alphaSlope, alphaCrop, (ietaBin ? table1 : table0));
}
#endif

#define fwdlinpuppi_x2w_table_size 1024
void fwdlinpuppi_init_w(ap_uint<9> table[fwdlinpuppi_x2w_table_size]) {
    const int xavg = fwdlinpuppi_x2w_table_size/2;
    for (int i = 0; i <= 1023; ++i) {
        int x2 = i - xavg;
        int val = 1.0/(1.0 + std::exp(- float(x2)/(1<<6))) * ( 1 << 8 ) + 0.5;
        table[i] = val;
    }
}

pt_t fwdlinpuppi_calc_wpt(pt_t pt, int x2) {
    static ap_uint<9> table[fwdlinpuppi_x2w_table_size];
#ifdef __SYNTHESIS__
    fwdlinpuppi_init_w(table);
#else // initialize the table only once, otherwise this is really slow
    static bool is_init = false;
    if (!is_init) { fwdlinpuppi_init_w(table); is_init = true; }
#endif

    const int xavg = fwdlinpuppi_x2w_table_size/2;
    int index;
    if (x2 < -xavg) index = 0;
    else if (x2 >= xavg) index = fwdlinpuppi_x2w_table_size-1;
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

    const int x2_bits = LINPUPPI_x2_bits;    // decimal bits the discriminator values
    const int ptSlope_bits = LINPUPPI_ptSlope_bits;    // decimal bits of the ptSlope values 
    const int weight_bits = LINPUPPI_weight_bits;

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
#ifdef HLS_pipeline_II
 #if HLS_pipeline_II == 1
    #pragma HLS pipeline II=1
 #elif HLS_pipeline_II == 2
    #pragma HLS pipeline II=2
 #elif HLS_pipeline_II == 3
    #pragma HLS pipeline II=3
 #elif HLS_pipeline_II == 4
    #pragma HLS pipeline II=4
 #elif HLS_pipeline_II == 6
    #pragma HLS pipeline II=6
 #endif
#else
    #pragma HLS pipeline II=2
#endif

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
#ifdef HLS_pipeline_II
 #if HLS_pipeline_II == 1
    #pragma HLS pipeline II=1
 #elif HLS_pipeline_II == 2
    #pragma HLS pipeline II=2
 #elif HLS_pipeline_II == 3
    #pragma HLS pipeline II=3
 #elif HLS_pipeline_II == 4
    #pragma HLS pipeline II=4
 #elif HLS_pipeline_II == 6
    #pragma HLS pipeline II=6
 #endif
#else
    #pragma HLS pipeline II=2
#endif

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
    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=outallch complete
    #pragma HLS LATENCY min=1
#ifdef HLS_pipeline_II
 #if HLS_pipeline_II == 1
    #pragma HLS pipeline II=1
 #elif HLS_pipeline_II == 2
    #pragma HLS pipeline II=2
 #elif HLS_pipeline_II == 3
    #pragma HLS pipeline II=3
 #elif HLS_pipeline_II == 4
    #pragma HLS pipeline II=4
 #elif HLS_pipeline_II == 6
    #pragma HLS pipeline II=6
 #endif
#else
    #pragma HLS pipeline II=2
#endif

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

    const int x2_bits = LINPUPPI_x2_bits;    // decimal bits the discriminator values
    const int ptSlope_bits = LINPUPPI_ptSlope_bits;    // decimal bits of the ptSlope values 
    const int weight_bits = LINPUPPI_weight_bits;

#ifndef LINPUPPI_etaBins
    const int ptSlopeNe = LINPUPPI_ptSlopeNe * (1 << ptSlope_bits);
    const int ptSlopePh = LINPUPPI_ptSlopePh * (1 << ptSlope_bits);
    const int ptZeroNe = LINPUPPI_ptZeroNe / LINPUPPI_ptLSB; // in pt scale
    const int ptZeroPh = LINPUPPI_ptZeroPh / LINPUPPI_ptLSB; // in pt scale
    const int priorNe = LINPUPPI_priorNe * (1 << x2_bits);
    const int priorPh = LINPUPPI_priorPh * (1 << x2_bits);
    const int ptCut = LINPUPPI_ptCut; 
#elif LINPUPPI_etaBins == 2
    const int ptSlopeNe_0 = LINPUPPI_ptSlopeNe * (1 << ptSlope_bits);
    const int ptSlopePh_0 = LINPUPPI_ptSlopePh * (1 << ptSlope_bits);
    const int ptZeroNe_0 = LINPUPPI_ptZeroNe / LINPUPPI_ptLSB; // in pt scale
    const int ptZeroPh_0 = LINPUPPI_ptZeroPh / LINPUPPI_ptLSB; // in pt scale
    const int priorNe_0 = LINPUPPI_priorNe * (1 << x2_bits);
    const int priorPh_0 = LINPUPPI_priorPh * (1 << x2_bits);
    const int ptCut_0 = LINPUPPI_ptCut; 
    const int ptSlopeNe_1 = LINPUPPI_ptSlopeNe_1 * (1 << ptSlope_bits);
    const int ptSlopePh_1 = LINPUPPI_ptSlopePh_1 * (1 << ptSlope_bits);
    const int ptZeroNe_1 = LINPUPPI_ptZeroNe_1 / LINPUPPI_ptLSB; // in pt scale
    const int ptZeroPh_1 = LINPUPPI_ptZeroPh_1 / LINPUPPI_ptLSB; // in pt scale
    const int priorNe_1 = LINPUPPI_priorNe_1 * (1 << x2_bits);
    const int priorPh_1 = LINPUPPI_priorPh_1 * (1 << x2_bits);
    const int ptCut_1 = LINPUPPI_ptCut_1; 
#endif

    ap_int<12>  x2a[NALLNEUTRALS], x2ptp[NALLNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=x2a complete    
    #pragma HLS ARRAY_PARTITION variable=x2ptp complete    

#if defined(LINPUPPI_etaBins) && LINPUPPI_etaBins == 2
    bool ietaBin[NALLNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=ietaBin complete    
    for (int in = 0; in < NALLNEUTRALS; ++in) {
        ietaBin[in] = (caloin[in].hwEta <= LINPUPPI_etaCut) ? LINPUPPI_invertEta : (1-LINPUPPI_invertEta);
    }
#endif

    for (int in = 0; in < NALLNEUTRALS; ++in) {
#ifndef LINPUPPI_etaBins
        x2a[in] = fwdlinpuppi_calc_x2a(sums[in]);
#elif LINPUPPI_etaBins == 2
        x2a[in] = fwdlinpuppi_calc_x2a(sums[in], ietaBin[in]);
#endif
    }

    for (int in = 0; in < NALLNEUTRALS; ++in) {
        if (caloin[in].hwId == PID_Photon) {
#ifndef LINPUPPI_etaBins
            int val = (ptSlopePh*caloin[in].hwPt - ptSlopePh*ptZeroPh) >> (ptSlope_bits + 2 - x2_bits);
            x2ptp[in] =  val < 2047 ? val - priorPh : 2047; // saturate
#elif LINPUPPI_etaBins == 2
            int val = ((ietaBin[in] ? ptSlopePh_1 : ptSlopePh_0)*caloin[in].hwPt - (ietaBin[in] ? ptSlopePh_1*ptZeroPh_1 : ptSlopePh_0*ptZeroPh_0)) >> (ptSlope_bits + 2 - x2_bits);
            x2ptp[in] =  val < 2047 ? val - (ietaBin[in] ? priorPh_1 : priorPh_0) : 2047; // saturate
#endif
        } else {
#ifndef LINPUPPI_etaBins
            int val = (ptSlopeNe*caloin[in].hwPt - ptSlopeNe*ptZeroNe) >> (ptSlope_bits + 2 - x2_bits);
            x2ptp[in] =  val < 2047 ? val - priorNe : 2047; // saturate
#elif LINPUPPI_etaBins == 2
            int val = ((ietaBin[in] ? ptSlopeNe_1 : ptSlopeNe_0)*caloin[in].hwPt - (ietaBin[in] ? ptSlopeNe_1*ptZeroNe_1 : ptSlopeNe_0*ptZeroNe_0)) >> (ptSlope_bits + 2 - x2_bits);
            x2ptp[in] =  val < 2047 ? val - (ietaBin[in] ? priorNe_1 : priorNe_0) : 2047; // saturate
#endif
        }
    }

    for (int in = 0; in < NALLNEUTRALS; ++in) {
        int x2 = x2a[in]+x2ptp[in];
        pt_t puppiPt = fwdlinpuppi_calc_wpt(caloin[in].hwPt, x2);
#ifndef LINPUPPI_etaBins
        if (puppiPt >= LINPUPPI_ptCut) {
#elif LINPUPPI_etaBins == 2
        if (puppiPt >= (ietaBin[in] ? LINPUPPI_ptCut_1 : LINPUPPI_ptCut)) {
#endif
            out[in] = caloin[in];
            out[in].hwPtPuppi = puppiPt;
        } else {
            clear(out[in]);
        }
#ifndef __SYNTHESIS__
        if (caloin[in].hwPt == 0) continue;
#ifndef LINPUPPI_etaBins
        if (gdebug_) printf("hw  candidate %02d pt %7.2f  em %1d: alpha %+7.2f   x2a %+5d = %+7.3f  x2pt %+5d = %+7.3f   x2 %+5d = %+7.3f  -->                       puppi pt %7.2f\n",
                   in, caloin[in].hwPt*LINPUPPI_ptLSB, int(caloin[in].hwId == PID_Photon), 
                   sums[in] > 0 ? std::log2(float(sums[in]) * LINPUPPI_pt2DR2_scale / (1<<15))*std::log(2.) : 0., 
                   int(x2a[in]), x2a[in]/float(1<<x2_bits), 
                   (int(x2ptp[in]) + (caloin[in].hwId == PID_Photon ? priorPh : priorNe) ), 
                   (int(x2ptp[in]) + (caloin[in].hwId == PID_Photon ? priorPh : priorNe) )/float(1<<x2_bits), 
                   x2, x2/float(1<<x2_bits), 
                   puppiPt*LINPUPPI_ptLSB);
#elif LINPUPPI_etaBins == 2
        if (gdebug_) printf("hw  candidate %02d pt %7.2f  em %1d  ieta %1d: alpha %+7.2f   x2a %+5d = %+7.3f  x2pt %+5d = %+7.3f   x2 %+5d = %+7.3f  -->                       puppi pt %7.2f\n",
                   in, caloin[in].hwPt*LINPUPPI_ptLSB, int(caloin[in].hwId == PID_Photon), int(ietaBin[in]),
                   sums[in] > 0 ? std::log2(float(sums[in]) * LINPUPPI_pt2DR2_scale / (1<<15))*std::log(2.) : 0., 
                   int(x2a[in]), x2a[in]/float(1<<x2_bits), 
                   (int(x2ptp[in]) + (caloin[in].hwId == PID_Photon ? (ietaBin[in] ? priorPh_1 : priorPh_0) : (ietaBin[in] ? priorNe_1 : priorNe_0)) ), 
                   (int(x2ptp[in]) + (caloin[in].hwId == PID_Photon ? (ietaBin[in] ? priorPh_1 : priorPh_0) : (ietaBin[in] ? priorNe_1 : priorNe_0)) )/float(1<<x2_bits), 
                   x2, x2/float(1<<x2_bits), 
                   puppiPt*LINPUPPI_ptLSB);
#endif // etaBins
#endif // synthesis
    }
}


void linpuppiNoCrop(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], PFNeutralObj outallne[NALLNEUTRALS]) {
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS ARRAY_PARTITION variable=outallne complete
#ifdef HLS_pipeline_II
 #if HLS_pipeline_II == 1
    #pragma HLS pipeline II=1
 #elif HLS_pipeline_II == 2
    #pragma HLS pipeline II=2
 #elif HLS_pipeline_II == 3
    #pragma HLS pipeline II=3
 #elif HLS_pipeline_II == 4
    #pragma HLS pipeline II=4
 #elif HLS_pipeline_II == 6
    #pragma HLS pipeline II=6
 #endif
#else
    #pragma HLS pipeline II=2
#endif

    ap_uint<32> sums[NALLNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=sums complete
    linpuppiSum(track, pvZ0, pfallne, sums);

    linpuppiSum2All(pfallne, sums, outallne);
}

void linpuppi(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], PFNeutralObj outselne[NNEUTRALS]) {
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS ARRAY_PARTITION variable=outselne complete
#ifdef HLS_pipeline_II
 #if HLS_pipeline_II == 1
    #pragma HLS pipeline II=1
 #elif HLS_pipeline_II == 2
    #pragma HLS pipeline II=2
 #elif HLS_pipeline_II == 3
    #pragma HLS pipeline II=3
 #elif HLS_pipeline_II == 4
    #pragma HLS pipeline II=4
 #elif HLS_pipeline_II == 6
    #pragma HLS pipeline II=6
 #endif
#else
    #pragma HLS pipeline II=2
#endif

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

#if defined(PACKING_DATA_SIZE) && defined(PACKING_NCHANN)
void packed_fwdlinpuppi(const ap_uint<PACKING_DATA_SIZE> input[PACKING_NCHANN], ap_uint<PACKING_DATA_SIZE> output[PACKING_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=output complete
#ifdef HLS_pipeline_II
 #if HLS_pipeline_II == 1
    #pragma HLS pipeline II=1
 #elif HLS_pipeline_II == 2
    #pragma HLS pipeline II=2
 #elif HLS_pipeline_II == 3
    #pragma HLS pipeline II=3
 #elif HLS_pipeline_II == 4
    #pragma HLS pipeline II=4
 #elif HLS_pipeline_II == 6
    #pragma HLS pipeline II=6
 #endif
#else
    #pragma HLS pipeline II=2
#endif

    HadCaloObj caloin[NCALO]; PFNeutralObj pfselne[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=pfselne complete
    l1pf_pattern_unpack<NCALO,0>(input, caloin);
    fwdlinpuppi(caloin, pfselne);
    l1pf_pattern_pack<NNEUTRALS,0>(pfselne, output);
}

void packed_fwdlinpuppiNoCrop(const ap_uint<PACKING_DATA_SIZE> input[PACKING_NCHANN], ap_uint<PACKING_DATA_SIZE> output[PACKING_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=output complete
#ifdef HLS_pipeline_II
 #if HLS_pipeline_II == 1
    #pragma HLS pipeline II=1
 #elif HLS_pipeline_II == 2
    #pragma HLS pipeline II=2
 #elif HLS_pipeline_II == 3
    #pragma HLS pipeline II=3
 #elif HLS_pipeline_II == 4
    #pragma HLS pipeline II=4
 #elif HLS_pipeline_II == 6
    #pragma HLS pipeline II=6
 #endif
#else
    #pragma HLS pipeline II=2
#endif

    HadCaloObj caloin[NCALO]; PFNeutralObj pfallne[NCALO];
    #pragma HLS ARRAY_PARTITION variable=caloin complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    l1pf_pattern_unpack<NCALO,0>(input, caloin);
    fwdlinpuppiNoCrop(caloin, pfallne);
    l1pf_pattern_pack<NCALO,0>(pfallne, output);
}

void packed_linpuppi_chs(const ap_uint<PACKING_DATA_SIZE> input[PACKING_NCHANN], ap_uint<PACKING_DATA_SIZE> output[PACKING_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=output complete
#ifdef HLS_pipeline_II
 #if HLS_pipeline_II == 1
    #pragma HLS pipeline II=1
 #elif HLS_pipeline_II == 2
    #pragma HLS pipeline II=2
 #elif HLS_pipeline_II == 3
    #pragma HLS pipeline II=3
 #elif HLS_pipeline_II == 4
    #pragma HLS pipeline II=4
 #elif HLS_pipeline_II == 6
    #pragma HLS pipeline II=6
 #endif
#else
    #pragma HLS pipeline II=2
#endif

    z0_t pvZ0; PFChargedObj pfch[NTRACK], outallch[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=outallch complete
    linpuppi_chs_unpack_in(input, pvZ0, pfch);
    linpuppi_chs(pvZ0, pfch, outallch);
    l1pf_pattern_pack<NTRACK,0>(outallch, output);
}

void packed_linpuppi(const ap_uint<PACKING_DATA_SIZE> input[PACKING_NCHANN], ap_uint<PACKING_DATA_SIZE> output[PACKING_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=output complete
#ifdef HLS_pipeline_II
 #if HLS_pipeline_II == 1
    #pragma HLS pipeline II=1
 #elif HLS_pipeline_II == 2
    #pragma HLS pipeline II=2
 #elif HLS_pipeline_II == 3
    #pragma HLS pipeline II=3
 #elif HLS_pipeline_II == 4
    #pragma HLS pipeline II=4
 #elif HLS_pipeline_II == 6
    #pragma HLS pipeline II=6
 #endif
#else
    #pragma HLS pipeline II=2
#endif

    TkObj track[NTRACK]; z0_t pvZ0; PFNeutralObj pfallne[NALLNEUTRALS], outselne[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS ARRAY_PARTITION variable=outselne complete
    linpuppi_unpack_in(input, track, pvZ0, pfallne);
    linpuppi(track, pvZ0, pfallne, outselne);
    l1pf_pattern_pack<NNEUTRALS,0>(outselne, output);
}

void packed_linpuppiNoCrop(const ap_uint<PACKING_DATA_SIZE> input[PACKING_NCHANN], ap_uint<PACKING_DATA_SIZE> output[PACKING_NCHANN]) {
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=output complete
#ifdef HLS_pipeline_II
 #if HLS_pipeline_II == 1
    #pragma HLS pipeline II=1
 #elif HLS_pipeline_II == 2
    #pragma HLS pipeline II=2
 #elif HLS_pipeline_II == 3
    #pragma HLS pipeline II=3
 #elif HLS_pipeline_II == 4
    #pragma HLS pipeline II=4
 #elif HLS_pipeline_II == 6
    #pragma HLS pipeline II=6
 #endif
#else
    #pragma HLS pipeline II=2
#endif

    TkObj track[NTRACK]; z0_t pvZ0; PFNeutralObj pfallne[NALLNEUTRALS], outallne[NALLNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS ARRAY_PARTITION variable=outallne complete
    linpuppi_unpack_in(input, track, pvZ0, pfallne);
    linpuppiNoCrop(track, pvZ0, pfallne, outallne);
    l1pf_pattern_pack<NALLNEUTRALS,0>(outallne, output);
}

void linpuppi_pack_in(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], ap_uint<PACKING_DATA_SIZE> input[PACKING_NCHANN]) {
    assert(NTRACK+1+NALLNEUTRALS <= PACKING_NCHANN);
    const int TK_OFFS = 0, PV_OFFS = TK_OFFS + NTRACK, PFNE_OFFS = PV_OFFS + 1;
    l1pf_pattern_pack<NTRACK, TK_OFFS>(track, input);
    linpuppi_pack_pv(pvZ0, input[PV_OFFS]);
    l1pf_pattern_pack<NALLNEUTRALS, PFNE_OFFS>(pfallne, input);
}

void linpuppi_unpack_in(const ap_uint<PACKING_DATA_SIZE> input[PACKING_NCHANN], TkObj track[NTRACK], z0_t & pvZ0, PFNeutralObj pfallne[NALLNEUTRALS]) {
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS inline
    assert(NTRACK+1+NALLNEUTRALS <= PACKING_NCHANN);
    const int TK_OFFS = 0, PV_OFFS = TK_OFFS + NTRACK, PFNE_OFFS = PV_OFFS + 1;
    l1pf_pattern_unpack<NTRACK, TK_OFFS>(input, track);
    linpuppi_unpack_pv(input[PV_OFFS], pvZ0);
    l1pf_pattern_unpack<NALLNEUTRALS, PFNE_OFFS>(input, pfallne);
}

void linpuppi_chs_pack_in(z0_t pvZ0, const PFChargedObj pfch[NTRACK], ap_uint<PACKING_DATA_SIZE> input[PACKING_NCHANN]) {
    assert(NTRACK+1 <= PACKING_NCHANN);
    linpuppi_pack_pv(pvZ0, input[0]);
    l1pf_pattern_pack<NTRACK, 1>(pfch, input);
}

void linpuppi_chs_unpack_in(const ap_uint<PACKING_DATA_SIZE> input[PACKING_NCHANN], z0_t & pvZ0, PFChargedObj pfch[NTRACK]) {
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS inline
    assert(NTRACK+1 <= PACKING_NCHANN);
    linpuppi_unpack_pv(input[0], pvZ0);
    l1pf_pattern_unpack<NTRACK,1>(input, pfch);
}     

void linpuppi_pack_pv(z0_t pvZ0, ap_uint<PACKING_DATA_SIZE> & word) {
    word = 0;
    word(z0_t::width-1,0) = pvZ0;
}
void linpuppi_unpack_pv(ap_uint<PACKING_DATA_SIZE> word, z0_t & pvZ0) {
    pvZ0 = word(z0_t::width-1,0);
}

#endif
