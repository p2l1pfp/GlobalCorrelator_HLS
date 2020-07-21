#include "tk_pt_inversion.h"

template<class pt_T>
void init_pt_inv_table(pt_T table_out[(1<<PT_INV_TAB_SIZE)]) {
    // index is a uint from 0 to 111...11=2^PT_INV_TAB_SIZE-1 that encodes 0.000 to 0.4999999
    // resulting value pt_T is a uint from 0 to 2^16-1
    table_out[0] = (1<<PT_INV_MAX_BITS);
    for (unsigned int i = 1; i < (1<<PT_INV_TAB_SIZE); i++) {
        float invpt = float(i)/(1<<PT_INV_TAB_SIZE) * 0.5; // in 1/GeV
        table_out[i] = PF_PT_SCALE / invpt;
    }
    return;
}

template<class pt_inv_T, class pt_T>
void convert_pt(pt_inv_T inv, pt_T &pt){
    //pt_inv_T is ap_fixed<> (signed)
    //pt_T is ap_uint<> !

    // Initialize the lookup tables
#ifdef __HLS_SYN__
    bool initialized = false;
    pt_t inv_table[(1<<PT_INV_TAB_SIZE)];
#else 
    static bool initialized = false;
    static pt_t inv_table[(1<<PT_INV_TAB_SIZE)];
#endif
    if (!initialized) {
        init_pt_inv_table<pt_T>(inv_table);
        initialized = true;
    }

    // if(inv<0) inv = -inv; assume input is positive (upstream)
    urinv_t uinv = inv;

    // cutoffs at high and low pt
    if(uinv >= urinv_t(0.5)){
        pt = 2. * PF_PT_SCALE;
        return;
    } else if (uinv <= urinv_t(1./(1<<PT_INV_MAX_BITS))){
        pt=(1<<PT_INV_MAX_BITS) * PF_PT_SCALE;
        return;
    }

    ap_uint<PT_INV_TAB_SIZE> index;
    const int offset = 1; // ignore the first bit since 0b0.01111.. = 0.499.. is largest value
    #pragma unroll
    for(int i=0; i<PT_INV_TAB_SIZE; i++){
        index[PT_INV_TAB_SIZE-1-i] = uinv[urinv_t::width-1-i-offset]; //msb down to lowest
    }

    pt = inv_table[index];
}
