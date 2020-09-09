//
#include "tk_input_converter.h"

template<class eta_T>
void init_eta_table(eta_T table_out[(1<<ETA_TAB_SIZE)]) {
    // index is a uint from 0 to 111...11=2^ETA_TAB_SIZE-1 that encodes 0.000 to 8
    // resulting value eta_T is a uint

    for (unsigned int i = 0; i < (1<<ETA_TAB_SIZE); i++) {
        // eta =  -ln(tan((pi/2 - arctan(TANLAM))/2))
        float tanlam = float(i)/(1<<ETA_TAB_SIZE) * 8.;
        float eta = -log(tan((M_PI/2 - atan(tanlam))/2));
        table_out[i] = eta * PF_ETAPHI_SCALE;
        // phi in -511,512 can only hold eta up to 2.23. else saturate for now
        if (eta * PF_ETAPHI_SCALE > (1<<(eta_T::width-1))-1) table_out[i] = (1<<(eta_T::width-1))-1;
    }
    return;
}

template<class tanlam_T, class eta_T>
void convert_eta(tanlam_T tanlam, eta_T &eta){
    // tanlam_T is ap_fixed<16,3>
    // eta_T is ap_int<10>

    // Initialize the lookup tables
#ifdef __HLS_SYN__
    bool initialized = false;
    etaphi_t eta_table[(1<<ETA_TAB_SIZE)];
#else 
    static bool initialized = false;
    static etaphi_t eta_table[(1<<ETA_TAB_SIZE)];
#endif
    if (!initialized) {
        init_eta_table<eta_T>(eta_table);
        initialized = true;
    }
    bool flip = false;
    if(tanlam<0){
        tanlam = -tanlam;
        flip=true;
    }
    utanlam_t utanlam = tanlam;

    ap_uint<ETA_TAB_SIZE> index;
    #pragma unroll
    for(int i=0; i<ETA_TAB_SIZE; i++){
        index[ETA_TAB_SIZE-1-i] = utanlam[utanlam_t::width-1-i]; //msb down to lowest
    }

    eta = eta_table[index];
    if(flip) eta = -eta;
}

