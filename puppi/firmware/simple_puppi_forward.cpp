#include "simple_puppi_forward.h"
#include "../../firmware/simple_forwardpfalgo.h"
#include <string.h>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

weight_t puppiweight(int iWeight){
  static weight_t _table[PUPPI_TABLE_SIZE];
  lut_puppiweight_init<weight_t,PUPPI_TABLE_SIZE>(_table);
  return _table[iWeight];
}

/*int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}*/


void _lut_shift15_invert_init(ap_uint<16> _table[512]) { // returns 2^15 / x
	_table[0] = 32768; // this is 2^15
	for (int i = 1; i <= 511; ++i) {
		_table[i] = (32768 / i);
	}
}
int _lut_shift15_divide(ap_uint<17> num, ap_uint<9> den) { // returns (num * 2^15) / den
	ap_uint<16> _table[512];
	_lut_shift15_invert_init(_table);
	return (num * _table[den]);
}

template <int DR2MAX>
int get_sum(int dr2[NNEUTRALS], ap_uint<17> pt2_shift[NNEUTRALS]) {

    #pragma HLS ARRAY_PARTITION variable=pt2_shift complete
    #pragma HLS ARRAY_PARTITION variable=dr2 complete

    int sum = 0;
    for (int it = 0; it < NNEUTRALS; ++it) {
        #pragma HLS LOOP UNROLL
        if (dr2[it] <= DR2MAX) {
            ap_uint<9> dr2short = dr2[it] >> 5; // why?
            int term = _lut_shift15_divide(pt2_shift[it], dr2short);
            sum += term;
        }
    }
    return sum;
}

//void simple_puppi_forward_hw(PFNeutralObj pfallne[NNEUTRALS], pt_t ptpuppi[NNEUTRALS]) {
void simple_puppi_forward_hw(PFNeutralObj pfallne[NNEUTRALS], pt_t ptpuppi[NNEUTRALS], em2calo_dr_t drvals[NPHOTON][NSELCALO]) {

    const int DR2MAX = 8404; // 0.4 cone

    #pragma HLS INTERFACE ap_none port=pfallne
    #pragma HLS INTERFACE ap_none port=ptpuppi
    #pragma HLS INTERFACE ap_none port=drvals

    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS ARRAY_PARTITION variable=ptpuppi complete
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0
    #pragma HLS pipeline II=HLS_pipeline_II

    ap_uint<32> eToAlphas[NNEUTRALS];
    ap_uint<8> weights[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=eToAlphas complete
    #pragma HLS ARRAY_PARTITION variable=weights complete    

    ap_uint<17> pt2_shift[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=pt2_shift complete
    int dr2[NNEUTRALS][NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=dr2 complete dim=0
    int dr2flat_lo[((NPHOTON)*(NPHOTON-1))/2];
    #pragma HLS ARRAY_PARTITION variable=dr2flat_lo complete dim=0
    int dr2flat_hi[((NSELCALO)*(NSELCALO-1))/2];
    #pragma HLS ARRAY_PARTITION variable=dr2flat_hi complete dim=0
    for (int in = 0; in < NNEUTRALS; ++in){ 
        eToAlphas[in] = 0; weights[in] = 0; 
        int mypt2 = (pfallne[in].hwPt*pfallne[in].hwPt) >> 5;
        pt2_shift[in] = (mypt2 < 131071 ? mypt2 : 131071);    
    }

    //for (int in = 0; in < NNEUTRALS; ++in){ 
    //    for (int it = 0; it < NNEUTRALS; ++it) {
    //        dr2[in][it] = dr2_int(pfallne[it].hwEta, pfallne[it].hwPhi, pfallne[in].hwEta, pfallne[in].hwPhi);
    //        if (in==it) dr2[in][it] = DR2MAX+1;
    //    }
    //}

    for (int in = 1; in < NPHOTON; ++in){ 
        for (int it = 0; it < in; ++it) {
            dr2flat_lo[((in*(in-1))/2)+it] = dr2_int(pfallne[it].hwEta, pfallne[it].hwPhi, pfallne[in].hwEta, pfallne[in].hwPhi);
        }
    }
    for (int in = 1; in < NSELCALO; ++in){ 
        for (int it = 0; it < in; ++it) {
            dr2flat_hi[((in*(in-1))/2)+it] = dr2_int(pfallne[it+NPHOTON].hwEta, pfallne[it+NPHOTON].hwPhi, pfallne[in+NPHOTON].hwEta, pfallne[in+NPHOTON].hwPhi);
        }
    }
    for (int in = 0; in < NPHOTON; ++in){ 
        for (int it = 0; it < NPHOTON; ++it) {
            if (it < in) dr2[in][it] = dr2flat_lo[((in*(in-1))/2)+it];
            else if (in==it) dr2[in][it] = DR2MAX+1;
            else dr2[in][it] = dr2flat_lo[((it*(it-1))/2)+in];
        }
    }
    for (int in = 0; in < NPHOTON; ++in){ 
        for (int it = 0; it < NSELCALO; ++it) {
            dr2[in][it+NPHOTON] = drvals[in][it];
            dr2[it+NPHOTON][in] = dr2[in][it+NPHOTON];
        }
    }
    for (int in = 0; in < NSELCALO; ++in){ 
        for (int it = 0; it < NSELCALO; ++it) {
            if (it < in) dr2[in+NPHOTON][it+NPHOTON] = dr2flat_hi[((in*(in-1))/2)+it];
            else if (in==it) dr2[in+NPHOTON][it+NPHOTON] = DR2MAX+1;
            else dr2[in+NPHOTON][it+NPHOTON] = dr2flat_hi[((it*(it-1))/2)+in];
        }
    }
    
    for (int in = 0; in < NNEUTRALS; ++in) {
        /*int sum = 0;
        for (int it = 0; it < NNEUTRALS; ++it) {
            int dr2temp = 0;
            if (in==it) continue;
            else if (it<in) {dr2temp = dr2[((in*(in-1))/2)+it];
            } else {dr2temp = dr2[((it*(it-1))/2)+in]; }
            //int dr2temp = dr2_int(pfallne[it].hwEta, pfallne[it].hwPhi, pfallne_cpy[in].hwEta, pfallne_cpy[in].hwPhi);
            if (dr2temp <= DR2MAX) {
                //std::cout<<"HW:  "<<dr2temp<<"\t "<<in<<" "<<it<<std::endl;
                ap_uint<9> dr2short = dr2temp >> 5; // why?
                int term = _lut_shift15_divide(pt2_shift[it], dr2short);
                sum += term;
            }
        }*/

        int sum = get_sum<DR2MAX>(dr2[in], pt2_shift);
        eToAlphas[in] = sum >> 10;
    }
    

    for (int in = 0; in < NNEUTRALS; ++in) {
    	int ptnew = pfallne[in].hwPt;
        if (eToAlphas[in] <= 0){
            ptnew = 0; // < e^10 where that is the median
        } 
        else if (eToAlphas[in] < PUPPI_TABLE_SIZE) {
            weights[in] = puppiweight(eToAlphas[in]); //(int) puppiweight_table[index];
    	    ptnew = (ptnew * weights[in]) >> 8;
        }
    	ptpuppi[in] = (pt_t) ptnew;
    }

}


