#ifndef multififo_regionizer_unpackers_h
#define multififo_regionizer_unpackers_h

#include <ap_int.h>

void unpack_track( ap_uint<96>   in,  bool   in_valid, 
                   ap_uint<64> & out, bool & out_valid);

void unpack_track_3to2( ap_uint<64>   in1,   bool   in1_valid, 
                        ap_uint<64>   in2,   bool   in2_valid, 
                        ap_uint<64>   in3,   bool   in3_valid, 
                        ap_uint<64> & out1,  bool & out1_valid,
                        ap_uint<64> & out2,  bool & out2_valid);


void unpack_hgcal( ap_uint<128>   in, bool   in_valid, 
                   ap_uint<64> & out, bool & out_valid);


void unpack_hgcal_3to1( ap_uint<64>   in1,   bool   in1_valid, 
                        ap_uint<64>   in2,   bool   in2_valid, 
                        ap_uint<64>   in3,   bool   in3_valid, 
                        ap_uint<64> & out1,  bool & out1_valid);

void unpack_mu( ap_uint<128>   in, bool   in_valid, 
                ap_uint<64> & out, bool & out_valid);

void unpack_mu_3to12(ap_uint<64>   in1,   bool   in1_valid, 
                     ap_uint<64>   in2,   bool   in2_valid, 
                     ap_uint<64>   in3,   bool   in3_valid, 
                     ap_uint<64> & out1,  bool & out1_valid,
                     ap_uint<64> & out2,  bool & out2_valid);

#endif
