#include "obj_unpackers.h"

void unpack_track( ap_uint<96>   in,  bool   in_valid, 
                   ap_uint<64> & out, bool & out_valid)
{
    #pragma HLS pipeline II=1
    #pragma HLS latency min=1

    out_valid = in_valid;
    out(63, 0) = in(63, 0); // TODO check if it should instead be (95,32)
}


void unpack_track_3to2( ap_uint<64>   in1,   bool   in1_valid, 
                        ap_uint<64>   in2,   bool   in2_valid, 
                        ap_uint<64>   in3,   bool   in3_valid, 
                        ap_uint<64> & out1,  bool & out1_valid,
                        ap_uint<64> & out2,  bool & out2_valid)
{
    #pragma HLS pipeline II=1
    #pragma HLS latency min=1

    ap_uint<96> w1, w2;
    w1(95,32) = in1(63, 0);
    w1(31, 0) = in2(63,32);
    w2(95,64) = in2(31, 0);
    w2(63, 0) = in3(63, 0);
    unpack_track(w1, in1_valid && in2_valid, out1, out1_valid);
    unpack_track(w2, in2_valid && in3_valid, out2, out2_valid);
}



void unpack_hgcal( ap_uint<128>   in, bool   in_valid, 
                   ap_uint<64> & out, bool & out_valid)
{
    #pragma HLS pipeline II=1
    #pragma HLS latency min=1

    out_valid = in_valid;
    out(63, 0) = in(63, 0); // TODO check if it should instead be (128,64)
}



void unpack_hgcal_3to1( ap_uint<64>   in1,   bool   in1_valid, 
                        ap_uint<64>   in2,   bool   in2_valid, 
                        ap_uint<64>   in3,   bool   in3_valid, 
                        ap_uint<64> & out1,  bool & out1_valid) 
{
    #pragma HLS pipeline II=1
    #pragma HLS latency min=1

    ap_uint<128> w;
    w(127,64) = in1;
    w( 63, 0) = in2;
    unpack_hgcal(w, in1_valid && in2_valid, out1, out1_valid);
}

void unpack_mu( ap_uint<128>   in, bool   in_valid, 
                ap_uint<64> & out, bool & out_valid) 
{
    #pragma HLS pipeline II=1
    #pragma HLS latency min=1

    out_valid = in_valid;
    out(63, 0) = in(63, 0); // TODO check if it should instead be (128,64)
}

void unpack_mu_3to12(ap_uint<64>   in1,   bool   in1_valid, 
                     ap_uint<64>   in2,   bool   in2_valid, 
                     ap_uint<64>   in3,   bool   in3_valid, 
                     ap_uint<64> & out1,  bool & out1_valid,
                     ap_uint<64> & out2,  bool & out2_valid)
{
    #pragma HLS pipeline II=1
    #pragma HLS latency min=1

    static ap_uint<64> queue = 0; 
    static bool queue_valid = 0;
    if (queue_valid) {
        ap_uint<128> w1, w2;
        w1(127,64) = queue(63, 0);
        w1( 63, 0) =   in1(63, 0);
        w2(127,64) =   in2(63, 0);
        w2( 63, 0) =   in3(63, 0);
        unpack_mu(w1, in1_valid,              out1, out1_valid);
        unpack_mu(w2, in2_valid && in3_valid, out2, out2_valid);
        queue_valid = 0;
    } else {
        ap_uint<128> w1;
        w1(127,64) =   in1(63, 0);
        w1( 63, 0) =   in2(63, 0);
        queue       = in3;
        queue_valid = in3_valid;
        unpack_mu(w1, in1_valid && in2_valid, out1, out1_valid);
        out2 = 0;
        out2_valid = false;
    }

}
