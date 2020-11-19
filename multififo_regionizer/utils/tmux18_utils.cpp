#include "tmux18_utils.h"

ap_uint<65> DelayQueue::operator()(const ap_uint<65> & in) 
{
    ap_uint<65> ret = data_[ptr_];
    data_[ptr_] = in;
    ptr_++; 
    if (ptr_ == n_) ptr_ = 0;
    return ret;
}

void DelayQueue::operator()(const ap_uint<64> & in,  const bool & in_valid,
                                  ap_uint<64> & out,       bool & out_valid) 
{
    ap_uint<65> in65  = in; in65[64] = in_valid;
    ap_uint<65> out65 = (*this)(in65);
    out = out65(63,0); out_valid = out65[64];
}
