#ifndef __BITONIC_SORT_HPP__
#define __BITONIC_SORT_HPP__

#include "../../firmware/data.h"

// avoid Cosim Error
//#include "/opt/Xilinx/Vivado/2018.1/include/gmp.h"
//#include "/opt/Xilinx/Vivado/2018.1/include/mpfr.h"
#include "/data/Xilinx/Vivado/Vivado/2018.1/include/gmp.h"
#include "/data/Xilinx/Vivado/Vivado/2018.1/include/mpfr.h"

#include "ap_int.h"
#include "hls_stream.h"

#define BIT_LOW 0
#define BIT_HIGH 15

#define DATA_SIZE 64
#define DATA_NONZERO 47
#define DATA_W 64

void bitonic_sort(hls::stream<PFOutputObj > axis_in[], hls::stream<PFOutputObj > axis_out[]);


#endif
