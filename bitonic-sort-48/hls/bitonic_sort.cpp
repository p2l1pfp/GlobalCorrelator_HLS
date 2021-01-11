#include "bitonic_sort.hpp"
#include "sorting_network.hpp"

void input(PFOutputObj datas[DATA_SIZE], hls::stream<PFOutputObj > axis_in[DATA_SIZE]) {
	// input data from stream
	//for (int i = 0; i < DATA_SIZE; ++i) {
	for (int i = 0; i < DATA_NONZERO; ++i) {
	#pragma HLS unroll
		datas[i] = axis_in[i].read();
	}
	for (int i = DATA_NONZERO; i < DATA_SIZE; ++i) {
	#pragma HLS unroll
		datas[i].hwPt = 0;
	}
}


void output(PFOutputObj datas[DATA_SIZE], hls::stream<PFOutputObj > axis_out[DATA_SIZE]) {
	// output data to stream
	for (int i = 0; i < DATA_SIZE; ++i) {
	#pragma HLS unroll
		axis_out[i].write(datas[i]);
	}
}


void bitonic_sort(hls::stream<PFOutputObj > axis_in[DATA_SIZE], hls::stream<PFOutputObj > axis_out[DATA_SIZE]) {
#pragma HLS interface ap_ctrl_hs port=return
#pragma HLS interface axis port=axis_in register
#pragma HLS interface axis port=axis_out register
#pragma HLS dataflow

	PFOutputObj datas[DATA_SIZE];
	#pragma HLS ARRAY_RESHAPE variable=datas complete dim=1


	input(datas, axis_in);

	sorting_network(datas);

	output(datas, axis_out);

}
