# GlobalCorrelator_HLS
Vivado High Level Synthesis framework for Global Correlator

Project settings:
period: 4 ns
board: xc7vx690tffg1927-2

Import global_correlator.cpp as source.  Choose the function you want to test as the top level.  No need to manually include global_correlator.h.

Import global_correlator_test.cpp as test bench.  Input data is in tb_data.  You can compare output with expectations inside the test bench using standard C++.
