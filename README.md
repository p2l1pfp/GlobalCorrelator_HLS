# GlobalCorrelator_HLS
Vivado High Level Synthesis framework for Global Correlator

Project Settings:
period: 4 ns,
board: xc7vx690tffg1927-2

# Batch mode using Tcl script
`vivado_hls -f run_hls.tcl`

# To open the project in the GUI
After creating the project using the tcl script (default name from tcl script is `proj0`)
`vivado_hls -p proj0`

# creating the project in the GUI 
n.b. this needs some modificiations to the files
Import simple_pflow.cpp as source.  Choose the function you want to test as the top level.  No need to manually include simple_pflow.h.
Import simple_pfow_test.cpp as test bench. You can compare output with expectations inside the test bench using standard C++.
