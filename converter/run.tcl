# get the configuration
#source config.tcl

set mainAlgo pf_input_track_conv
set topFunc ${mainAlgo}_hw

# open the project, don't forget to reset
open_project -reset "track-conv-test"
set_top ${topFunc}
# design
add_files firmware/tk_input_converter.cpp
# reference
add_files -tb reference/tk_input_converter_ref.cpp
# tests
add_files -tb tests/test.cpp
add_files -tb tests/test_tanlambda_to_eta.cpp
add_files -tb tests/test_tk_resolution.cpp
add_files -tb tests/test_prop_tanlambda.cpp
add_files -tb tests/test_prop_phi.cpp
add_files -tb tests/test_pt_inversion.cpp
add_files -tb tests/test_pack_l1tk.cpp
add_files -tb tests/test_pack_pf.cpp

# reset the solution
open_solution -reset "solution"
set_part {xcvu9p-flgb2104-2-i}
create_clock -period 3.125 -name default
set_clock_uncertainty 1.5

config_interface -trim_dangling_port
# do stuff
csim_design
# csynth_design
# cosim_design -trace_level all
# export_design -format ip_catalog -vendor "cern-cms" -version 0.1 -description "${mainAlgo}"

# exit Vivado HLS
exit
