############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2016 Xilinx, Inc. All Rights Reserved.
############################################################

# open the project, don't forget to reset
open_project -reset proj_pfpuppi_forward
# set_top compute_puppi_weight_hw
set_top simple_puppi_forward_hw
add_files firmware/simple_puppi_forward.cpp -cflags "-DTESTMP7 -DHLS_pipeline_II=2"
#add_files ../firmware/simple_fullpfalgo.cpp 
add_files -tb simple_pfpuppi_forward_test.cpp -cflags "-DTESTMP7"
add_files -tb ../simple_forwardpfalgo_ref.cpp -cflags "-DTESTMP7"
add_files -tb simple_puppi_forward_ref.cpp -cflags "-DTESTMP7"
add_files -tb ../utils/pattern_serializer.cpp -cflags "-std=c++0x"
add_files -tb ../utils/test_utils.cpp -cflags "-DTESTMP7"
add_files -tb ../data/regions_TTbar_PU140.dump

# reset the solution
open_solution -reset "solution1"
#set_part {xcku9p-ffve900-2-i-EVAL}
#set_part {xc7vx690tffg1927-2}
set_part {xcvu9p-flgb2104-2-i}
#set_part {xcku5p-sfvb784-3-e}
#set_part {xcku115-flvf1924-2-i}
create_clock -period 4.16667 -name default
set_clock_uncertainty 1.0
#source "./nb1/solution1/directives.tcl"
config_array_partition -maximum_size 8192

config_interface -trim_dangling_port
# do stuff
csim_design
csynth_design
#cosim_design -trace_level all
#export_design -format ip_catalog

cosim_design -trace_level all
export_design -format ip_catalog -vendor "cern-cms" -version 1.0 -description "pfpuppi"

# exit Vivado HLS
exit
