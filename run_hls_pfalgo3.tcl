############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2016 Xilinx, Inc. All Rights Reserved.
############################################################

# open the project, don't forget to reset
open_project -reset proj3-mp7-fast
#open_project -reset proj3
#set_top pfalgo3_full
#set_top mp7wrapped_pfalgo3_full
set_top mp7wrapped_pfalgo3_fast
add_files src/simple_pfalgo3.cpp
#add_files -tb simple_pfalgo3_test.cpp  -cflags "-DTESTFULL"
add_files -tb simple_pfalgo3_test.cpp  -cflags "-DTESTMP7 -DTESTMP7FAST"
add_files -tb simple_pfalgo3_ref.cpp
add_files -tb pattern_serializer.cpp
add_files -tb DiscretePFInputs.h -cflags "-std=c++0x"
add_files -tb DiscretePFInputs_IO.h -cflags "-std=c++0x"
add_files -tb data/regions_TTbar_PU140.dump

# reset the solution
open_solution -reset "solution1"
#set_part {xcku9p-ffve900-2-i-EVAL}
set_part {xc7vx690tffg1927-2}
#set_part {xcku5p-sfvb784-3-e}
#set_part {xcku115-flvf1924-2-i}
create_clock -period 4.16667 -name default
set_clock_uncertainty 1.5
#source "./nb1/solution1/directives.tcl"

config_interface -trim_dangling_port
# do stuff
csim_design
csynth_design
#cosim_design -trace_level all
export_design -format ip_catalog -vendor "cern-cms" -version 1.3 -description "mp7wrapped_pfalgo3_fast"

# exit Vivado HLS
exit
