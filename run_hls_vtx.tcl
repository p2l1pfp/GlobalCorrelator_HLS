############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2016 Xilinx, Inc. All Rights Reserved.
############################################################

# open the project, don't forget to reset
open_project -reset proj1
set_top simple_vtx_hwopt
add_files firmware/simple_vtx.cpp
add_files -tb simple_vtx_test.cpp 
add_files -tb simple_vtx_ref.cpp

# reset the solution
open_solution -reset "solution1"
set_part {xc7k160tfbg484-2} -tool vivado
#set_part {xcku9p-ffve900-2-i-EVAL}
#set_part {xc7vx690tffg1927-2}
create_clock -period 5 -name default
#source "./nb1/solution1/directives.tcl"

# do stuff
csim_design
csynth_design
#cosim_design -trace_level all
#export_design -format ip_catalog

# exit Vivado HLS
exit
