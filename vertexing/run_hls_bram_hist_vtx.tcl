############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2016 Xilinx, Inc. All Rights Reserved.
############################################################

# open the project, don't forget to reset
open_project -reset proj_bhv
#set_top bhv_add_track
set_top bhv_find_pv
add_files firmware/bram_hist_vtx.cpp
add_files -tb bram_hist_vtx_ref.cpp -cflags "-std=c++0x"
add_files -tb bram_hist_vtx_test.cpp  -cflags "-std=c++0x"
add_files -tb ../utils/DiscretePFInputs_IO.h -cflags "-std=c++0x"
add_files -tb ../utils/pattern_serializer.cpp
add_files -tb ../regionizer/data/barrel_sectors_1x12_TTbar_PU140.dump
add_files -tb ../regionizer/data/barrel_alltracks_sectors_1x12_TTbar_PU140.dump

# reset the solution
open_solution -reset "solution1"
#set_part {xc7k160tfbg484-2} -tool vivado
#set_part {xcku9p-ffve900-2-i-EVAL}
#set_part {xc7vx690tffg1927-2}
set_part {xcku115-flvf1924-2-i}
create_clock -period 5 -name default

# do stuff
csim_design
#csynth_design
#cosim_design -trace_level all
#export_design -format ip_catalog

# exit Vivado HLS
exit
