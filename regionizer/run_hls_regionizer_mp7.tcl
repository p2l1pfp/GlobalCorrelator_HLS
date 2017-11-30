############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2016 Xilinx, Inc. All Rights Reserved.
############################################################

# open the project, don't forget to reset
open_project -reset proj-regionizer-mp7
set_top merge_hadcalo
add_files firmware/regionizer.cpp -cflags "-DHLS_pipeline_II=1"
add_files -tb regionizer_test.cpp -cflags "-DMP7"
add_files -tb regionizer_ref.cpp
add_files -tb ../utils/pattern_serializer.cpp
add_files -tb ../utils/test_utils.cpp
add_files -tb ../DiscretePFInputs.h -cflags "-std=c++0x"
add_files -tb ../utils/DiscretePFInputs_IO.h -cflags "-std=c++0x"
add_files -tb data/barrel_sectors_1x12_TTbar_PU140.dump

# reset the solution
open_solution -reset "solution"
set_part {xc7vx690tffg1927-2}
#set_part {xcku9p-ffve900-2-i-EVAL}
#set_part {xcku115-flvf1924-2-i}
create_clock -period 4.16667 -name default
set_clock_uncertainty 1.0

config_interface -trim_dangling_port
# do stuff
csim_design
#csynth_design
#cosim_design -trace_level all
#export_design -format ip_catalog

# exit Vivado HLS
exit
