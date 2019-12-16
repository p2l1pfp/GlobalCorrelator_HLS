# get the configuration
#source config_hls_fullpfalgo_mp7.tcl

# open the project, don't forget to reset
open_project -reset "proj_pfbarrel"
set_top pfalgo3
add_files firmware/pfalgo3.cpp -cflags "-std=c++0x -DREG_Barrel"
add_files -tb pfalgo3_test.cpp -cflags "-std=c++0x -DREG_Barrel"
add_files -tb pfalgo3_ref.cpp  -cflags "-std=c++0x -DREG_Barrel"
add_files -tb pfalgo_common_ref.cpp  -cflags "-std=c++0x -DREG_Barrel"
add_files -tb utils/pattern_serializer.cpp -cflags "-std=c++0x -DREG_Barrel"
add_files -tb utils/test_utils.cpp -cflags "-std=c++0x -DREG_Barrel"
add_files -tb data/TTbar_PU200_Barrel.dump

# reset the solution
open_solution -reset "solution"
set_part {xcvu9p-flga2104-2L-e}
create_clock -period 3.0 -name default

config_interface -trim_dangling_port
# do stuff
csim_design
csynth_design
cosim_design -trace_level all
#export_design -format ip_catalog -vendor "cern-cms" -version ${l1pfIPVersion} -description "${l1pfTopFunc}"

# exit Vivado HLS
exit
