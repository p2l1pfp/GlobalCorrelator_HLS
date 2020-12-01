# get the configuration
#set pfBoard "none"
set pfBoard "VCU118"
set pfReg "HGCal"
set hlsIPVersion 25.4.0

set cflags "-std=c++0x -DREG_${pfReg} -DBOARD_${pfBoard} -DHLS_pipeline_II=4 -DL1PF_DSP_LATENCY3"

# open the project, don't forget to reset
open_project -reset "proj_pf${pfReg}_${pfBoard}_240MHz_II4"
if { $pfBoard == "none" } {
    set hlsTopFunc pfalgo2hgc
} else {
    set hlsTopFunc packed_pfalgo2hgc
}
set_top ${hlsTopFunc}
add_files firmware/pfalgo2hgc.cpp -cflags "${cflags}"
add_files -tb pfalgo2hgc_test.cpp  -cflags "${cflags}"
add_files -tb ref/pfalgo2hgc_ref.cpp  -cflags "${cflags}"
add_files -tb ref/pfalgo_common_ref.cpp  -cflags "${cflags}"
add_files -tb utils/pattern_serializer.cpp -cflags "${cflags}"
add_files -tb utils/test_utils.cpp -cflags "${cflags}"
add_files -tb data/TTbar_PU200_HGCal.dump

# reset the solution
open_solution -reset "solution"
set_part {xcvu9p-flga2104-2L-e}
create_clock -period 2.5 -name default

config_interface -trim_dangling_port
# do stuff
csim_design
csynth_design
#cosim_design -trace_level all
export_design -format ip_catalog -vendor "cern-cms" -version ${hlsIPVersion} -description "${hlsTopFunc}"

# exit Vivado HLS
exit
