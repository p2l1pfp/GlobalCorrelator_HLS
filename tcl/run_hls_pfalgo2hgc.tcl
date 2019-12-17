# get the configuration
#source config_hls_fullpfalgo_mp7.tcl
#set pfBoard "none"
set pfBoard "VCU118"
set pfReg "HGCal"
set cflags "-std=c++0x -DREG_${pfReg} -DBOARD_${pfBoard}"

# open the project, don't forget to reset
open_project -reset "proj_pf${pfReg}_${pfBoard}"
if { $pfBoard == "none" } {
    set_top pfalgo2hgc
} else {
    set_top packed_pfalgo2hgc
}
add_files firmware/pfalgo2hgc.cpp -cflags "${cflags}"
add_files -tb pfalgo2hgc_test.cpp  -cflags "${cflags}"
add_files -tb pfalgo2hgc_ref.cpp  -cflags "${cflags}"
add_files -tb pfalgo_common_ref.cpp  -cflags "${cflags}"
add_files -tb utils/pattern_serializer.cpp -cflags "${cflags}"
add_files -tb utils/test_utils.cpp -cflags "${cflags}"
add_files -tb data/TTbar_PU200_HGCal.dump

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
