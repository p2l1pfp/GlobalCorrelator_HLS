set puppiReg "Barrel"
#set puppiReg "HGCal"

open_project -reset proj_linpuppi_${puppiReg}

set_top linpuppi_chs
add_files firmware/linpuppi.cpp  -cflags "-DREG_${puppiReg} -std=c++0x"
add_files -tb linpuppi_ref.cpp   -cflags "-DREG_${puppiReg} -std=c++0x"
add_files -tb ../utils/test_utils.cpp  -cflags "-DREG_${puppiReg}"
add_files -tb ../utils/pattern_serializer.cpp -cflags "-std=c++0x -DREG_${puppiReg}"
add_files -tb linpuppi_test.cpp   -cflags "-DREG_${puppiReg} -DTEST_PUPPI_NOCROP -DTEST_PT_CUT=120"
#add_files -tb linpuppi_test.cpp   -cflags "-DREG_${puppiReg} -DTEST_PT_CUT=120"
if { $puppiReg == "Barrel" } {
    add_files -tb ../pfalgo3_ref.cpp   -cflags "-DREG_${puppiReg} -std=c++0x"
} elseif { $puppiReg == "HGCal" } {
    add_files -tb ../pfalgo2hgc_ref.cpp   -cflags "-DREG_${puppiReg} -std=c++0x"
}
add_files -tb ../pfalgo_common_ref.cpp   -cflags "-DREG_${puppiReg} -std=c++0x"
add_files -tb ../data/TTbar_PU200_${puppiReg}.dump

# reset the solution
open_solution -reset "solution"
set_part {xcvu9p-flga2104-2L-e}
create_clock -period 3.0 -name default

config_interface -trim_dangling_port
# do stuff
csim_design
#csynth_design
#cosim_design -trace_level all
#export_design -format ip_catalog

# exit Vivado HLS
exit
