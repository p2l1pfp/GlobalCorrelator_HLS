set puppiReg "HGCalNoTK"
#set puppiReg "HF"
#set puppiBoard "none"
set puppiBoard "VCU118"
set cflags "-std=c++0x -DREG_${puppiReg} -DBOARD_${puppiBoard}" 

open_project -reset proj_linpuppi_${puppiReg}_${puppiBoard}

if { $puppiBoard == "none" } {
    set_top fwdlinpuppiNoCrop
} else {
    set_top packed_fwdlinpuppiNoCrop

}

add_files firmware/linpuppi.cpp  -cflags "${cflags}"
add_files -tb linpuppi_ref.cpp   -cflags "${cflags}"
add_files -tb ../utils/test_utils.cpp  -cflags "${cflags}"
add_files -tb ../utils/pattern_serializer.cpp -cflags "${cflags}"
add_files -tb fwlinpuppi_test.cpp   -cflags "${cflags} -DTEST_PUPPI_NOCROP -DTEST_PT_CUT=120"
#add_files -tb fwlinpuppi_test.cpp   -cflags "${cflags}  -DTEST_PT_CUT=120"
#add_files -tb ../data/TTbar_PU200_${puppiReg}.dump
add_files -tb ../data/VBFHToBB_PU200_${puppiReg}.dump

# reset the solution
open_solution -reset "solution"
set_part {xcvu9p-flga2104-2L-e}
create_clock -period 3.0 -name default

config_interface -trim_dangling_port
# do stuff
csim_design
csynth_design
#cosim_design -trace_level all
#export_design -format ip_catalog

# exit Vivado HLS
exit
