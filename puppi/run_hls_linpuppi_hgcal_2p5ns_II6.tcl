#set puppiReg "Barrel"
set puppiReg "HGCal"
set puppiBoard "none"
#set puppiBoard "VCU118"
set cflags "-std=c++0x -DREG_${puppiReg} -DBOARD_${puppiBoard} -DHLS_pipeline_II=6" 

open_project -reset "proj_linpuppi_${puppiReg}_${puppiBoard}_2.5ns_II6_30"

if { $puppiBoard == "none" } {
    set_top linpuppiNoCrop
    #set_top linpuppi
    #set_top linpuppi_chs
} else {
    set_top packed_linpuppiNoCrop
    #set_top packed_linpuppi
    #set_top packed_linpuppi_chs
}
add_files firmware/linpuppi.cpp  -cflags "${cflags}"
add_files -tb linpuppi_ref.cpp   -cflags "${cflags}"
add_files -tb ../utils/test_utils.cpp  -cflags "${cflags}"
add_files -tb ../utils/pattern_serializer.cpp -cflags "${cflags}"
add_files -tb linpuppi_test.cpp   -cflags "${cflags} -DTEST_PUPPI_NOCROP -DTEST_PT_CUT=80" 
#add_files -tb linpuppi_test.cpp   -cflags "${cflags} -DTEST_PT_CUT=80"
if { $puppiReg == "Barrel" } {
    add_files -tb ../ref/pfalgo3_ref.cpp   -cflags "${cflags}"
} elseif { $puppiReg == "HGCal" } {
    add_files -tb ../ref/pfalgo2hgc_ref.cpp   -cflags "${cflags}"
}
add_files -tb ../ref/pfalgo_common_ref.cpp   -cflags "${cflags}"
add_files -tb ../data/TTbar_PU200_${puppiReg}.dump

# reset the solution
open_solution -reset "solution"
set_part {xcvu9p-flga2104-2L-e}
create_clock -period 2.5 -name default

config_interface -trim_dangling_port
# do stuff
csim_design
csynth_design
#cosim_design -trace_level all
#export_design -format ip_catalog

# exit Vivado HLS
exit
