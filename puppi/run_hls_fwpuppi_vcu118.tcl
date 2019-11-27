# open the project, don't forget to reset
open_project -reset proj_fwdpuppi

# set_top compute_puppi_weight_hw
#set_top fwdlinpuppiPt_hw
#set_top fwdlinpuppiNoCrop_hw
set_top fwdlinpuppi_hw
add_files firmware/linpuppi.cpp  -cflags "-DREG_HGCALNOTK"
add_files -tb linpuppi_ref.cpp   -cflags "-DREG_HGCALNOTK"
add_files -tb ../utils/test_utils.cpp  -cflags "-DREG_HGCALNOTK"
#add_files -tb fwlinpuppi_test.cpp   -cflags "-DREG_HGCALNOTK -DTEST_PUPPI_PT -DTEST_PT_CUT=120"
#add_files -tb fwlinpuppi_test.cpp   -cflags "-DREG_HGCALNOTK -DTEST_PUPPI_NOCROP -DTEST_PT_CUT=120"
add_files -tb fwlinpuppi_test.cpp   -cflags "-DREG_HGCALNOTK -DTEST_PT_CUT=120"
#add_files -tb ../data/TTbar_PU200_HGCalNoTK.dump
add_files -tb ../data/VBFHToBB_PU200_HGCalNoTK.dump

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
