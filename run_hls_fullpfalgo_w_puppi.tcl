# get the configuration
source config_hls_fullpfalgo_mp7.tcl

# open the project, don't forget to reset
open_project -reset "l1pfpuppi-resource-test"
set_top ${l1pfTopFunc}
#set_top pfalgo3_full
add_files firmware/simple_fullpfalgo.cpp -cflags "-DTESTMP7 -DHLS_pipeline_II=2"
add_files puppi/firmware/simple_puppi.cpp -cflags "-DTESTMP7 -DHLS_pipeline_II=2"
add_files bitonic-sort-48/hls/swap.hpp -cflags "-DTESTMP7"
add_files bitonic-sort-48/hls/swap.cpp -cflags "-DTESTMP7"
add_files bitonic-sort-48/hls/sorting_network_corr.hpp -cflags "-DTESTMP7"
add_files bitonic-sort-48/hls/sorting_network_corr.cpp -cflags "-DHLS_pipeline_II=2 -DTESTMP7"
add_files sorting/sort_group.cpp -cflags "-DHLS_pipeline_II=2 -DTESTMP7"
add_files sorting/sort_group.h -cflags "-DHLS_pipeline_II=2 -DTESTMP7"
add_files -tb fullpfalgo_test.cpp  -cflags "-DTESTMP7 -DHLS_pipeline_II=2 -DMP7_TOP_FUNC=${l1pfTopFunc} -DMP7_REF_FUNC=${l1pfRefFunc} -DMP7_VALIDATE=${l1pfValidate}"
add_files -tb simple_fullpfalgo_ref.cpp -cflags "-DTESTMP7"
add_files -tb utils/pattern_serializer.cpp -cflags "-DTESTMP7"
add_files -tb utils/test_utils.cpp -cflags "-DTESTMP7"
add_files -tb DiscretePFInputs.h    -cflags "-DTESTMP7 -std=c++0x"
add_files -tb utils/DiscretePFInputs_IO.h -cflags "-DTESTMP7 -std=c++0x"
add_files -tb data/regions_TTbar_PU140.dump
add_files -tb puppi/simple_puppi_ref.cpp -cflags "-DTESTMP7"
add_files -tb vertexing/simple_vtx_ref.cpp -cflags "-DTESTMP7"

# reset the solution
open_solution -reset "solution"
#set_part {xc7vx690tffg1927-2}
set_part {xcvu9p-flgb2104-2-i}
#set_part {xcvu9p-flga2577-2-e}
#set_part {xcvu13p-fhgb2104-2-i}
#create_clock -period 3.125 -name default
create_clock -period 2.857 -name default
set_clock_uncertainty 1.5

config_interface -trim_dangling_port
# do stuff
csim_design
csynth_design
cosim_design -trace_level all
export_design -format ip_catalog -vendor "cern-cms" -version ${l1pfIPVersion} -description "${l1pfTopFunc}"
#export_design -flow syn -format ip_catalog -vendor "cern-cms" -version ${l1pfIPVersion} -description "${l1pfTopFunc}"

# exit Vivado HLS
exit
