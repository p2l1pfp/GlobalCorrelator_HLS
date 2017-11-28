# get the configuration
source config_hls_fullpfalgo_mp7.tcl

# open the project, don't forget to reset
open_project -reset "proj3-mp7-full"
set_top ${l1pfTopFunc}
add_files firmware/simple_fullpfalgo.cpp -cflags "-DTESTMP7 -DHLS_pipeline_II=2"
add_files -tb simple_fullpfalgo_test.cpp  -cflags "-DTESTMP7 -DHLS_pipeline_II=2 -DMP7_TOP_FUNC=${l1pfTopFunc} -DMP7_REF_FUNC=${l1pfRefFunc} -DMP7_VALIDATE=${l1pfValidate}"
add_files -tb simple_fullpfalgo_ref.cpp -cflags "-DTESTMP7"
add_files -tb utils/pattern_serializer.cpp -cflags "-DTESTMP7"
add_files -tb utils/test_utils.cpp -cflags "-DTESTMP7"
add_files -tb DiscretePFInputs.h    -cflags "-DTESTMP7 -std=c++0x"
add_files -tb utils/DiscretePFInputs_IO.h -cflags "-DTESTMP7 -std=c++0x"
add_files -tb data/regions_TTbar_PU140.dump

# reset the solution
open_solution -reset "solution"
set_part {xc7vx690tffg1927-2}
create_clock -period 4.16667 -name default
set_clock_uncertainty 1.5

config_interface -trim_dangling_port
# do stuff
csim_design
csynth_design
cosim_design -trace_level all
export_design -format ip_catalog -vendor "cern-cms" -version ${l1pfIPVersion} -description "${l1pfTopFunc}"

# exit Vivado HLS
exit
