# get the configuration
source config_hls_pfalgo3.tcl

# open the project, don't forget to reset
open_project -reset "proj3-mp7-fast"

set_top ${l1pfTopFunc}
add_files src/simple_pfalgo3.cpp
add_files -tb simple_pfalgo3_test.cpp  -cflags "-DTESTMP7 -DMP7_TOP_FUNC=${l1pfTopFunc} -DMP7_REF_FUNC=${l1pfRefFunc} -DMP7_VALIDATE=${l1pfValidate}"
add_files -tb simple_pfalgo3_ref.cpp
add_files -tb pattern_serializer.cpp
add_files -tb DiscretePFInputs.h    -cflags "-std=c++0x"
add_files -tb DiscretePFInputs_IO.h -cflags "-std=c++0x"
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
#cosim_design -trace_level all
export_design -format ip_catalog -vendor "cern-cms" -version ${l1pfIPVersion} -description "${l1pfTopFunc}"

# exit Vivado HLS
exit
