# move the configuration in here for ctp7, copy the mp7 project which just serializes the inputs

set l1pfAlgo "pfalgo3_full"
# PF implementation, with MP7 data wrapper around it
set l1pfTopFunc mp7wrapped_${l1pfAlgo}
# reference (non-wrapped) PF implementation
set l1pfRefFunc ${l1pfAlgo}_ref
# set to zero to turn off C validation (runs but does not check against the reference implementation)
set l1pfValidate 1
## version of the IP Core output
set l1pfIPVersion 2.0

# open the project, don't forget to reset
open_project -reset "proj3-ctp7-full"
set_top ${l1pfTopFunc}
add_files firmware/simple_fullpfalgo.cpp -cflags "-DTESTCTP7 -DHLS_pipeline_II=2"
add_files -tb simple_fullpfalgo_test.cpp  -cflags "-DTESTCTP7 -DHLS_pipeline_II=2 -DMP7_TOP_FUNC=${l1pfTopFunc} -DMP7_REF_FUNC=${l1pfRefFunc} -DCTP7_VALIDATE=${l1pfValidate}"
add_files -tb simple_fullpfalgo_ref.cpp -cflags "-DTESTCTP7"
add_files -tb utils/pattern_serializer.cpp -cflags "-DTESTCTP7"
add_files -tb utils/test_utils.cpp -cflags "-DTESTCTP7"
add_files -tb DiscretePFInputs.h    -cflags "-DTESTCTP7 -std=c++0x"
add_files -tb utils/DiscretePFInputs_IO.h -cflags "-DTESTCTP7 -std=c++0x"
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
