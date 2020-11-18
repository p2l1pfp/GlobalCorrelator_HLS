source configIP.tcl
# create a project
open_project -reset "project"
# specify the name of the function to synthetize
set_top ${hlsTopFunc}
# load source code for synthesis
add_files firmware/tdemux.cpp
# load source code for the testbench
add_files -tb testbench_tdemux.cpp
add_files -tb tdemux_ref.cpp

# create a solution (i.e. a hardware configuration for synthesis)
open_solution -reset "solution"
# set the FPGA (VU9P on VCU118), and a 360 MHz clock (2.78ns) with some extra margin
set_part {xcvu9p-flga2104-2L-e}
create_clock -period 2.2

# end here, so that we can then open the project interactively in the gui
csim_design
csynth_design
cosim_design -trace_level all
export_design -format ip_catalog -vendor "cern-cms" -version ${hlsIPVersion} -description ${hlsTopFunc}
exit
