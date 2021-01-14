source configIP.tcl

set cflags "-std=c++0x -DREG_${pfReg} -DBOARD_${pfBoard} -DROUTER_NOSTREAM -DNO_VALIDATE"
open_project -reset "project_csim_pf_puppi_tm18"

set_top ${hlsTopFunc}

set sample TTbar_PU200

add_files -tb regionizer_pf_puppi_test_tm18.cpp -cflags "${cflags}"
add_files -tb regionizer_ref.cpp -cflags "${cflags}"
add_files -tb utils/readMC.cpp -cflags "${cflags}"
add_files -tb tdemux/tdemux_ref.cpp   -cflags "${cflags}"
add_files -tb utils/tmux18_utils.cpp -cflags "${cflags}"
add_files -tb firmware/obj_unpackers.cpp -cflags "${cflags}"
add_files -tb utils/obj_packers.cpp -cflags "${cflags}"
add_files -tb ../utils/pattern_serializer.cpp -cflags "${cflags}"
add_files -tb ../utils/test_utils.cpp -cflags "${cflags}"
add_files -tb ../ref/pfalgo_common_ref.cpp   -cflags "${cflags}"
add_files -tb ../firmware/pfalgo2hgc.cpp   -cflags "${cflags}"
add_files -tb ../ref/pfalgo2hgc_ref.cpp   -cflags "${cflags}"
add_files -tb ../puppi/linpuppi_ref.cpp   -cflags "${cflags}"
add_files -tb data/caloDump_hgcal.${sample}.txt
add_files -tb data/trackDump_hgcalPos.${sample}.txt
add_files -tb data/muonDump_all.${sample}.txt
add_files -tb data/vertexDump_all.${sample}.txt

open_solution -reset "solution"
set_part {xcvu9p-flga2104-2L-e}
create_clock -period 2.5

csim_design
exit
