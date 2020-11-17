set pfBoard "VCU118"
set pfReg "HGCal"
set hlsIPVersion 1.0

set cflags "-std=c++0x -DREG_${pfReg} -DBOARD_${pfBoard}"
set funcs { "unpack_track_3to2" "unpack_hgcal_3to1" "unpack_mu_3to12" }
foreach func ${funcs} {
    open_project -reset "project_${func}"
    set_top ${func}

    add_files firmware/obj_unpackers.cpp -cflags "${cflags}"
    add_files -tb regionizer_ref.cpp -cflags "${cflags}"
    add_files -tb readMC.cpp -cflags "${cflags}"
    add_files -tb ../utils/pattern_serializer.cpp -cflags "${cflags}"
    add_files -tb obj_unpackers_test.cpp -cflags "${cflags}"
    add_files -tb data/caloDump_hgcal.TTbar_PU200.txt
    add_files -tb data/trackDump_hgcalPos.TTbar_PU200.txt
    add_files -tb data/muonDump_all.TTbar_PU200.txt

    open_solution -reset "solution"
    set_part {xcvu9p-flga2104-2L-e}
    create_clock -period 2.2

    csim_design
    csynth_design
    cosim_design -trace_level all
}

exit
