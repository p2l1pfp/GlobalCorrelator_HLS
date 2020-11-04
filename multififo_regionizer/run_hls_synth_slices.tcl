source configIP.tcl

set cflags "-std=c++0x -DREG_${pfReg} -DBOARD_${pfBoard} -DROUTER_NOMUX"
set dets { "tk" "calo" "mu" }
foreach det ${dets} {
    if {${det} == "tk"} {
        set slices { "input" "fifo" "merge2" "merge3" "output" }
    } elseif {${det} == "calo"} {
        set slices { "input" "fifo" "merge2" "merge4" "merge" "output" }
    } else {
        set slices { "input" "fifo" "merge" "output" }
    }

    foreach slice ${slices} {
        open_project -reset "project_synth_${det}_${slice}"
        set hlsTopFunc ${det}_router_${slice}_slice
        set_top ${hlsTopFunc}

        add_files firmware/calo_regionizer.cpp -cflags "${cflags}"
        add_files firmware/tk_regionizer.cpp -cflags "${cflags}"
        add_files firmware/mu_regionizer.cpp -cflags "${cflags}"
        add_files -tb regionizer_ref.cpp -cflags "${cflags}"
        add_files -tb readMC.cpp -cflags "${cflags}"
        add_files -tb regionizer_test.cpp -cflags "${cflags}"
        add_files -tb ../utils/pattern_serializer.cpp -cflags "${cflags}"
        add_files -tb ../utils/test_utils.cpp -cflags "${cflags}"
        add_files -tb data/caloDump_hgcal.TTbar_PU200.txt
        add_files -tb data/trackDump_hgcalPos.TTbar_PU200.txt
        add_files -tb data/muonDump_all.TTbar_PU200.txt

        open_solution -reset "solution"
        set_part {xcvu9p-flga2104-2L-e}
        create_clock -period 2.5

        csim_design
        csynth_design
        #cosim_design -trace_level all
    }
}

exit
