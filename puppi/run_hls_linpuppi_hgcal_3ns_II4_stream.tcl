set puppiReg "HGCal"
set puppiBoard "VCU118"
set hlsIPVersion 30.4.0

set cflags "-std=c++0x -DREG_${puppiReg} -DBOARD_${puppiBoard} -DHLS_pipeline_II=4" 

set kinds { "stream_prep" "stream_one" "stream_chs"  }

foreach kind $kinds {
    open_project -reset "proj_linpuppi_${puppiReg}_${puppiBoard}_3ns_II4_${kind}"

    set test "TEST_PUPPI_NOCROP"

    if { $puppiBoard == "none" } {
        if { $kind == "neutral" } {
            set hlsTopFunc linpuppiNoCrop
            #set hlsTopFunc linpuppi
        } else {
            set hlsTopFunc linpuppi_chs
        }
    } else {
        if { $kind == "stream_prep" } {
            set test "TEST_PUPPI_STREAM"
            set hlsTopFunc packed_linpuppi_prepare_track
        } elseif { $kind == "stream_one" } {
            set test "TEST_PUPPI_STREAM"
            set hlsTopFunc packed_linpuppi_one
        } elseif { $kind == "stream_chs" } {
            set test "TEST_PUPPI_STREAM"
            set hlsTopFunc packed_linpuppi_chs_one
        } elseif { $kind == "neutral" } {
            set hlsTopFunc packed_linpuppiNoCrop
            #set hlsTopFunc packed_linpuppi
        } else {
            set hlsTopFunc packed_linpuppi_chs
        }
    }

    set_top ${hlsTopFunc}
    add_files firmware/linpuppi.cpp  -cflags "${cflags}"
    add_files -tb linpuppi_ref.cpp   -cflags "${cflags}"
    add_files -tb ../utils/test_utils.cpp  -cflags "${cflags}"
    add_files -tb ../utils/pattern_serializer.cpp -cflags "${cflags}"
    add_files -tb linpuppi_test.cpp   -cflags "${cflags} -D${test} -DTEST_PT_CUT=80" 
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
    create_clock -period 3.0 -name default

    config_interface -trim_dangling_port
    if [ string match "stream*" $kind ] {
        config_rtl -reset none  
    }
    # do stuff
    csim_design
    csynth_design
    #cosim_design -trace_level all
    export_design -format ip_catalog -vendor "cern-cms" -version ${hlsIPVersion} -description "${hlsTopFunc}"

}
# exit Vivado HLS
exit
