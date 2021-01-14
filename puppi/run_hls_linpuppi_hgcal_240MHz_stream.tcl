source run_hls_linpuppi_hgcal_common.tcl
set hlsIPVersion 25.1.0

set cflags "-std=c++0x -DREG_${puppiReg} -DBOARD_${puppiBoard} -DHLS_pipeline_II=1 -DLINPUPPI_DR2_LATENCY3" 

set kinds { "stream_prep" "stream_one" "stream_chs"  }

foreach kind $kinds {
    make_puppi ${puppiReg} ${puppiBoard} 2.5 "240MHz" ${kind} ${cflags} ${hlsIPVersion}
}

exit
