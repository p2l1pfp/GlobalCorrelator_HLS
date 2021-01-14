source run_hls_linpuppi_hgcal_common.tcl
set hlsIPVersion 22.6

set cflags "-std=c++0x -DREG_${puppiReg} -DBOARD_${puppiBoard} -DHLS_pipeline_II=6 -DLINPUPPI_DR2_LATENCY3" 

set kinds { "neutral" "charged" }

foreach kind $kinds {
    make_puppi ${puppiReg} ${puppiBoard} 2.2 "360MHz_II6" ${kind} ${cflags} ${hlsIPVersion}
}

exit
