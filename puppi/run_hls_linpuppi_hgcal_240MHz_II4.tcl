source run_hls_linpuppi_hgcal_common.tcl
set hlsIPVersion 30.4.0

set cflags "-std=c++0x -DREG_${puppiReg} -DBOARD_${puppiBoard} -DHLS_pipeline_II=4" 
set kinds { "neutral" "charged" }

foreach kind $kinds {
    make_puppi ${puppiReg} ${puppiBoard} 3.0 "240MHz_II4" ${kind} ${cflags} ${hlsIPVersion}
}

exit
