# GlobalCorrelator_HLS
Vivado High Level Synthesis framework for Global Correlator

## Batch mode using Tcl script
```
vivado_hls -f run_hls_fullpfalgo_w_puppi.tcl #barrel
vivado_hls -f run_hls_hgcpfalgo_w_puppi.tcl #hgcal
vivado_hls -f run_hls_forwardpfalgo_test.tcl #hf
```
## Generating text files from .dump files
```
vivado_hls -f make_tmux_inputs.tcl > inputs.txt #generates text file of inputs for gen-x demonstrator or regionizer
vivado_hls -f make_tmux_outputs.tcl > outputs.txt #generates text file of expected outputs of regionizer (regionized inputs to pf+puppi)
vivado_hls -f make_tmux_layer1.tcl > layer1.txt #generates text file of expected outputs after regionizer and pf+puppi
```
n.b these scripts could be cleaned up, currently need some manual pruning of vivado messages etc.

## Generating custom .dump files
```
g++ test_dump.cc #modify test_dump.cc to alter content
./a.out
```
