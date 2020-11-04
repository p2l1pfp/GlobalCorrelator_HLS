#!/bin/bash


#VHDLS="<IMPL>/route_link2fifo.vhd <IMPL>/router_monolythic_fifos_data_V_0.vhd <IMPL>/router_monolythic.vhd phi_regionizer_tb.vhd"
FW="../hdl"
HLS="../../.."
if [[ "$1" == "hls_nomerge" ]]; then
    VHDLS="<IMPL>/route_link2fifo.vhd <IMPL>/router_nomerge_fifos_data_V_0.vhd <IMPL>/router_nomerge.vhd phi_regionizer_nomerge_tb.vhd"
    HLSPROJ="project_nomerge"
#elif [[ "$1" == "vhdl_nomerge" ]]; then
#    VHDLS="${FW}/regionizer_data.vhd ${FW}/tk_router_element.vhd ${FW}/tk_router.vhd ${FW}/rolling_fifo.vhd ${FW}/phi_regionizer_nomerge.vhd phi_regionizer_nomerge_vhdl_tb.vhd"
#    HLSPROJ="project_nomerge"
#elif [[ "$1" == "vhdl_m2" ]]; then
#    VHDLS="${FW}/regionizer_data.vhd ${FW}/tk_router_element.vhd ${FW}/tk_router.vhd ${FW}/rolling_fifo.vhd ${FW}/fifo_merge2.vhd ${FW}/phi_regionizer_m2.vhd phi_regionizer_m2_vhdl_tb.vhd"
#    HLSPROJ="project_m2_input"
#elif [[ "$1" == "hls_m2_slices" ]]; then
#    VHDLS="$VHDLS project_m2_input/solution/syn/vhdl/route_link2fifo.vhd  project_m2_input/solution/syn/vhdl/router_input_slice.vhd"
#    VHDLS="$VHDLS project_m2_fifo/solution/syn/vhdl/router_fifo_slice_fifos_data_V_0.vhd  project_m2_fifo/solution/syn/vhdl/router_fifo_slice.vhd"
#    VHDLS="$VHDLS project_m2_merge2/solution/syn/vhdl/router_m2_merge2_slice.vhd"
#    VHDLS="$VHDLS project_m2_output/solution/syn/vhdl/router_m2_output_slice.vhd"
#    VHDLS="$VHDLS ${FW}/regionizer_data_stdlogic.vhd ${FW}/phi_regionizer_m2_hls_slices.vhd phi_regionizer_m2_vhdl_tb.vhd"
#    HLSPROJ="project_m2_input"
#elif [[ "$1" == "hls_slices" ]]; then
#    VHDLS="$VHDLS project_full_input/solution/syn/vhdl/route_link2fifo.vhd  project_full_input/solution/syn/vhdl/router_input_slice.vhd"
#    VHDLS="$VHDLS project_full_fifo/solution/syn/vhdl/router_fifo_slice_fifos_data_V_0.vhd  project_full_fifo/solution/syn/vhdl/router_fifo_slice.vhd"
#    VHDLS="$VHDLS project_full_merge2/solution/syn/vhdl/router_merge2_slice.vhd"
#    VHDLS="$VHDLS project_full_merge3/solution/syn/vhdl/router_merge3_slice.vhd"
#    VHDLS="$VHDLS project_full_output/solution/syn/vhdl/router_full_output_slice.vhd"
#    VHDLS="$VHDLS ${FW}/regionizer_data_stdlogic.vhd ${FW}/phi_regionizer_hls_slices.vhd phi_regionizer_hls_tb.vhd"
#    HLSPROJ="project_full_input"
elif [[ "$1" == "vhdl-tk" ]]; then
    VHDLS="${FW}/regionizer_data.vhd ${FW}/tk_router_element.vhd ${FW}/tk_router.vhd ${FW}/rolling_fifo.vhd ${FW}/fifo_merge2_full.vhd ${FW}/fifo_merge3.vhd"
    VHDLS="${VHDLS} ${FW}/tk_regionizer.vhd tk_regionizer_vhdl_tb.vhd"
    HLSPROJ="project_csim"
    DET="tk"
#elif [[ "$1" == "vhdl_sort" ]]; then
#    VHDLS="${FW}/regionizer_data.vhd ${FW}/tk_router_element.vhd ${FW}/tk_router.vhd ${FW}/rolling_fifo.vhd ${FW}/fifo_merge2_full.vhd ${FW}/fifo_merge3.vhd ${FW}/stream_sort.vhd ${FW}/phi_regionizer_sort.vhd phi_regionizer_sorted_vhdl_tb.vhd"
#    HLSPROJ="project_mux"
#elif [[ "$1" == "vhdl_mux" ]]; then
#    VHDLS="${FW}/regionizer_data.vhd ${FW}/tk_router_element.vhd ${FW}/tk_router.vhd ${FW}/rolling_fifo.vhd ${FW}/fifo_merge2_full.vhd ${FW}/fifo_merge3.vhd ${FW}/stream_sort.vhd ${FW}/region_mux_stream.vhd ${FW}/phi_regionizer_mux.vhd phi_regionizer_mux_vhdl_tb.vhd"
#    HLSPROJ="project_mux"
#elif [[ "$1" == "calo_vhdl_nomerge" ]]; then
#    VHDLS="${FW}/regionizer_data.vhd ${FW}/calo_router.vhd ${FW}/rolling_fifo.vhd ${FW}/calo_phi_regionizer_nomerge.vhd calo_phi_regionizer_nomerge_vhdl_tb.vhd"
#    HLSPROJ="project_nomergeMC_calo"
elif [[ "$1" == "vhdl-calo" ]]; then
    VHDLS="${FW}/regionizer_data.vhd ${FW}/calo_router.vhd ${FW}/rolling_fifo.vhd ${FW}/fifo_merge2_full.vhd ${FW}/fifo_merge2.vhd "
    VHDLS="${VHDLS} ${FW}/calo_regionizer.vhd calo_regionizer_vhdl_tb.vhd"
    HLSPROJ="project_csim"
    DET="calo"
#elif [[ "$1" == "calo_hls_slices" ]]; then
#    VHDLS="$VHDLS project_MC_calo_input/solution/syn/vhdl/calo_router_input_slice.vhd"
#    VHDLS="$VHDLS project_MC_calo_fifo/solution/syn/vhdl/calo_router_fifo_slice_fifos_data_V_0.vhd"
#    VHDLS="$VHDLS project_MC_calo_fifo/solution/syn/vhdl/calo_router_fifo_slice.vhd"
#    VHDLS="$VHDLS project_MC_calo_merge2/solution/syn/vhdl/calo_router_merge2_slice.vhd"
#    VHDLS="$VHDLS project_MC_calo_merge4/solution/syn/vhdl/calo_router_merge4_slice.vhd"
#    VHDLS="$VHDLS project_MC_calo_merge/solution/syn/vhdl/calo_router_merge_slice.vhd"
#    VHDLS="$VHDLS project_MC_calo_output/solution/syn/vhdl/calo_router_full_output_slice.vhd"
#    VHDLS="$VHDLS ${FW}/regionizer_data_stdlogic.vhd ${FW}/calo_phi_regionizer_hls_slices.vhd calo_phi_regionizer_hls_slices_tb.vhd"
#    HLSPROJ="project_MC_calo_input"
elif [[ "$1" == "vhdl-mu" ]]; then
    VHDLS="${FW}/regionizer_data.vhd ${FW}/mu_router.vhd ${FW}/rolling_fifo.vhd ${FW}/fifo_merge2.vhd "
    VHDLS="${VHDLS} ${FW}/mu_regionizer.vhd mu_regionizer_vhdl_tb.vhd"
    HLSPROJ="project_csim"
    DET="mu"

fi


CSIM=$HLS/$HLSPROJ/solution/csim/build
if test -f $CSIM/input-${DET}.txt; then
    echo " ## Getting C simulation inputs from $CSIM";
    cp -v $CSIM/*-${DET}.txt .
else
    echo "Couldn't find C simulation inputs in $CSIM.";
    echo "Run vivado_hls in the parent directory before.";
    exit 1;
fi;

# cleanup
rm -r xsim* xelab* webtalk* xvhdl* test.wdb 2> /dev/null || true;

echo " ## Compiling VHDL files: $VHDLS";
for V in $VHDLS; do
    xvhdl ${V/<IMPL>/$IMPL} || exit 2;
    grep -q ERROR xvhdl.log && exit 2;
done;

echo " ## Elaborating: ";
xelab testbench -s test -debug all || exit 3;
grep -q ERROR xelab.log && exit 3;

if [[ "$1" == "--gui" ]]; then
    echo " ## Running simulation in the GUI: ";
    xsim test --gui
else
    echo " ## Running simulation in batch mode: ";
    xsim test -R || exit 4;
    grep -q ERROR xsim.log && exit 4;

    test -f output_vhdl_tb.txt && echo " ## Output produced in output_vhdl_tb.txt ";
fi
