#!/bin/sh

PF_DIR=/data/drankin/GlobalCorrelator_HLS/

echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t "

if [[ "$1x" == "x" ]]; then
    echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t usage: ./run_barebones.sh <NTRACKS> <NCALOS> <NEMCALOS> <II> <CLOCK_PERIOD> <CLOCK_FREQ> <PARTID>" 
    exit
fi

SCRIPT_DIR="$( 
  cd "$(dirname "$(readlink "$0" || printf %s "$0")")"
  pwd -P 
)"


echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t Script Path = $SCRIPT_DIR"

NTRACKS=$1
echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t NTRACKS = $NTRACKS"
NCALOS=$2
echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t NCALOS = $NCALOS"
NEMCALOS=$3
echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t NEMCALOS = $NEMCALOS"
II=$4
echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t II = $II"
CLOCK_PERIOD=$5
echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t CLOCK_PERIOD = $CLOCK_PERIOD"
CLOCK_FREQ=$6
echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t CLOCK_FREQ = $CLOCK_FREQ"
PARTID=$7
echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t PARTID = $PARTID"

# setup Vivado
. /home/drankin/setup_vivado.2018.1.sh #setup vivado and vivado_hls (2018.1)


###################### HLS

echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t Preparing for HLS..."

SRC_HLS_DIR=${PF_DIR}firmware/data.h.bak
DEST_HLS_DIR=${PF_DIR}firmware/data.h

echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t Source = $SRC_HLS_DIR"
echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t Destination = $DEST_HLS_DIR"

if true; then #for debugging

	cp $SRC_HLS_DIR $DEST_HLS_DIR

	sed -i "s/NTRACK [0-9]*$/NTRACK ${NTRACKS}/g" $DEST_HLS_DIR
	sed -i "s/NCALO [0-9]*$/NCALO ${NCALOS}/g" $DEST_HLS_DIR
	sed -i "s/NEMCALO [0-9]*$/NEMCALO ${NEMCALOS}/g" $DEST_HLS_DIR
	sed -i "s/NSELCALO [0-9]*$/NSELCALO `expr ${NCALOS} - 5`/g" $DEST_HLS_DIR
	sed -i "s/MP7_NCHANN [0-9]*$/MP7_NCHANN `expr ${NTRACKS} + ${NTRACKS} + ${NCALOS} + ${NCALOS} + ${NCALOS} + ${NCALOS} + ${NCALOS} + ${NCALOS} + ${NEMCALOS} + ${NEMCALOS} - 10 + 4`/g" $DEST_HLS_DIR
	
	echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t Running HLS..."
	echo
	
        sed "s/II=[0-9]*/II=${II}/g" run_hls_fullpfalgo_w_puppi.tcl > run_hls_test.tcl
        sed -i "s/create_clock -period [0-9\.]*/create_clock -period ${CLOCK_PERIOD}/g" run_hls_test.tcl
        sed -i "s/set_part {.*}/set_part {${PARTID}}/g" run_hls_test.tcl

	cd ${PF_DIR}
        #rm -rf l1pf-resource-test
	#vivado_hls -f run_hls_fullpfalgo_test.tcl
        rm -rf l1pfpuppi-resource-test
        #vivado_hls -f run_hls_fullpfalgo_w_puppi.tcl
        vivado_hls -f run_hls_test.tcl

	#cd ${PF_DIR}/puppi
        #rm -rf proj_pfpuppi
	#vivado_hls -f run_hls_pfpuppi.tcl

	cd ${PF_DIR}
	cp $SRC_HLS_DIR $DEST_HLS_DIR

fi

###################### Build vivado project
if true; then #for debugging

	echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t Preparing for Vivado project..."

	SRC_BUILD_DIR=$SCRIPT_DIR/build.tcl
	#SRC_BUILD_DIR=$SCRIPT_DIR/build_forward.tcl
	DEST_BUILD_DIR=$SCRIPT_DIR/modbuild.tcl

	echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t Source = $SRC_BUILD_DIR"
	echo -e `date +"%h%y %T"` "BB [${LINENO}]  \t Destination = $DEST_BUILD_DIR"

	cp $SRC_BUILD_DIR $DEST_BUILD_DIR

#------CENTRAL
	sed -i s/XX_NINPUTS_XX/$((2*(${NTRACKS} + (${NEMCALOS} + ${NCALOS})) + 4))/g $DEST_BUILD_DIR
	sed -i s/XX_NOUTPUTS_XX/$((2*(${NTRACKS} + (${NEMCALOS} + ${NCALOS})) - 6))/g $DEST_BUILD_DIR
	sed -i s/XX_NPFCH_XX/${NTRACKS}/g $DEST_BUILD_DIR
	sed -i s/XX_NNEUTRAL_XX/$(((${NEMCALOS} + ${NCALOS}) - 5))/g $DEST_BUILD_DIR

#-------FORWARD
#	sed -i s/XX_NINPUTS_XX/$((2*(${NTRACKS} + (2 * ${NCALOS})) + 4))/g $DEST_BUILD_DIR
#	sed -i s/XX_NOUTPUTS_XX/$((2*(${NTRACKS} + (2 * ${NCALOS})) + 4))/g $DEST_BUILD_DIR
#	sed -i s/XX_NPFCH_XX/${NTRACKS}/g $DEST_BUILD_DIR
#	sed -i s/XX_NNEUTRAL_XX/$((2 * ${NCALOS}))/g $DEST_BUILD_DIR

        sed -i s/XX_CLOCK_XX/${CLOCK_FREQ}/g $DEST_BUILD_DIR
        sed -i s/XX_PART_XX/${PARTID}/g $DEST_BUILD_DIR

	cd $SCRIPT_DIR

	rm -rf pftest_tieoff vivado*.log vivado*.jou
	vivado -mode batch -source $DEST_BUILD_DIR

fi

