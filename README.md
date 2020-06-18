# GlobalCorrelator_HLS
Vivado High Level Synthesis framework for Global Correlator

## Structure of this repository:
 * `firmware` directory: core firmware for the L1PF, and associated dataformat header files
 * `ref` directory: reference L1PF implementation for bitwise validation in the HLS testbench
 * `run_hls_*.tcl` files: files to run the c++ simulation and synthesis of L1PF
 * `*_test.cpp` files: c++ testbench files for L1PF
 * `data` directory: input data files used by the testbench
 * `utils` directory: utilities used by the testbench (e.g. to read input data files or do validation)

Additional subdirectories contain other projects. In particular, `puppi` contains the L1 linearized puppi implementation, with a structure similar to the current one

## Setup of the code

The compilation of the firmware depends on some compile-time constants defining the detector region, input object multiplicities, parameters of the algorithm, and test setup used.

The detector region is set with `-DREG_(region)` with region being one of `Barrel`, `HGCal`, `HGCalNoTK`, `HF`. This drives the configuration of the object multiplicities (set in `firmware/data.h` and potentially further customized depending on the test board used), and the parameters of the algorithms (set in the corresponding header files; currently this is the case only for puppi, not for PF).

The board used is set with `-DBOARD_(board)`, which drives the code used to serialize the inputs and outputs (word size and number of channels), and possibly customizes the multiplicies, all in `firmware/data.h`. If no board is defined, or `-DBOARD_none`, the testbench will just run the algorithm without input and output data packing.

## Implementation status for L1PF & L1 Linearized Puppi

* PFAlgo3 (i.e. using emcalo, hadcalo, tracks, muons) is implemented and tested but only in the Barrel region. &Delta;R cuts are synchronized with CMSSW, but other parameters have not been checked recently.
* PFAlgo2HGC (i.e. using calo, tracks, muons) is implemented and tested but only in the HGCal region. &Delta;R cuts are synchronized with CMSSW, but other parameters have not been checked recently.
* Linearized puppi inside the tracker coverage is implemented and tested both in the Barrel and in the HGCal (where it has parameters split in 2 eta bins). Parameters are also synchronized with CMSSW
* In the forward region outside tracker coverage, linearized puppi goes in one step from calo clusters to puppi candidates, and is implemented and tested for both HF and HGCalNoTK. Parameters are also synchronized with CMSSW.

## Pending items

In random order:
* Investigate whether we can define integer coordinates in a clever way to have 2&pi; be equal to a power of 2, so that the wrapping of &Delta;&phi; would happen automatically, and we could use global coordinates.
* Separate in the inputs the tracks regionized using calorimeter (eta, phi) - for PF - and the tracks regionized using vertex (eta, phi) - for Puppi.
* Puppi implementation can probably be cleaned up using `ap_fixed` instead of bitshifts by hand
* LinearPuppi firmware implementation with eta bins can probably be cleaned up
* Introduce options to customize the pattern file layout depending on the board (now it can be done only by changing the options passed in the constructor in the testbench source)
* Cleanup & resurrect the remaining elements of the TDR demonstrator (vertexing, regionizer, layer 2, ...)
