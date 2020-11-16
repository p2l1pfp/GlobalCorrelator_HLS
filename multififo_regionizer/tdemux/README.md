# `tdemux`: a time-demultiplexer from TM 18 to TM 6

A module that reads from 3 links corresponding to 3 different time slices at TM 18 with time offsets 0 BX, 6 BX, 12 BX, and reassembles the frames outputing events at TM 6.
Internally, it uses 6 BRAM36, to have the bandwith to store 3 x 64 bits input data (it could also be implemented with 3 URAMs, though I didn't try).
