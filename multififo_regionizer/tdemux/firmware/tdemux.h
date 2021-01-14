#ifndef multififo_regionizer_tdemux_h
#define multififo_regionizer_tdemux_h

#include <ap_int.h>

typedef ap_uint<65> w65; // bit 64 is used for the valid bit

#define TMUX_IN 18
#define TMUX_OUT 6
#define NLINKS   3 
#define NCLK     9 // clocks per BX (8 = 320 MHz, 9 = 360 MHz)
#define BLKSIZE  (NCLK*TMUX_OUT)
#define PAGESIZE (NCLK*TMUX_IN)

bool tdemux(bool newEvent, const w65 links[NLINKS], w65 out[NLINKS]) ;

#endif
