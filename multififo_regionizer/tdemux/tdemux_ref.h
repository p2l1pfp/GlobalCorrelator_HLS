#ifndef multififo_regionizer_tdemux_ref_h
#define multififo_regionizer_tdemux_ref_h

#include "firmware/tdemux.h"

typedef ap_uint<64> w64;

class TDemuxRef {
    public:
        TDemuxRef() ;
        bool operator()(bool newEvent, const w65 links[NLINKS], w65 out[NLINKS]) ;
        bool operator()(bool newEvent, const w64 links[NLINKS], const bool valid[NLINKS], 
                                             w64 out[NLINKS],         bool vout[NLINKS]) ;
    private:
        static const unsigned int MEMSIZE = 2*PAGESIZE;
        int   counter; // must count up to TMUX_IN * NLCK
        bool  fold[NLINKS];    // we use twice as much memory, to avoid contentions
        int   offs[NLINKS];
        int   robin;
        int   toread, readrobin, readcount;

        w65 buffer[NLINKS][MEMSIZE];
};

#endif
