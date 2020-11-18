#include "firmware/tdemux.h"
#include "tdemux_ref.h"
#include <algorithm>

TDemuxRef::TDemuxRef() {
    counter = 0;
    robin = 0;
    for (int i = 0; i < NLINKS; ++i) {
        fold[i] = (i == 0 ? 1 : 0);
        offs[i]  =  i * BLKSIZE + (fold[i] ? PAGESIZE : 0);
        std::fill_n(&buffer[i][0], MEMSIZE, 0);
    }
    toread = -1;
    readcount = 0;
}

bool TDemuxRef::operator()(bool newEvent, const w65 links[NLINKS], w65 out[NLINKS]) {
    if (newEvent) {
        counter  = 0;
        robin = 0;
        for (int i = 0; i < NLINKS; ++i) {
            fold[i] = (i == 0 ? 1 : 0);
            offs[i]  =  i * BLKSIZE + (fold[i] ? PAGESIZE : 0);
        }
        toread = -1;
        //printf("Initialized\n");
    }
    for (int i = 0; i < NLINKS; ++i) {
        buffer[(robin+i)%NLINKS][offs[i]] = links[i];
    }
    robin++;
    counter++;

    /*printf("counter %4d, robin %d, toread %6d, readrobin %d, readcount %3d\n", counter, robin, toread, readrobin, readcount);
    for (unsigned int j = 0; j < NLINKS; ++j) {
        printf("M%d  fold %d  offs %5d : ", j, fold[j] ? 1 : 0, offs[j]);
        for (unsigned int i = 0; i < MEMSIZE; ++i) {
            printf("%5d%c | ", int(buffer[j][i]), i == offs[(NLINKS-robin+1+j)%NLINKS] ? '!' : ' ');
        }
        printf("\n");
    }*/

    if (robin == NLINKS) {
        robin = 0;
        if ((counter % BLKSIZE) == 0) {
            int completed = (counter/BLKSIZE)%NLINKS; // index of the link whose buffer we have filled
            //printf("Swap block for %d as counter %d is multiple of %d, %d\n", completed, counter, BLKSIZE, (counter % BLKSIZE));
            fold[completed] = !fold[completed];  // fold that one around
            for (int i = 0; i < NLINKS; ++i) {
                offs[i]  =  (i == completed ? i * BLKSIZE + (fold[i] ? PAGESIZE : 0) : offs[i]+1);
            }
        } else  {
            //printf("Increase offset as counter %d is not multiple of %d\n", counter, BLKSIZE);
            for (int i = 0; i < NLINKS; ++i) {
                offs[i]++;
            }
        }
    }
    if (toread == -1) {
        //printf("counter %4d, threshold %4d, toread %6d, readrobin %d, readcount %3d\n", counter, (NLINKS-1)*BLKSIZE, toread, readrobin, readcount);
        if (counter == (NLINKS-1)*BLKSIZE + 1) { 
            toread = PAGESIZE;
            readrobin = 0; readcount = 0;
        }
        for (int i = 0; i < NLINKS; ++i) {
            out[i] = 0;
        }
        return false;
    } else {
        for (int i = 0; i < NLINKS; ++i) {
            out[i] = buffer[(i+readrobin)%NLINKS][toread];
        }
        toread    = (toread+1)    % MEMSIZE;
        readcount = (readcount+1) % BLKSIZE;
        if (readcount == 0) readrobin = (readrobin+1) % NLINKS;
        return true;
    }
}

bool TDemuxRef::operator()(bool newEvent, const w64 links[NLINKS], const bool valid[NLINKS], 
                                                w64 out[NLINKS],         bool vout[NLINKS]) {
    w65 links65[NLINKS], out65[NLINKS];
    for (int i = 0; i < NLINKS; ++i) {
        links65[i](63,0) = links[i];
        links65[i][64]   = valid[i];
    }
    (*this)(newEvent, links65, out65);
    for (int i = 0; i < NLINKS; ++i) {
        out[i]  = out65[i](63,0);
        vout[i] = out65[i][64];
    }
}
