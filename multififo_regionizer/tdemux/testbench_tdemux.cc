#include "firmware/tdemux.h"
#include <cstdio>
#include <cstdlib>

bool tdemux_ref(bool newEvent, const w65 links[NLINKS], w65 out[NLINKS]) ;
bool tdemux_ref(bool newEvent, const w65 links[NLINKS], w65 out[NLINKS]) {
    static int   counter; // must count up to TMUX_IN * NLCK
    static bool  fold[NLINKS];    // we use twice as much memory, to avoid contentions
    static int   offs[NLINKS];
    static int   robin;
    static int   toread = -1, readrobin = 0, readcount = 0;

    const unsigned int MEMSIZE = 2*PAGESIZE;
    static w65 buffer[NLINKS][MEMSIZE];

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

int main() {
    srand(125);
    const int NDATA = TMUX_IN*NCLK*4+10;
    w65 data[NLINKS][NDATA], in[NLINKS], out[NLINKS], refout[NLINKS];

    FILE * f_patterns_in, * f_patterns_out; char fnbuff[25]; 

    for (unsigned int itest = 0, ntest = 20; itest <= ntest; ++itest) {
        // create some input data
        bool isok = true;
        for (unsigned int j = 0; j < NLINKS; ++j) {
            for (unsigned int i = 0; i < NDATA; ++i) {
                if (itest == 0){ // special case, human readable pattern
                    int iclock = i - j * TMUX_OUT * NCLK;
                    if (iclock >= 0) {
                        if (NCLK > 1) {
                            int sub = iclock % NCLK;
                            int bx  = (iclock / NCLK) % TMUX_IN;
                            int ev  = (iclock / (TMUX_IN * NCLK)) * NLINKS + j;
                            data[j][i](63,0) = 1 + sub + 10*bx + 1000 * ev;
                            data[j][i][64] = 1;
                        } else {
                            int bx  = iclock % TMUX_IN;
                            int ev  = iclock / TMUX_IN * NLINKS + j;
                            data[j][i](63,0) = 1  + bx + 100 * ev;
                            data[j][i][64] = 1;
                        }
                    } else {
                        data[j][i] = 0;
                    }
                } else {
                    data[j][i][64] = rand() & 0x1;
                    data[j][i](63,0) = ap_uint<65>(rand() & 0xFFFFFF);
                }
            }
        }
        if (itest == 0) {
            for (unsigned int j = 0; j < NLINKS; ++j) {
                printf("L[%d]: ", j);
                for (unsigned int i = 0; i < NDATA; ++i) printf("%5d | ", int(data[j][i]));
                printf("\n");
            }
        }

        snprintf(fnbuff, 25, "patterns-in-%d.txt", itest);
        f_patterns_in = fopen(fnbuff, "w");
        snprintf(fnbuff, 25, "patterns-out-%d.txt", itest);
        f_patterns_out = fopen(fnbuff, "w");


        for (unsigned int iclock = 0; iclock <  NDATA; ++iclock) {

            fprintf(f_patterns_in,  "Frame %04u :", iclock);
            fprintf(f_patterns_out, "Frame %04u :", iclock);

            for (unsigned int j = 0; j < NLINKS; ++j) {
                in[j] = data[j][iclock];
                fprintf(f_patterns_in, " %1dv%016llx", int(in[j][64]), in[j](63,0).to_uint64());
            }

            bool ret = tdemux(iclock == 0, in, out);
            bool ref = tdemux_ref(iclock == 0, in, refout);

            for (unsigned int j = 0; j < NLINKS; ++j) {
                fprintf(f_patterns_out, " %1dv%016llx", int(in[j][64]), refout[j](63,0).to_uint64());
            }

            bool ok = (ret == ref);
            if (ok) {
                for (unsigned int j = 0; j < NLINKS; ++j) {
                    ok = ok && (out[j] == refout[j]);
                }
            }

            if (itest == 0) {
                printf("%04d |  ", iclock);
                for (unsigned int j = 0; j < NLINKS; ++j) printf("%6d ", int(in[j]));
                printf(" | v%d  ", ret ? 1 : 0);
                for (unsigned int j = 0; j < NLINKS; ++j) printf("%6d ", int(out[j]));
                printf(" | v%d  ", ref ? 1 : 0);
                for (unsigned int j = 0; j < NLINKS; ++j) printf("%6d ", int(refout[j]));
                printf(isok ? "\n" : "   <=== ERROR \n");
            }

            if (!ok) isok = false;

            fprintf(f_patterns_in, "\n");
            fprintf(f_patterns_out, "\n");
        }
        if (!isok) {
            printf("\ntest %d failed\n", itest);
            return 1;
        } else {
            printf("\ntest %d passed\n", itest);
        }
        fclose(f_patterns_in);
        fclose(f_patterns_out);
    }
    return 0;
}
