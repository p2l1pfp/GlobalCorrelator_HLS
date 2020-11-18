#include "firmware/tdemux.h"
#include <cstdio>
#include <cstdlib>
#include "tdemux_ref.h"


int main() {
    srand(125);
    const int NDATA = TMUX_IN*NCLK*4+10;
    w65 data[NLINKS][NDATA], in[NLINKS], out[NLINKS], refout[NLINKS];

    FILE * f_patterns_in, * f_patterns_out; char fnbuff[25]; 
    TDemuxRef tdemux_ref;

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
