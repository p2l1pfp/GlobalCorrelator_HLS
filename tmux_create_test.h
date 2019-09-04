#include <cstdio>
#include <iomanip>
#include "firmware/simple_fullpfalgo.h"
#include "vertexing/firmware/simple_vtx.h"
#include "puppi/firmware/simple_puppi.h"
#include "utils/random_inputs.h"
#include "utils/DiscretePFInputs_IO.h"
#include "utils/pattern_serializer.h"
#include "utils/test_utils.h"
#include "firmware/mp7pf_encoding.h"

#define NTEST 6
#define NLINKS_APX_GEN0 96
#define NFRAMES_APX_GEN0 3

#define NLINKS_PER_TRACK 10
#define NLINKS_PER_CALO 10
#define NLINKS_PER_EMCALO 10
#define NLINKS_PER_MU 2
#define NLINKS_PER_REG (NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO+NLINKS_PER_MU)

#define MAXETA_INT 243
#define MAXPHI_INT 512
#define NPHI_INT 1024

#define ETA_BUFFER 32
#define PHI_BUFFER 32

void mp7wrapped_pack_in_reorder(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], MP7DataWord data[MP7_NCHANN]) {
    // pack inputs
    assert(2*NEMCALO + 2*NTRACK + 2*NCALO + 2*NMU <= MP7_NCHANN);
    #define EMOFFS 2*NTRACK
    #define HADOFFS 2*NEMCALO+EMOFFS
    #define MUOFFS 2*NCALO+HADOFFS
    mp7_pack<NTRACK,0>(track, data);
    mp7_pack<NEMCALO,EMOFFS>(emcalo, data);
    mp7_pack<NCALO,HADOFFS>(hadcalo, data);
    mp7_pack<NMU,MUOFFS>(mu, data);
}

template<typename T, int NIn, int NOut>
void ptsort_out(T in[NIn], T out[NOut]) {
    T tmp[NOut];

    for (int iout = 0; iout < NOut; ++iout) {
        tmp[iout].hwPt = 0;
    }

    for (int it = 0; it < NIn; ++it) {
        for (int iout = NOut-1; iout >= 0; --iout) {
            if (tmp[iout].hwPt <= in[it].hwPt) {
                if (iout == 0 || tmp[iout-1].hwPt > in[it].hwPt) {
                    tmp[iout] = in[it];
                } else {
                    tmp[iout] = tmp[iout-1];
                }
            }
        }

    }
    for (int iout = 0; iout < NOut; ++iout) {
        out[iout] = tmp[iout];
    }

}

void mp7wrapped_pfalgo3_full_sort(MP7DataWord input[MP7_NCHANN], MP7DataWord output[MP7_NCHANN], z0_t Z0) {
    
    PFChargedObj pfch[NTRACK]; PFNeutralObj pfpho[NPHOTON]; PFNeutralObj pfne[NSELCALO]; PFChargedObj pfmu[NMU];
    PFChargedObj pfch_sort[NTRACK]; PFNeutralObj pfpho_sort[NPHOTON]; PFNeutralObj pfne_sort[NSELCALO]; PFChargedObj pfmu_sort[NMU];

    mp7wrapped_pfalgo3_full(input, output, Z0);
    mp7wrapped_unpack_out_necomb(output, pfch, pfpho, pfne, pfmu);
    
    ptsort_out<PFChargedObj,NTRACK,NTRACK>(pfch, pfch_sort);
    ptsort_out<PFNeutralObj,NPHOTON,NPHOTON>(pfpho, pfpho_sort);
    ptsort_out<PFNeutralObj,NSELCALO,NSELCALO>(pfne, pfne_sort);
    ptsort_out<PFChargedObj,NMU,NMU>(pfmu, pfmu_sort);

    mp7wrapped_pack_out(pfch_sort, pfpho_sort, pfne_sort, pfmu_sort, output);

}

// LEQ LE
bool isInPhiRegion(int test, int loBound, int hiBound){
    // check whether "lowBound <= test < hiBound" while accounting for
    // "wraparound" when hi/loBound are outside of the (MAXPHI_INT-NPHI_INT,MAXPHI_INT) range    
    if(loBound<=MAXPHI_INT-NPHI_INT && hiBound>MAXPHI_INT) return true;
    if(loBound<MAXPHI_INT-NPHI_INT) return (test < hiBound) || (test >= MAXPHI_INT-((MAXPHI_INT-NPHI_INT)-loBound));
    if(hiBound>MAXPHI_INT) return (test >= loBound) || (test < (MAXPHI_INT-NPHI_INT)+(hiBound-MAXPHI_INT));
    return (test < hiBound) && (test >= loBound);
}
// LE LEQ
/* bool isInPhiRegion(int test, int loBound, int hiBound){ */
/*     // check whether "lowBound < test <= hiBound" while accounting for */
/*     // "wraparound" when hi/loBound are outside of the (MAXPHI_INT-NPHI_INT,MAXPHI_INT) range     */
/*     if(loBound<=MAXPHI_INT-NPHI_INT && hiBound>MAXPHI_INT) return true; */
/*     if(loBound<MAXPHI_INT-NPHI_INT) return (test<= hiBound) || (test>MAXPHI_INT-((MAXPHI_INT-NPHI_INT)-loBound)); */
/*     if(hiBound>MAXPHI_INT) return (test>loBound) || (test<= (MAXPHI_INT-NPHI_INT)+(hiBound-MAXPHI_INT)); */
/*     return (test<=hiBound) && (test>loBound); */
/* } */
