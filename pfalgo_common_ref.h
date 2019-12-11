#ifndef PFALGO_COMMON_REF_H
#define PFALGO_COMMON_REF_H

#include "firmware/data.h"
#include "firmware/pfalgo_common.h"

template <typename T> inline int sqr(const T & t) { return t*t; }

template<int NCAL, int DR2MAX, typename CO_t>
int best_match_with_pt_ref(const CO_t calo[NCAL], const TkObj & track) ;
 
template<typename T, int NIn, int NOut>
void ptsort_ref(const T in[NIn], T out[NOut]) ;

template<unsigned int NTrack, unsigned int NMu>
void pfalgo_mu_ref(const TkObj track[NTrack], const MuObj mu[NMu], bool isMu[NTrack], PFChargedObj outmu[NMu], bool debug) ;

#endif
