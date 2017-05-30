#ifndef _GLOBAL_CORRELATOR_H_
#define _GLOBAL_CORRELATOR_H_

#include <stdio.h>
#include "ap_int.h"

#define NUM_TRANS 40
#define ieta_max 10
#define iphi_max 10
#define ncalo_max 20
#define ntrack_max 20
#define npf_neutral_max 20

typedef ap_uint<10> int10;
typedef int din_t;
typedef int dint_t;
typedef int dout_t;

void hier_func(din_t A, din_t B, dout_t *C, dout_t *D);

void calo_track_linking_grid(int10 calos_pt[ieta_max][iphi_max],
							 int10 calos_h_over_e[ieta_max][iphi_max],
							 int10 tracks_pt[ieta_max][iphi_max],
							 int10 pf_neutral_pt[ieta_max][iphi_max]);

void calo_track_linking_list(int10 track_list_pt[ntrack_max],
							 int10 track_list_ieta[ntrack_max],
							 int10 track_list_iphi[ntrack_max],
							 int10 calo_list_pt[ncalo_max],
							 int10 calo_list_h_over_e[ncalo_max],
							 int10 calo_list_ieta[ncalo_max],
							 int10 calo_list_iphi[ncalo_max],
							 int10 pf_neutral_list_pt[npf_neutral_max],
							 int10 pf_neutral_list_ieta[npf_neutral_max],
							 int10 pf_neutral_list_iphi[npf_neutral_max]);


#endif
