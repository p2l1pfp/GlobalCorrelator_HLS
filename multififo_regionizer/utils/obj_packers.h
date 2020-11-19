#ifndef multififo_regionizer_obj_packers_h
#define multififo_regionizer_obj_packers_h

#include "../firmware/regionizer.h"
#include <vector>

std::vector<ap_uint<64>> pack_tracks(const std::vector<TkObj> & tracks) ; 
std::vector<ap_uint<64>> pack_hgcal(const std::vector<HadCaloObj> & calo) ;
std::vector<ap_uint<64>> pack_muons(const std::vector<GlbMuObj> & mu) ; 

#endif
