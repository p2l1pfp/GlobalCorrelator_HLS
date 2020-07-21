//

#include "tk_input_converter.h"


void propagate_tanlam(tkz0_t z0, tanlam_t tanlam, tanlam_t &tanlam_at_det);

inline float tanlam_to_eta(float tanlam){return -log(tan((M_PI/2 - atan(tanlam))/2));}

