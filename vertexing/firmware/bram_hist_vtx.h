#include "../../firmware/data.h"
#include "../../regionizer/firmware/regionizer.h"

typedef ap_uint<9> ptsum_t;
typedef ap_uint<18> twoptsums_t;
typedef ap_uint<10> zbin_t;
struct zbin_vt {
    zbin_t bin;
    bool   valid;
};
// the current z0 scale is 20 counts per CM (1 bit = 0.5 mm)
// we divide down by 8 to have 1 bit = 0.4 cm. 
// with 71 bins, we cover +/- 14 cm, i.e. about +/- 3 sigma(Z)

// the standard pt scale is 0.25 GeV per unit
// if we max the track pt to 50 GeV, it would be 200 units
// we divide down by two
#define BHV_MAXPT 100
//#define BHV_MAXBIN 511
#define BHV_MAXBIN 511

#define BNV_SHIFT 3
#define BHV_NBINS 72
#define BHV_NHALFBINS (BHV_NBINS/2)
#define BHV_NSECTORS 2*N_IN_SECTORS
#define BHV_NTRACKS 18


inline zbin_vt fetch_bin_ref(z0_t z0) {
    int zbin = (z0 >> BNV_SHIFT) + BHV_NHALFBINS; 
    bool valid = true;
    if (zbin < 0) { zbin = 0; valid = false; }
    if (zbin > BHV_NBINS-1) { zbin = 0; valid = false; }
    zbin_vt ret;
    ret.bin = zbin;
    ret.valid = valid;
    return ret;
}
inline z0_t bin_center_ref(zbin_t iz) {
    int z = int(iz) - BHV_NHALFBINS;
    return z0_t((z << BNV_SHIFT) + ( 1 << (BNV_SHIFT-1) ));
}
void bhv_add_track(zbin_vt zbin, pt_t tkpt, ptsum_t hist[BHV_NBINS]) ;
zbin_t bhv_find_pv(twoptsums_t hist[BHV_NSECTORS][BHV_NHALFBINS], pt_t *sumpt) ;
void bhv_find_pv_ref(TkObj tracks[BHV_NSECTORS][BHV_NTRACKS], zbin_t & pvbin, z0_t & pv, int &pvsum) ;
bool dummy(z0_t z0) ;
