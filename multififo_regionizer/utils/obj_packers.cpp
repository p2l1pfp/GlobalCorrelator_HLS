#include "obj_packers.h"

#ifdef TRIVIAL_ENCODING_64

template<typename T>
std::vector<ap_uint<64>> pack_objs(const std::vector<T> & objs) {
    std::vector<ap_uint<64>> ret;
    for (unsigned int i = 0, n = objs.size(); i < n; ++i) {
        ret.push_back(l1pf_pattern_pack_one(objs[i]));
    }
    return ret;
}

std::vector<ap_uint<64>> pack_tracks(const std::vector<TkObj> & tracks) { return pack_objs(mu); } 
std::vector<ap_uint<64>> pack_hgcal(const std::vector<HadCaloObj> & calo) return pack_objs(mu); }
std::vector<ap_uint<64>> pack_muons(const std::vector<GlbMuObj> & mu) { return pack_objs(mu);  }

#else

std::vector<ap_uint<64>> pack_tracks(const std::vector<TkObj> & tracks) {
    std::vector<ap_uint<64>> ret;
    for (unsigned int i = 0, n = tracks.size(); i < n; ++i) {
        // simulate 96 bit objects
        ap_uint<96> packedtk = l1pf_pattern_pack_one(tracks[i]);
        if (i % 2 == 0) {
            ret.emplace_back(packedtk(95,32));
            ret.emplace_back((packedtk(31,0), ap_uint<32>(0)));
        } else {
            ret.back()(31,0) = packedtk(95,64);
            ret.emplace_back(packedtk(63,0));
        }
    }
    return ret;
}

std::vector<ap_uint<64>> pack_hgcal(const std::vector<HadCaloObj> & calo) {
    std::vector<ap_uint<64>> ret;
    for (unsigned int i = 0, n = calo.size(); i < n; ++i) {
        // simulate 128 bit objects on a 16 G link
        ap_uint<128> packed = l1pf_pattern_pack_one(calo[i]);
        ret.emplace_back(packed(127,64));
        ret.emplace_back(packed( 63, 0));
        ret.push_back(0); // third zero frame, that will have the strobe bit off
    }
    return ret;

}

std::vector<ap_uint<64>> pack_muons(const std::vector<GlbMuObj> & mu) {
    std::vector<ap_uint<64>> ret;
    for (unsigned int i = 0, n = mu.size(); i < n; ++i) {
        // simulate 128 bit objects on a 16 G link
        ap_uint<128> packed = l1pf_pattern_pack_one(mu[i]);
        ret.emplace_back(packed(127,64));
        ret.emplace_back(packed( 63, 0));
    }
    return ret;
}
#endif


