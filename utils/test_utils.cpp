#include "test_utils.h"
#include <cstdio>

bool had_equals(const HadCaloObj &out_ref, const HadCaloObj &out, const char *what, int idx) {
    bool ret;
    if (out_ref.hwPt == 0) {
        ret = (out.hwPt == 0);
    } else {
        ret = (out_ref.hwPt == out.hwPt && out_ref.hwEmPt == out.hwEmPt && out_ref.hwEta == out.hwEta && out_ref.hwPhi == out.hwPhi && out_ref.hwIsEM  == out.hwIsEM);
    }
    if  (!ret) {
        printf("Mismatch at %s[%d] ref vs test, hwPt % 7d % 7d   hwEmPt % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwId %1d %1d\n", what, idx,
                int(out_ref.hwPt), int(out.hwPt), int(out_ref.hwEmPt), int(out.hwEmPt), int(out_ref.hwEta), int(out.hwEta), int(out_ref.hwPhi), int(out.hwPhi), int(out_ref.hwIsEM), int(out.hwIsEM));
    }
    return ret;
}
bool em_equals(const EmCaloObj &out_ref, const EmCaloObj &out, const char *what, int idx) {
    bool ret;
    if (out_ref.hwPt == 0) {
        ret = (out.hwPt == 0);
    } else {
        ret = (out_ref.hwPt == out.hwPt && out_ref.hwPtErr == out.hwPtErr && out_ref.hwEta == out.hwEta && out_ref.hwPhi == out.hwPhi);
    }
    if  (!ret) {
        printf("Mismatch at %s[%d] ref vs test, hwPt % 7d % 7d   hwPtErr % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d\n", what, idx,
                int(out_ref.hwPt), int(out.hwPt), int(out_ref.hwPtErr), int(out.hwPtErr), int(out_ref.hwEta), int(out.hwEta), int(out_ref.hwPhi), int(out.hwPhi));
    }
    return ret;
}

bool track_equals(const TkObj &out_ref, const TkObj &out, const char *what, int idx) {
    bool ret;
    if (out_ref.hwPt == 0) {
        ret = (out.hwPt == 0);
    } else {
        ret = (out_ref.hwPt == out.hwPt && out_ref.hwPtErr == out.hwPtErr && out_ref.hwEta == out.hwEta && out_ref.hwPhi == out.hwPhi && out_ref.hwZ0  == out.hwZ0);
    }
    if  (!ret) {
        printf("Mismatch at %s[%d] ref vs test, hwPt % 7d % 7d   hwPtErr % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwZ0 %+7d %+7d\n", what, idx,
                int(out_ref.hwPt), int(out.hwPt), int(out_ref.hwPtErr), int(out.hwPtErr), int(out_ref.hwEta), int(out.hwEta), int(out_ref.hwPhi), int(out.hwPhi), int(out_ref.hwZ0), int(out.hwZ0));
    }
    return ret;
}
bool mu_equals(const MuObj &out_ref, const MuObj &out, const char *what, int idx) {
    bool ret;
    if (out_ref.hwPt == 0) {
        ret = (out.hwPt == 0);
    } else {
        ret = (out_ref.hwPt == out.hwPt && out_ref.hwPtErr == out.hwPtErr && out_ref.hwEta == out.hwEta && out_ref.hwPhi == out.hwPhi);
    }
    if  (!ret) {
        printf("Mismatch at %s[%d] ref vs test, hwPt % 7d % 7d   hwPtErr % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d\n", what, idx,
                int(out_ref.hwPt), int(out.hwPt), int(out_ref.hwPtErr), int(out.hwPtErr), int(out_ref.hwEta), int(out.hwEta), int(out_ref.hwPhi), int(out.hwPhi));
    }
    return ret;
}

bool pf_equals(const PFChargedObj &out_ref, const PFChargedObj &out, const char *what, int idx) {
    bool ret;
    if (out_ref.hwPt == 0) {
        ret = (out.hwPt == 0);
    } else {
        ret = (out_ref.hwPt == out.hwPt && out_ref.hwEta == out.hwEta && out_ref.hwPhi == out.hwPhi && out_ref.hwId  == out.hwId && out_ref.hwZ0  == out.hwZ0);
    }
    if  (!ret) {
        printf("Mismatch at %s[%d] ref vs test, hwPt % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwId %1d %1d      hwZ0 %+7d %+7d   \n", what, idx,
                int(out_ref.hwPt), int(out.hwPt),
                int(out_ref.hwEta), int(out.hwEta),
                int(out_ref.hwPhi), int(out.hwPhi),
                int(out_ref.hwId), int(out.hwId),
                int(out_ref.hwZ0), int(out.hwZ0));
    }
    return ret;
}
bool pf_equals(const PFNeutralObj &out_ref, const PFNeutralObj &out, const char *what, int idx) {
    bool ret;
    if (out_ref.hwPt == 0) {
        ret = (out.hwPt == 0);
    } else {
        ret = (out_ref.hwPt == out.hwPt && out_ref.hwEta == out.hwEta && out_ref.hwPhi == out.hwPhi && out_ref.hwId  == out.hwId);
    }
    if  (!ret) {
        printf("Mismatch at %s[%d] ref vs test, hwPt % 7d % 7d   hwEta %+7d %+7d   hwPhi %+7d %+7d   hwId %1d %1d \n", what, idx,
                int(out_ref.hwPt), int(out.hwPt),
                int(out_ref.hwEta), int(out.hwEta),
                int(out_ref.hwPhi), int(out.hwPhi),
                int(out_ref.hwId), int(out.hwId));
    }
    return ret;
}

