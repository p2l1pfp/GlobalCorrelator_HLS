#include "firmware/regionizer.h"

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cassert>
#include <vector>

bool readEventTk(FILE *file, std::vector<TkObj> inputs[NTKSECTORS][NTKFIBERS], uint32_t &irun, uint32_t &ilumi, uint64_t &ievent) {
    if (feof(file)) return false;

    uint32_t run, lumi; uint64_t event;
    if (fscanf(file, "event %u %u %lu\n", &run, &lumi, &event) != 3) return false;
    if (irun == 0 && ilumi == 0 && ievent == 0) { 
        irun = run; ilumi = lumi; ievent = event; 
    } else if (irun != run || ilumi != lumi || ievent != event) {
        printf("event number mismatch: read  %u %u %lu, expected  %u %u %lu\n", run, lumi, event, irun, ilumi, ievent);
        return false;
    }
    //printf("reading event  %u %u %lu\n", run, lumi, event);

    int nfound = 0, maxfib = 0, maxsec = 0;
    for (int s = 0; s < NTKSECTORS; ++s) {
        int sec; uint64_t ntracks;
        if (fscanf(file, "sector %d tracks %lu\n", &sec, &ntracks) != 2) return false;
        assert(sec == s);
        //printf("reading sector %d -> %d tracks\n", sec, int(ntracks));
        for (int f = 0; f < NTKFIBERS; ++f) inputs[s][f].clear();
        for (int i = 0, n = ntracks; i < n; ++i) {
            int hwPt, hwEta, hwPhi, hwCaloPtErr, hwZ0, hwCharge, hwTight;
            //printf("read track %d/%d of sector %d\n", i, n, sec);
            int ret = fscanf(file, "track ipt %d ieta %d iphi %d ipterr %d iz0 %d icharge %d iqual %d\n",
                                &hwPt, &hwEta, &hwPhi, &hwCaloPtErr, &hwZ0, &hwCharge, &hwTight);
            if (ret != 7) return false;
            TkObj t;
            t.hwPt = hwPt; t.hwEta = hwEta; t.hwPhi = hwPhi;
            t.hwPtErr = hwCaloPtErr; t.hwZ0 = hwZ0; t.hwCharge = hwCharge; t.hwTightQuality = hwTight;
            inputs[s][i % NTKFIBERS].push_back(t);
            nfound++;
            maxfib = std::max<int>(maxfib, inputs[s][i % NTKFIBERS].size());
        }
        maxsec = std::max<int>(maxsec, ntracks);
    }
    //printf("read %d tracks for this event. max %d tracks/sector, %d tracks/fiber\n", nfound, maxsec, maxfib);
    return true;
}


bool readEventCalo(FILE *file, std::vector<HadCaloObj> inputs[NCALOSECTORS][NCALOFIBERS], bool zside, uint32_t &irun, uint32_t &ilumi, uint64_t &ievent) {
    if (feof(file)) return false;

    uint32_t run, lumi; uint64_t event;
    if (fscanf(file, "event %u %u %lu\n", &run, &lumi, &event) != 3) return false;
    if (irun == 0 && ilumi == 0 && ievent == 0) { 
        irun = run; ilumi = lumi; ievent = event; 
    } else if (irun != run || ilumi != lumi || ievent != event) {
        printf("event number mismatch: read  %u %u %lu, expected  %u %u %lu\n", run, lumi, event, irun, ilumi, ievent);
        return false;
    }
    //printf("reading event  %u %u %lu\n", run, lumi, event);

    for (int s = 0; s < NCALOSECTORS; ++s) {
        for (int f = 0; f < NCALOFIBERS; ++f) inputs[s][f].clear();
    }

    int nfound = 0, maxfib = 0, maxsec = 0;
    for (int s = 0; s < 2*NCALOSECTORS; ++s) {
        int zs, sec; uint64_t nclusters;
        if (fscanf(file, "zside %d sector %d cluster %lu\n", &zs, &sec, &nclusters) != 3) return false;
        //printf("reading zside %d sector %d -> %d clusters\n", zs, sec, int(nclusters));
        for (int i = 0, n = nclusters; i < n; ++i) {
            int hwPt, hwEta, hwPhi, hwPtErr, hwEmPt, hwIsEM;
            int ret = fscanf(file, "cluster ipt %d ieta %d iphi %d iempt %d ipterr %d isem %1d\n",
                                &hwPt, &hwEta, &hwPhi, &hwEmPt, &hwPtErr, &hwIsEM);
            if (ret != 6) return false;
            if (zs == zside) {
                HadCaloObj t;
                t.hwPt = hwPt; t.hwEta = hwEta; t.hwPhi = hwPhi; 
                t.hwEmPt = hwEmPt; t.hwIsEM = hwIsEM;
                inputs[sec][i % NCALOFIBERS].push_back(t);
                nfound++;
                maxfib = std::max<int>(maxfib, inputs[sec][i % NCALOFIBERS].size());
            }
        }
        maxsec = std::max<int>(maxsec, nclusters);
    }
    //printf("read %d clusters for this event. max %d clusters/sector, %d clusters/fiber\n", nfound, maxsec, maxfib);
    return true;
}



bool readEventMu(FILE *file, std::vector<GlbMuObj> inputs[NMUFIBERS], uint32_t &irun, uint32_t &ilumi, uint64_t &ievent) {
    if (feof(file)) return false;

    uint32_t run, lumi; uint64_t event;
    if (fscanf(file, "event %u %u %lu\n", &run, &lumi, &event) != 3) return false;
    if (irun == 0 && ilumi == 0 && ievent == 0) { 
        irun = run; ilumi = lumi; ievent = event; 
    } else if (irun != run || ilumi != lumi || ievent != event) {
        printf("event number mismatch: read  %u %u %lu, expected  %u %u %lu\n", run, lumi, event, irun, ilumi, ievent);
        return false;
    }
    //printf("reading event  %u %u %lu\n", run, lumi, event); fflush(stdout);

    for (int f = 0; f < NMUFIBERS; ++f) inputs[f].clear();

    int nfound = 0;
    uint64_t nmuons;
    if (fscanf(file, "muons %lu\n", &nmuons) != 1) return false;
    //printf("reading -> %d muons\n", int(nmuons)); fflush(stdout);
    for (int i = 0, n = nmuons; i < n; ++i) {
        int hwPt, hwEta, hwPhi, hwCharge, hwQual;
        int ret = fscanf(file, "muon ipt %d ieta %d iphi %d icharge %d iqual %d\n",
                                &hwPt, &hwEta, &hwPhi, &hwCharge, &hwQual);
        if (ret != 5) return false;
        GlbMuObj t;
        t.hwPt = hwPt; t.hwEta = hwEta; t.hwPhi = hwPhi; 
        t.hwPtErr = 0;
        //t.hwCharge = hwCharge; t.hwQual = hwQual;
        inputs[i % NMUFIBERS].push_back(t);
        nfound++;
    }
    //printf("read %d muons for this event\n", nfound); fflush(stdout);
    return true;
}

