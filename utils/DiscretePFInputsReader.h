#ifndef utils_DiscretePFInputsReader_h
#define utils_DiscretePFInputsReader_h

#include <vector>
#include <cassert>
#include "../firmware/data.h"
#include "../DiscretePFInputs.h"
#include "DiscretePFInputsIO.h"
#include "DiscretePF2Firmware.h"
#include <iostream>

typedef l1tpf_impl::InputRegion Region;

struct Event {
    uint32_t run, lumi; uint64_t event;
    float z0, genZ0;
    std::vector<float> puGlobals;
    std::vector<Region> regions;
    
    Event() : run(0), lumi(0), event(0), z0(0.), regions(), puGlobals() {}
    bool readFromFile(FILE *fRegionDump) {
        if (!fread(&run, sizeof(uint32_t), 1, fRegionDump)) {
            std::cout << "Cannot read from file (eof?)" << std::endl;
            return false;
        }
        fread(&lumi, sizeof(uint32_t), 1, fRegionDump);
        fread(&event, sizeof(uint64_t), 1, fRegionDump);

        //std::cout << "Reading run " << run << " lumi " << lumi <<  " event " << event << std::endl;

        l1tpf_impl::readManyFromFile(regions, fRegionDump); 

        fread(&z0, sizeof(float), 1, fRegionDump);
        fread(&genZ0, sizeof(float), 1, fRegionDump);

        l1tpf_impl::readManyFromFile(puGlobals, fRegionDump); 

        return true;
    }
};

class DiscretePFInputsReader {
    public:
        DiscretePFInputsReader(const char *fileName) : file_(fopen(fileName,"rb")), iregion_(0) {
            if (!file_) { std::cout << "ERROR: cannot read '" << fileName << "'" << std::endl; }
            assert(file_);
            std::cout << "INFO: opening "  << fileName << std::endl;
        }
        ~DiscretePFInputsReader() { fclose(file_); }
        // for region-by-region approach
        bool nextRegion(HadCaloObj calo[NCALO], EmCaloObj emcalo[NEMCALO], TkObj track[NTRACK], MuObj mu[NMU], z0_t & hwZPV) {
            if (!nextRegion()) return false;
            const Region &r = event_.regions[iregion_];

            dpf2fw::convert<NTRACK>(r.track, track);
            dpf2fw::convert<NCALO>(r.calo, calo);
            dpf2fw::convert<NEMCALO>(r.emcalo, emcalo);
            dpf2fw::convert<NMU>(r.muon, mu);

#ifdef __GXX_EXPERIMENTAL_CXX0X__
            hwZPV = event_.z0 * l1tpf_impl::InputTrack::Z0_SCALE;
#else
            hwZPV = event_.z0 * 20;
#endif

            //printf("Read region %u with %lu tracks, %lu em calo, %lu had calo, %lu muons\n", iregion_, r.track.size(), r.emcalo.size(), r.calo.size(), r.muon.size());
            iregion_++;
            return true;
        }
        // for full event approach (don't mix with the above)
        bool nextEvent() {
            if (feof(file_)) return false;
            if (!event_.readFromFile(file_)) return false;
            return true;
        }
        const Event & event() { return event_; }

    private:
        bool nextRegion() {
            while(true) {
                if (event_.event == 0 || iregion_ == event_.regions.size()) {
                    if (feof(file_)) { std::cout << "END OF FILE" << std::endl; return false; }
                    if (!event_.readFromFile(file_)) return false;
                    //printf("Beginning of run %u, lumi %u, event %lu \n", event_.run, event_.lumi, event_.event);
                    iregion_ = 0;
                }
               return true;
            }
        }

        FILE *file_;
        Event event_;
        unsigned int iregion_;
};
#endif
