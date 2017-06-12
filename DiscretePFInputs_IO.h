#ifndef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_IO_H
#define FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_IO_H

#include "DiscretePFInputs.h"
#include <vector>
#include <cassert>
#include "data.h"

struct Region {
	float etaCenter, etaMin, etaMax, phiCenter, phiHalfWidth, etaExtra, phiExtra;
	
	std::vector<l1tpf_int::CaloCluster> calo, emcalo;
	std::vector<l1tpf_int::PropagatedTrack> track;
	std::vector<l1tpf_int::Muon> muon;
	
	void readFromFile(FILE *file) {
	    fread(&etaCenter, sizeof(float), 1, file);
	    fread(&etaMin,    sizeof(float), 1, file);
	    fread(&etaMax,    sizeof(float), 1, file);
	    fread(&phiCenter, sizeof(float), 1, file);
	    fread(&phiHalfWidth, sizeof(float), 1, file);
	    fread(&etaExtra, sizeof(float), 1, file);
	    fread(&phiExtra, sizeof(float), 1, file);

	    uint32_t number; uint32_t size;
	    fread(&number, sizeof(uint32_t), 1, file);
	    fread(&size,   sizeof(uint32_t), 1, file);
	    calo.resize(number); assert(size == sizeof(l1tpf_int::CaloCluster));
	    fread(&calo[0], size, number, file);

	    fread(&number, sizeof(uint32_t), 1, file);
	    fread(&size,   sizeof(uint32_t), 1, file);
	    emcalo.resize(number); assert(size == sizeof(l1tpf_int::CaloCluster));
	    fread(&emcalo[0], size, number, file);

	    fread(&number, sizeof(uint32_t), 1, file);
	    fread(&size,   sizeof(uint32_t), 1, file);
	    track.resize(number); assert(size == sizeof(l1tpf_int::PropagatedTrack));
	    fread(&track[0], size, number, file);

	    fread(&number, sizeof(uint32_t), 1, file);
	    fread(&size,   sizeof(uint32_t), 1, file);
	    muon.resize(number); assert(size == sizeof(l1tpf_int::Muon));
	    fread(&muon[0], size, number, file);
	}
};

struct Event {
	uint32_t run, lumi; uint64_t event;
	float z0;
	float alphaCMed, alphaCRms, alphaFMed, alphaFRms;
	std::vector<Region> regions;
	
	Event() : run(0), lumi(0), event(0), z0(0.), alphaCMed(0.), alphaCRms(0.), alphaFMed(0.), alphaFRms(0.), regions() {}
	void readFromFile(FILE *fRegionDump) {
		fread(&run, sizeof(uint32_t), 1, fRegionDump);
		fread(&lumi, sizeof(uint32_t), 1, fRegionDump);
		fread(&event, sizeof(uint64_t), 1, fRegionDump);
		fread(&z0, sizeof(float), 1, fRegionDump);
		fread(&alphaCMed, sizeof(float), 1, fRegionDump);
		fread(&alphaCRms, sizeof(float), 1, fRegionDump);
		fread(&alphaFMed, sizeof(float), 1, fRegionDump);
		fread(&alphaFRms, sizeof(float), 1, fRegionDump);
		uint32_t nregions;
		fread(&nregions, sizeof(uint32_t), 1, fRegionDump);
		regions.resize(nregions);
		for (unsigned int i = 0; i < nregions; ++i) {
			regions[i].readFromFile(fRegionDump);
		}
	}
};

class DiscretePFInputs {
	public:
		DiscretePFInputs(const char *fileName) : file_(fopen(fileName,"rb")), iregion_(0) {}
		~DiscretePFInputs() { fclose(file_); }
		bool nextRegion(CaloObj calo[NCALO], TkObj track[NTRACK], z0_t & hwZPV) {
			while(true) {
				if (event_.event == 0 || iregion_ == event_.regions.size()) {
					if (feof(file_)) return false;
					event_.readFromFile(file_);
					printf("Beginning of run %u, lumi %u, event %lu \n", event_.run, event_.lumi, event_.event);
					iregion_ = 0;
				}
				const Region &r = event_.regions[iregion_];

				if (fabs(r.etaCenter) > 1.5) {
					iregion_++;
					continue; // use only regions in the barrel for now
				}

				hwZPV = event_.z0 * l1tpf_int::InputTrack::Z0_SCALE;
				for (unsigned int i = 0; i < std::min<unsigned>(NTRACK,r.track.size()); ++i) {
					track[i].hwPt = r.track[i].hwPt;
					track[i].hwPtErr = r.track[i].hwCaloPtErr;
					track[i].hwEta = r.track[i].hwEta; // @calo
					track[i].hwPhi = r.track[i].hwPhi; // @calo
					track[i].hwZ0 = r.track[i].hwZ0;
				}
				for (unsigned int i = 0; i < std::min<unsigned>(NCALO, r.calo.size()); ++i) {
					calo[i].hwPt = r.calo[i].hwPt;
					calo[i].hwEta = r.calo[i].hwEta;
					calo[i].hwPhi = r.calo[i].hwPhi;
				}
				printf("Read region %u with %lu tracks and %lu calo\n", iregion_, r.track.size(), r.calo.size());

				iregion_++;
				return true;
			}
		}
	private:
		FILE *file_;
		Event event_;
		unsigned int iregion_;
};
#endif
