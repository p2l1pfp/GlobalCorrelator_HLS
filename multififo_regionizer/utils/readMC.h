#ifndef multififo_regionizer_readMC_h
#define multififo_regionizer_readMC_h

bool readEventTk(FILE *file, std::vector<TkObj> inputs[NTKSECTORS], uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventTk(FILE *file, std::vector<TkObj> inputs[NTKSECTORS][NTKFIBERS], uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventCalo(FILE *file, std::vector<HadCaloObj> inputs[NCALOSECTORS*NCALOFIBERS], bool zside, uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventCalo(FILE *file, std::vector<HadCaloObj> inputs[NCALOSECTORS][NCALOFIBERS], bool zside, uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventMu(FILE *file, std::vector<GlbMuObj> &inputs, uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventMu(FILE *file, std::vector<GlbMuObj> inputs[NMUFIBERS], uint32_t &run, uint32_t &lumi, uint64_t &event) ;
bool readEventVtx(FILE *file, std::vector<std::pair<z0_t,pt_t>> & inputs, uint32_t &irun, uint32_t &ilumi, uint64_t &ievent) ;
 
#endif
