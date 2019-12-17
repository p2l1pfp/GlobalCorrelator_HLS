#include "pattern_serializer.h"
#include <cassert>
#include <string>
#include <cstdlib>
#include <cassert>

#if defined(PACKING_DATA_SIZE) and defined(PACKING_NCHANN)

PatternSerializer::PatternSerializer(const std::string &fname, unsigned int nmux, unsigned int nzero, bool zero_valid, unsigned int nprefix, unsigned int npostfix, const std::string &boardName) :
    fname_(fname), nin_(PACKING_NCHANN), nout_(PACKING_NCHANN/nmux), nmux_(nmux), nzero_(nzero), nprefix_(nprefix), npostfix_(npostfix), zerovalid_(zero_valid), file_(nullptr), ipattern_(0) 
{
    if (!fname.empty()) {
        const unsigned int extra_space = (PACKING_DATA_SIZE-32)/4;
        char extra_spacer[extra_space+1];
        std::fill(extra_spacer, &extra_spacer[extra_space], ' ');
        extra_spacer[extra_space] = '\0';
        file_ = fopen(fname.c_str(), "w");
        fprintf(file_, "Board %s\n", boardName.c_str());
        fprintf(file_, " Quad/Chan :    ");
        for (unsigned int i = 0; i < nlinks_; ++i) fprintf(file_, "q%02dc%1d%s      ", i/4, i % 4, extra_spacer);
        fprintf(file_, "\n      Link :     ");
        for (unsigned int i = 0; i < nlinks_; ++i) fprintf(file_, "%02d%s         ", i, extra_spacer);
        fprintf(file_, "\n");
    }
    if (nmux_ > 1) {
        assert(PACKING_NCHANN % nmux_ == 0);
    }

    if (nprefix_ || npostfix_ || nzero_) {
        zeroframe_.resize(nout_);
        for (unsigned int j = 0; j < nout_; ++j) zeroframe_[j] = 0;
    }
    for (unsigned int j = 0; j < nprefix_; ++j) {
      print(zeroframe_, false);
    }
}

PatternSerializer::~PatternSerializer() 
{
    for (unsigned int j = 0; j < npostfix_; ++j) {
      print(zeroframe_, false);
    }
    if (file_) {
        fclose(file_); file_ = nullptr;
        printf("Saved %u patterns to %s.\n", ipattern_, fname_.c_str());
    }
}

void PatternSerializer::operator()(const ap_uint<PACKING_DATA_SIZE> event[PACKING_NCHANN]) 
{
    if (!file_) return;
    if (nmux_ == 1) {
    for (unsigned int j = 0; j < nmux_; ++i) {
        print(event, true, j, nmux_);
    }
    for (unsigned int j = 0; j < nzero_; ++j) {
      print(zeroframe_, zerovalid_);
    }
}

template<typename T> void PatternSerializer::print(unsigned int iframe, const T & event, bool valid, unsigned int ifirst, unsigned int stride) 
{
    assert(PACKING_DATA_SIZE == 32 || PACKING_DATA_SIZE == 64);

    fprintf(file_, "Frame %04u :", iframe);
    for (unsigned int i = 0, j = ifirst; i < nlinks_; ++i, j += stride) {
#if PACKING_DATA_SIZE == 32
        fprintf(file_, " %dv%08x", int(valid), event[j].to_uint32());
#else 
        fprintf(file_, " %dv%016llx", int(valid), event[j].to_uint64());
#endif
    }
    fprintf(file_, "\n");
    ipattern_++;
}

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

HumanReadablePatternSerializer::HumanReadablePatternSerializer(const std::string &fname, bool zerosuppress) :
    fname_(fname), zerosuppress_(zerosuppress), file_(nullptr), ipattern_(0) 
{
    if (!fname.empty()) {
        if (fname == "-") {
            file_ = stdout;
        } else {
            file_ = fopen(fname.c_str(), "w");
        }
    }
}

HumanReadablePatternSerializer::~HumanReadablePatternSerializer() 
{
    if (file_ && (file_ != stdout)) {
        fclose(file_); 
        printf("Saved %u human readable patterns to %s.\n", ipattern_, fname_.c_str());
    }
}

bool HumanReadablePatternSerializer::startframe() {
    if (!file_) return false;
    fprintf(file_, "Frame %04u:\n", ipattern_);
}
void HumanReadablePatternSerializer::endframe() {
    fprintf(file_, "\n");
    if (file_ == stdout) fflush(file_);
    ipattern_++;
}
void HumanReadablePatternSerializer::operator()(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], const PFChargedObj outch[NTRACK], const PFNeutralObj outpho[NPHOTON], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU]) 
{
    if (!startframe()) return;
    dump_inputs(emcalo,hadcalo,track,mu);
    dump_outputs(outch,outpho,outne,outmu);
    endframe();
}

void HumanReadablePatternSerializer::operator()(const HadCaloObj calo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], const PFChargedObj outch[NTRACK], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU]) 
{
    if (!startframe()) return;
    dump_inputs(calo,track,mu);
    dump_outputs(outch,outne,outmu);
    endframe();
}

void HumanReadablePatternSerializer::operator()(const PFChargedObj inch[NTRACK], const PFChargedObj inem[NPHOTON], const PFChargedObj inne[NSELCALO], const PFChargedObj inmu[NMU], unsigned int DATA_SIZE, const PFChargedObj outpart[/*DATA_SIZE*/]) 
{
    if (!startframe()) return;
    fprintf(file_, "Frame %04u:\n", ipattern_);
    dump_puppi(NTRACK,    "in ch", inch); // FIXME: do we want to dump PF or Puppi here?
    dump_puppi(NPHOTON,   "in em", inem); 
    dump_puppi(NSELCALO,  "in ne", inne);
    dump_puppi(NMU,       "in mu", inmu);
    dump_puppi(DATA_SIZE, "out  ", outpart); // FIXME: do we want to dump PF or Puppi here?
    endframe();
}


void HumanReadablePatternSerializer::operator()(unsigned int DATA_SIZE, const PFChargedObj outpart[/*DATA_SIZE*/], unsigned int NTAU, const PFChargedObj outtau[/*NTAU*/]) 
{
    if (!startframe()) return;
    fprintf(file_, "Frame %04u:\n", ipattern_);
    dump_puppi(DATA_SIZE, "input puppi", outpart);
    dump_puppi(NTAU,      "output taus", outtau); // FIXME: do we want to dump PF or Puppi here?
    endframe();
}

void HumanReadablePatternSerializer::dump_inputs(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU]) {
    dump_hadcalo(hadcalo);
    dump_emcalo(emcalo);
    dump_track(track);
    dump_mu(mu);
}
void HumanReadablePatternSerializer::dump_inputs(const HadCaloObj calo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU]) {
    dump_hadcalo(calo);
    dump_track(track);
    dump_mu(mu);
}

void HumanReadablePatternSerializer::dump_outputs(const PFChargedObj outch[NTRACK], const PFNeutralObj outpho[NPHOTON], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU]) 
{
    dump_pf(NTRACK,   "charged pf", outch);
    dump_pf(NPHOTON,  "photon  pf", outpho);
    dump_pf(NSELCALO, "neutral pf", outne);
    dump_pf(NMU,      "muon    pf", outmu);
}

void HumanReadablePatternSerializer::dump_outputs(const PFChargedObj outch[NTRACK], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU]) 
{
    dump_pf(NTRACK,   "charged pf", outch);
    dump_pf(NSELCALO, "neutral pf", outne);
    dump_pf(NMU,      "muon    pf", outmu);
}


void HumanReadablePatternSerializer::dump_hadcalo(const HadCaloObj hadcalo[NCALO], unsigned int N) {
    for (int i = 0; i < N; ++i) {
        if (zerosuppress_ && !hadcalo[i].hwPt) continue;
        fprintf(file_, "   calo  %3d, hwPt % 7d   hwEmPt  % 7d    hwEta %+7d   hwPhi %+7d   hwIsEM %1d\n", i, int(hadcalo[i].hwPt), int(hadcalo[i].hwEmPt), int(hadcalo[i].hwEta), int(hadcalo[i].hwPhi), int(hadcalo[i].hwIsEM));
    }
    if (file_ == stdout) fflush(file_);
}
void HumanReadablePatternSerializer::dump_emcalo(const EmCaloObj emcalo[NEMCALO], unsigned int N) {
    for (int i = 0; i < N; ++i) {
        if (zerosuppress_ && !emcalo[i].hwPt) continue;
        fprintf(file_, "   em    %3d, hwPt % 7d   hwPtErr % 7d    hwEta %+7d   hwPhi %+7d\n", i, int(emcalo[i].hwPt), int(emcalo[i].hwPtErr), int(emcalo[i].hwEta), int(emcalo[i].hwPhi));
    }
    if (file_ == stdout) fflush(file_);
}
void HumanReadablePatternSerializer::dump_track(const TkObj track[NTRACK], unsigned int N) {
    for (int i = 0; i < N; ++i) {
        if (zerosuppress_ && !track[i].hwPt) continue;
        fprintf(file_, "   track %3d, hwPt % 7d   hwPtErr % 7d    hwEta %+7d   hwPhi %+7d     hwZ0 %+7d\n", i, int(track[i].hwPt), int(track[i].hwPtErr), int(track[i].hwEta), int(track[i].hwPhi), int(track[i].hwZ0));
    }
    if (file_ == stdout) fflush(file_);
}
void HumanReadablePatternSerializer::dump_mu(const MuObj mu[NMU], unsigned int N) {
    for (int i = 0; i < N; ++i) {
        if (zerosuppress_ && !mu[i].hwPt) continue;
        fprintf(file_, "   muon  %3d, hwPt % 7d   hwPtErr % 7d    hwEta %+7d   hwPhi %+7d\n", i, int(mu[i].hwPt), int(mu[i].hwPtErr), int(mu[i].hwEta), int(mu[i].hwPhi));
    }
    if (file_ == stdout) fflush(file_);
}

void HumanReadablePatternSerializer::dump_pf(unsigned int N, const char *label, const PFChargedObj outch[/*N*/]) 
{
    for (int i = 0; i < N; ++i) {
        if (zerosuppress_ && !outch[i].hwPt) continue;
        fprintf(file_, "   %s %3d, hwPt % 7d   hwEta %+7d   hwPhi %+7d   hwId %1d      hwZ0 %+7d\n", label, i, 
                int(outch[i].hwPt), int(outch[i].hwEta), int(outch[i].hwPhi), int(outch[i].hwId), int(outch[i].hwZ0));
    }
}
void HumanReadablePatternSerializer::dump_pf(unsigned int N, const char *label, const PFNeutralObj outne[/*N*/]) 
{
    for (int i = 0; i < N; ++i) {
        if (zerosuppress_ && !outne[i].hwPt) continue;
        fprintf(file_, "   %s %3d, hwPt % 7d   hwEta %+7d   hwPhi %+7d   hwId %1d\n", label, i,
                int(outne[i].hwPt), int(outne[i].hwEta), int(outne[i].hwPhi), int(outne[i].hwId));
    }
    if (file_ == stdout) fflush(file_);
}
void HumanReadablePatternSerializer::dump_puppi(unsigned int N, const char *label, const PFChargedObj outpuppi[/*N*/]) 
{
    //FIXME: PFCharged doesn't have a Puppi PT
    for (int i = 0; i < N; ++i) {
        if (zerosuppress_ && !outpuppi[i].hwPt) continue; 
        fprintf(file_, "   %s %3d, hwPtPuppi % 7d   hwPt % 7d   hwEta %+7d   hwPhi %+7d   hwId %1d      hwZ0 %+7d\n", label, i, 
                int(outpuppi[i].hwPt), int(outpuppi[i].hwPt), int(outpuppi[i].hwEta), int(outpuppi[i].hwPhi), int(outpuppi[i].hwId), int(outpuppi[i].hwZ0));
    }
}
void HumanReadablePatternSerializer::dump_puppi(unsigned int N, const char *label, const PFNeutralObj outpuppi[/*N*/]) 
{
    for (int i = 0; i < N; ++i) {
        if (zerosuppress_ && !outpuppi[i].hwPtPuppi) continue;
        fprintf(file_, "   %s %3d, hwPtPuppi % 7d   hwPt % 7d   hwEta %+7d   hwPhi %+7d   hwId %1d\n", label, i,
                int(outpuppi[i].hwPtPuppi), int(outpuppi[i].hwPt), int(outpuppi[i].hwEta), int(outpuppi[i].hwPhi), int(outpuppi[i].hwId));
    }
    if (file_ == stdout) fflush(file_);
}


