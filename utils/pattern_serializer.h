#ifndef UTILS_PATTERNSERIALIZER_H
#define UTILS_PATTERNSERIALIZER_H

#include <cstdio>
#include <vector>
#include "../firmware/data.h"


#if defined(PACKING_DATA_SIZE) && defined(PACKING_NCHANN)
class PatternSerializer {
    public:
        /** 
         *  each event corresponds to N=PACKING_NCHANN words of PACKING_DATA_SIZE bits
         *
         *  each event is serialized as <NMUX> frames, where the first <NMUX> words are send on the first link, etc.
         *  e.g. if event[i] = [ w0[i] w1[i] ... w(N-1)[i] ] then
         *  at NMUX = 1, the output uses N channels
         *      payload[t = 0] = [ w0[0] w1[0] ... w(N-1)[0] ]
         *      payload[t = 1] = [ w0[1] w1[1] ... w(N-1)[1] ]
         *  at NMUX = 2, the output uses N/2 channels
         *      payload[t = 0] = [ w0[0] w2[0] ... w(N-2)[0] ]
         *      payload[t = 1] = [ w1[0] w3[0] ... w(N-1)[0] ]
         *      payload[t = 2] = [ w0[1] w2[1] ... w(N-2)[1] ]                       
         *      payload[t = 3] = [ w1[1] w3[1] ... w(N-1)[1] ]
         *  at NMUX = 4, the output uses N/4 channels
         *      payload[t = 0] = [ w0[0] w4[0] ... w(N-4)[0] ]   
         *      payload[t = 1] = [ w1[0] w5[0] ... w(N-3)[0] ]
         *      payload[t = 2] = [ w2[0] w6[0] ... w(N-2)[0] ]   
         *      payload[t = 3] = [ w3[0] w7[0] ... w(N-1)[0] ]
         *      payload[t = 4] = [ w0[1] w4[1] ... w(N-4)[1] ]   
         *      payload[t = 5] = [ w1[1] w5[1] ... w(N-3)[1] ]
         *
         * In addition, it's possible to insert <NZERO> empty frames of zeros after the last frame of an event.
         * The "valid" bit in these patterns can be set to (either 0 or 1, depending on <ZERO_VALID>)
         * NMUX = 1, NZERO = 1
         *      payload[t = 0] = [ w0[0] w1[0] ... w(N-1)[0] ]
         *      payload[t = 1] = [  0     0    ...     0     ]
         *      payload[t = 2] = [ w0[1] w1[1] ... w(N-1)[1] ]
         *      payload[t = 3] = [  0     0    ...     0     ]
         * NMUX = 1, NZERO = 2
         *      payload[t = 0] = [ w0[0] w1[0] ... w(N-1)[0] ]
         *      payload[t = 1] = [  0     0    ...     0     ]
         *      payload[t = 2] = [  0     0    ...     0     ]
         *      payload[t = 3] = [ w0[1] w1[1] ... w(N-1)[1] ]
         *      payload[t = 4] = [  0     0    ...     0     ]
         *      payload[t = 5] = [  0     0    ...     0     ]
         * NMUX = 2, NZERO = 1
         *      payload[t = 0] = [ w0[0] w2[0] ... w(N-2)[0] ]
         *      payload[t = 1] = [ w1[0] w3[0] ... w(N-1)[0] ]
         *      payload[t = 2] = [  0     0    ...     0     ]
         *      payload[t = 3] = [ w0[1] w2[1] ... w(N-2)[1] ]                       
         *      payload[t = 4] = [ w1[1] w3[1] ... w(N-1)[1] ]
         *      payload[t = 5] = [  0     0    ...     0     ]
         *
         *
         * The serializer can also add <NPREFIX> frames of zero values with valid bit off at the beginning and <NPOSTFIX> frames and end of the file
         *
        */

        typedef ap_uint<PACKING_DATA_SIZE> Word;

        PatternSerializer(const std::string &fname, unsigned int nmux=1, unsigned int nzero=0, bool zero_valid=true, unsigned int nprefix=0, unsigned int npostfix=0, const std::string &boardName = "Board L1PF") ;
        ~PatternSerializer() ;
        
        void operator()(const Word event[PACKING_NCHANN], bool valid=true) ;
    
        template<typename T> void print(const T & event, bool valid = true, unsigned int ifirst = 0, unsigned int stride = 1);
        
    protected:
        const std::string fname_;
        const unsigned int nin_, nout_, nmux_, nzero_, nprefix_, npostfix_;
        const bool zerovalid_;
        FILE *file_;
        unsigned int ipattern_;
        std::vector<Word> zeroframe_;
};
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class HumanReadablePatternSerializer {
    public:
        HumanReadablePatternSerializer(const std::string &fname, bool zerosuppress=false) ;
        ~HumanReadablePatternSerializer() ;
        
        void operator()(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], const PFChargedObj outch[NTRACK], const PFNeutralObj outpho[NPHOTON], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU]) ;
        void operator()(const HadCaloObj calo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], const PFChargedObj outch[NTRACK], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU]) ;
        // layer 2 examples // FIXME: why PFCharged for the neutrals? 
        void operator()(const PFChargedObj inch[NTRACK], const PFChargedObj inem[NPHOTON], const PFChargedObj inne[NSELCALO], const PFChargedObj inmu[NMU], unsigned int DATA_SIZE, const PFChargedObj outpart[/*DATA_SIZE*/]) ;
        void operator()(unsigned int DATA_SIZE, const PFChargedObj outpart[/*DATA_SIZE*/], unsigned int NTAU, const PFChargedObj outtau[/*NTAU*/]) ;
        // PF all-in-one
        void dump_inputs(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU]) ;
        void dump_inputs(const HadCaloObj calo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU]) ;
        void dump_outputs(const PFChargedObj outch[NTRACK], const PFNeutralObj outpho[NPHOTON], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU]) ;
        void dump_outputs(const PFChargedObj outch[NTRACK], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU]) ;
        // Individual pieces
        void dump_emcalo(const EmCaloObj emcalo[NEMCALO], unsigned int N = NEMCALO) ;
        void dump_hadcalo(const HadCaloObj hadcalo[NCALO], unsigned int N = NCALO) ;
        void dump_track(const TkObj track[NTRACK], unsigned int N = NTRACK) ;
        void dump_mu(const MuObj mu[NMU], unsigned int N = NMU) ;
        void dump_pf(unsigned int N, const char *label, const PFChargedObj outch[/*N*/]) ;
        void dump_pf(unsigned int N, const char *label, const PFNeutralObj outch[/*N*/]) ;
        void dump_puppi(unsigned int N, const char *label, const PFChargedObj outpuppi[/*N*/]) ;
        void dump_puppi(unsigned int N, const char *label, const PFNeutralObj outpuppi[/*N*/]) ;
        bool startframe();
        void endframe();
    protected:
        const std::string fname_;
        bool zerosuppress_;
        FILE *file_; // may be stdout
        unsigned int ipattern_;
    
};

#endif

