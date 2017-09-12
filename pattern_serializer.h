#include "src/data.h"
#include <cstdio>
#include <vector>

class MP7PatternSerializer {
    public:
        /** 
            if nmux > 1, we assume that each event is sent in nmux clock cycles using nchann = MP7_NCHANN/nmux channels.
            there is no time offsetting in the inputs, i.e. the first nmux events all start at clock 0 end end at clock nmux-1
              at clock = 0,      channel[ 0      ...   nchann-1 ] = first nchann words of region (or event) 1
              at clock = 0,      channel[ nchann ... 2*nchann-1 ] = first nchann words of region (or event) 2
              at clock = 0,      channel .....
              at clock = 1,      channel[ 0      ...   nchann-1 ] = second nchann words of region (or event) 1
              at clock = 1,      channel[ nchann ... 2*nchann-1 ] = second nchann words of region (or event) 2
              at clock = 1,      channel .....
              at clock = nmux-1, channel[ 0      ...   nchann-1 ] = last nchann words of region (or event) 1
        */
        MP7PatternSerializer(const std::string &fname, unsigned int nmux=1, const std::string &boardName = "Board MP7_L1PF") ;
        ~MP7PatternSerializer() ;
        
        void operator()(const MP7DataWord event[MP7_NCHANN]) ;
        
    protected:
        const std::string fname_;
        const unsigned int nmux_, nchann_;
        FILE *file_;
        unsigned int ipattern_;
        class Pattern {
            public:
                Pattern() {}
                Pattern(const Pattern & other) { for (unsigned int i = 0; i < MP7_NCHANN; ++i) words[i] = other.words[i]; }
                MP7DataWord & operator[](int i) { return words[i]; }
                const MP7DataWord & operator[](int i) const { return words[i]; }
                void operator=(const Pattern & other) { for (unsigned int i = 0; i < MP7_NCHANN; ++i) words[i] = other.words[i]; }
            private:
                MP7DataWord words[MP7_NCHANN];
        };
        std::vector<Pattern> buffer_; // for muxing; holds the next patterns in output format. will fill nmux events, first, then print them all out

        template<typename T> void print(unsigned int iframe, const T & event);
        void push(const MP7DataWord event[MP7_NCHANN]);
        void flush();
        void zero();
    
};

class HumanReadablePatternSerializer {
    public:
        HumanReadablePatternSerializer(const std::string &fname) ;
        ~HumanReadablePatternSerializer() ;
        
        void operator()(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], const PFChargedObj outch[NTRACK], const PFNeutralObj outpho[NPHOTON], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU]) ;
        void dump_inputs(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU]) ;
        void dump_outputs(const PFChargedObj outch[NTRACK], const PFNeutralObj outpho[NPHOTON], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU]) ;

    protected:
        const std::string fname_; 
        FILE *file_; // may be stdout
        unsigned int ipattern_;
    
};
