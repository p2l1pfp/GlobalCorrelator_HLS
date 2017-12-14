#include <cstdio>
#include <vector>
#include "../firmware/data.h"

class MP7PatternSerializer {
    public:
        /** 
            if nmux > 1, we assume that each event is sent in nmux clock cycles using nchann = MP7_NCHANN/nmux channels.
            there is no time offsetting in the inputs, i.e. the first nmux events all start at clock 0 end end at clock nmux-1
            data is sent by column, i.e. the first tmux words are sent on the first fiber in tmux clock cycles
            suppose tmux = 4 and event[i] = [ w0[i] w1[i] ... w71[i] ]
            then    payload[t = 0] = [ w0[0] w4[0] ... w64[0] w68[0]    w0[1] w4[1] ... w64[1] w68[1]    w0[2] ... w68[2]     w0[3] ... w68[3] ]
                    payload[t = 1] = [ w1[0] w5[0] ... w65[0] w69[0]    w1[1] w5[1] ... w65[1] w69[1]    w1[2] ... w69[2]     w0[3] ... w69[3] ]
                    payload[t = 2] = [ w2[0] w6[0] ... w66[0] w70[0]    w2[1] w6[1] ... w66[1] w70[1]    w2[2] ... w70[2]     w0[3] ... w70[3] ]
                    payload[t = 3] = [ w3[0] w7[0] ... w67[0] w71[0]    w3[1] w7[1] ... w67[1] w71[1]    w3[2] ... w71[2]     w0[3] ... w71[3] ]
                    payload[t = 0] = [ w0[4] w4[4] ... w64[4] w68[4]    w0[5] w4[5] ... w64[5] w68[5]    w0[6] ... w68[6]     w0[7] ... w68[7] ]
                    payload[t = 1] = [ w1[4] w5[4] ... w65[4] w69[4]    w1[5] w5[5] ... w65[5] w69[5]    w1[6] ... w69[6]     w0[7] ... w69[7] ]
                    payload[t = 2] = [ w2[4] w6[4] ... w66[4] w70[4]    w2[5] w6[5] ... w66[5] w70[5]    w2[6] ... w70[6]     w0[7] ... w70[7] ]
                    payload[t = 3] = [ w3[4] w7[4] ... w67[4] w71[4]    w3[5] w7[5] ... w67[5] w71[5]    w3[6] ... w71[6]     w0[7] ... w71[7] ]

            if nempty > 0, put more N frames of empty events in the pattern file in between each good event.
            e.g. nempty = 1, nmux = 1 will give
                  | event 1 ......    |
                  | zeros             |
                  | event 2 ......    |
                  | zeros             |
            while nempty = 1, nmux = 2 will give
                  | event 1 (part 1)...   zero .... |
                  | event 1 (part 2)...   zero .... |
                  | event 2 (part 1)...   zero .... |
                  | event 2 (part 2)...   zero .... |

            if nempty < 0, will do like nempty > 0 but filling with "magic" data instead of zeros

        */
        MP7PatternSerializer(const std::string &fname, unsigned int nmux=1, int nempty=0, unsigned int nlinks=0, const std::string &boardName = "Board MP7_L1PF") ;
        ~MP7PatternSerializer() ;
        
        void operator()(const MP7DataWord event[MP7_NCHANN]) ;
        
    protected:
        const std::string fname_;
        const unsigned int nlinks_, nmux_, nchann_, nempty_;
        const bool fillmagic_;
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class HumanReadablePatternSerializer {
    public:
        HumanReadablePatternSerializer(const std::string &fname) ;
        ~HumanReadablePatternSerializer() ;
        
        void operator()(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU], const PFChargedObj outch[NTRACK], const PFNeutralObj outpho[NPHOTON], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU]) ;
        void dump_inputs(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const MuObj mu[NMU]) ;
        void dump_emcalo(const EmCaloObj emcalo[NEMCALO], unsigned int N = NEMCALO) ;
        void dump_hadcalo(const HadCaloObj hadcalo[NCALO], unsigned int N = NCALO) ;
        void dump_track(const TkObj track[NTRACK], unsigned int N = NTRACK) ;
        void dump_mu(const MuObj mu[NMU], unsigned int N = NMU) ;
        void dump_outputs(const PFChargedObj outch[NTRACK], const PFNeutralObj outpho[NPHOTON], const PFNeutralObj outne[NSELCALO], const PFChargedObj outmu[NMU]) ;

    protected:
        const std::string fname_; 
        FILE *file_; // may be stdout
        unsigned int ipattern_;
    
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CTP7PatternSerializer {
    
    public:
        CTP7PatternSerializer(const std::string &fname, unsigned int nchann_max, bool isInput, unsigned int nmux=1, int nempty=0, const std::string &boardName = "Board_CTP7_L1PF") ;
        ~CTP7PatternSerializer() ;
        
        void operator()(const MP7DataWord event[MP7_NCHANN], unsigned int nchann);
        
    protected:
        
        const std::string fname_;
        const unsigned int nmux_, nchann_, nempty_;
        const bool fillmagic_;
        FILE *file_;
        unsigned int ipattern_;
        bool isInput_;
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

        template<typename T> void print(unsigned int iframe, const T & event, unsigned int nchann);
        // void push(const MP7DataWord event[MP7_NCHANN]);
        // void flush();
        // void zero();
    
};

