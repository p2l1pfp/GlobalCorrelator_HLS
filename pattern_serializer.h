#include "src/data.h"
#include <cstdio>

class MP7PatternSerializer {
    public:
        MP7PatternSerializer(const std::string &fname, const std::string &boardName = "Board MP7_L1PF") ;
        ~MP7PatternSerializer() ;
        
        void operator()(const MP7DataWord event[MP7_NCHANN]) ;

    protected:
        const std::string fname_;
        FILE *file_;
        unsigned int ipattern_;
    
};

class HumanReadablePatternSerializer {
    public:
        HumanReadablePatternSerializer(const std::string &fname) ;
        ~HumanReadablePatternSerializer() ;
        
        void operator()(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK], const PFChargedObj outch[NTRACK], const PFNeutralObj outpho[NPHOTON], const PFNeutralObj outne[NSELCALO]) ;
        void dump_inputs(const EmCaloObj emcalo[NEMCALO], const HadCaloObj hadcalo[NCALO], const TkObj track[NTRACK]) ;
        void dump_outputs(const PFChargedObj outch[NTRACK], const PFNeutralObj outpho[NPHOTON], const PFNeutralObj outne[NSELCALO]) ;

    protected:
        const std::string fname_; 
        FILE *file_; // may be stdout
        unsigned int ipattern_;
    
};
