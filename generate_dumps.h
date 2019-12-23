//
// generate_dumps.h
//
#include <iostream>
#include <string>
#include <string.h>     /* strlen */
#include <stdlib.h>     /* strtoul */
#include <vector>
#include <stdio.h>
#include <stdint.h>
#include <sstream>
#include <chrono>

#include "DiscretePFInputs.h"

using namespace std;


typedef l1tpf_int::InputRegion Region;
typedef l1tpf_int::CaloCluster Calo;
typedef l1tpf_int::Muon Muon;
typedef l1tpf_int::PropagatedTrack Track;

/* Calo GetCalo(int16_t hwPt, int16_t hwEta, int16_t hwPhi, */
/*              bool isEM=true, int16_t hwEmPt=0, int16_t hwPtErr=0){ */
/*     Calo x{hwPt,hwEmPt,hwPtErr,hwEta,hwPhi,0,isEM}; */
/*     return x; */
/* } */

#define NINTYPE 4
enum InputType{kCa,kEM,kTk,kMu};
#define NPFTYPE 5
enum PFType{kPFCh,kPFNe,kPFPh,kPFEl,kPFMu};

class Config {
 public:
    long num_events  = 1;

    // avg PF inputs/outputs
    float avg_nObj[NINTYPE];
    float avg_nPF[NPFTYPE];

    // generate random inputs
    inline int16_t Pt() {return (*rng_pt )(*rng);}
    inline int16_t Eta(){return (*rng_eta)(*rng);}
    inline int16_t Phi(){return (*rng_phi)(*rng);}
    
    // gen a poisson number of objects or fixed
    void Poisson();
    // generate events using PF inputs quantities or outputs
    // i.e. pfch -> one track and calo, etc...
    bool genFromPF = false;

    /* void Init(); */
    /* bool _initialized=false; */
    uint GetNCalo();
    uint GetNObj(InputType kObj);
    /* template<typename T> T GetObj(InputType kObj); */
    /* template<typename T> std::vector<T> GetObjs(InputType kObj); */
    template<typename T> void FillObj(T& obj, InputType kObj);
    template<typename T> void FillObjs(std::vector<T>& objs, InputType kObj);
    
    Config();
    ~Config();
 private:
    std::default_random_engine *rng;

    std::uniform_int_distribution<int> *rng_pt ;
    std::uniform_int_distribution<int> *rng_eta;
    std::uniform_int_distribution<int> *rng_phi;
    std::uniform_int_distribution<int> *rng_high_pt;
    
    bool _doPoisson;
    std::poisson_distribution<int> * rng_nObj[NINTYPE];
    /* std::poisson_distribution<int> *rng_ncalo  ; */
    /* std::poisson_distribution<int> *rng_nemcalo; */
    /* std::poisson_distribution<int> *rng_ntrack ; */
    /* std::poisson_distribution<int> *rng_nmuon  ; */
};
/* void Config::Init(){ */
/*     _initialized=true; */

/* } */

void Config::Poisson(){
    _doPoisson=true;
    for(int i=0;i< NINTYPE;i++)
        rng_nObj[i] = new std::poisson_distribution<int>(avg_nObj[i]);
}

Config::Config(){
    // PF Inputs
    avg_nObj[kCa] = 15.;
    avg_nObj[kEM] = 9.;
    avg_nObj[kTk] = 16.;
    avg_nObj[kMu] = 1.;
    // PF outputs
    avg_nPF[kPFCh] = 8.;
    avg_nPF[kPFNe] = 6.;
    avg_nPF[kPFPh] = 6.;
    avg_nPF[kPFEl] = 1.;
    avg_nPF[kPFMu] = 1.;

    //initialize random number services
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    rng = new std::default_random_engine(seed); // seed
    //rng = new std::default_random_engine(1776); // seed

    rng_pt  = new std::uniform_int_distribution<int>(1,100); // pt
    rng_eta = new std::uniform_int_distribution<int>(-270,270); // eta, with a bit beyond range of +/- 243
    rng_phi = new std::uniform_int_distribution<int>(-512,511); // phi val in 10 bits
    rng_high_pt = new std::uniform_int_distribution<int>(100,10*1000);

    /* if(genFromPF){ */
    /*     avg_ncalo   = avg_pfch+avg_pfne; */
    /*     avg_nemcalo = avg_pfel+avg_pfph+avg_pfch*0.3; */
    /*     avg_ntrack  = avg_pfch+avg_pfel+avg_pfmu; */
    /*     avg_nmuon   = avg_pfmu; */
    /* } */

    /* rng_nObj[kCa] = new std::poisson_distribution<int>(avg_nObj[kCa]); */
    /* rng_nObj[kEM] = new std::poisson_distribution<int>(avg_nObj[kEM]); */
    /* rng_nObj[kTk] = new std::poisson_distribution<int>(avg_nObj[kTk]); */
    /* rng_nObj[kMu] = new std::poisson_distribution<int>(avg_nObj[kMu]); */
}

Config::~Config(){
    /* if (!_initialized) return; */
    delete rng;
    delete rng_pt ;
    delete rng_eta;
    delete rng_phi;
    delete rng_high_pt;

    if(!_doPoisson) return;
    for(int i=0;i< NINTYPE;i++)
        delete rng_nObj[i];
    /* delete rng_ncalo  ; */
    /* delete rng_nemcalo; */
    /* delete rng_ntrack ; */
    /* delete rng_nmuon  ; */
}

uint Config::GetNObj(InputType kObj){
    if(_doPoisson) return (*rng_nObj[kObj])(*rng);
    else return uint(avg_nObj[kObj]);
}


template<typename T>
void Config::FillObj(T &obj, InputType kObj){
    int16_t pt=Pt();
    int16_t eta=Eta();
    int16_t phi=Phi();
    
    obj.hwPt = pt;
    obj.hwEta = eta;
    obj.hwPhi = phi;

    /* if (is_same<T,Calo>::value){ */
    /*     if (kObj==kCa){ */
    /*         obj.isEM=false; */
    /*         obj.hwEmPt=0; */
    /*         obj.hwPtErr=pt/8+1; */
    /*     } else { */
    /*         obj.isEM=true; */
    /*         obj.hwEmPt=pt; */
    /*         obj.hwPtErr=pt/16+1; */
    /*     } */
    /* } */
    /* else if (is_same<T,Track>::value){ */
    /*     obj.isEM=(kObj==kEM); */
    /*     obj.hwInvpt = pow(2*10)/pt; */
    /*     obj.hwVtxEta = eta; */
    /*     obj.hwVtxPhi = phi; */
    /*     obj.hwZ0 = phi/5; // don't remember the range off-hand */
    /*     obj.hwChi2 = 1000/pt;// just pick some convenient number */
    /*     obj.hwPtErr = pt/16+1; */
    /*     obj.hwCaloPtErr = pt/8+1; */
    /* } */
    return;
}

/* template<typename T> */
/* T Config::GetObj(InputType kObj){ */
/*     int16_t pt=Pt(); */
/*     int16_t eta=Eta(); */
/*     int16_t phi=Phi(); */
/*     switch (kObj) */
/*     case kCa: { */
/*         //Calo x{pt, 0, int16_t(pt/4), eta, phi, 0, false}; */
/*         T x{pt, 0, int16_t(pt/4), eta, phi, 0, false}; */
/*         return x; */
/*     } */
/*     case kEM: { */
/*         //Calo x{pt, 0, int16_t(pt/8), eta, phi, 0, true}; */
/*         T x{pt, 0, int16_t(pt/8), eta, phi, 0, true}; */
/*         return x; */
/*     } */
/*     /\* case kTk: { *\/ */
/*  case (is_same<T,Track>::value): { */
/*         //Track x{pt, int16_t(pt/16), int16_t(pt/16), Eta(), Phi()}; */
/*         //Track x;//{16, 1, 2, 23, 49}; */
/*         T x;//{16, 1, 2, 23, 49}; */
/*         x.hwInvpt = pow(2*10)/pt; */
/*         x.hwVtxEta = eta; */
/*         x.hwVtxPhi = phi; */
/*         x.hwZ0 = phi/5; // don't remember the range off-hand */
/*         x.hwChi2 = 1000/pt;// just pick some convenient number */
/*         x.hwPt = pt; */
/*         x.hwPtErr = pt/16; */
/*         x.hwCaloPtErr = pt/8; */
/*         x.hwEta = eta; */
/*         x.hwPhi = phi; */
/*         return x; */
/*     } */
/*     case kMu: { */
/*         //Muon x{pt, eta, phi}; */
/*         T x{pt, eta, phi}; */
/*         return x; */
/*     } */
/* } */

/* template<typename T> */
/* std::vector<T> Config::GetObjs(InputType kObj){ */
/*     std::vector<T> x{}; */
/*     int N = GetNObj(kObj); */
/*     cout << N << endl; */
/*     while (x.size()<N){x.push_back(GetObj<T>(kObj));} */
/*     return x; */
/* } */
template<typename T>
void Config::FillObjs(std::vector<T>& objs, InputType kObj){
    /* int N = ; */
    /* cout << N << endl; */
    /* for(int i=0;i<N;i++)  */
    objs.resize(GetNObj(kObj));
    for(auto& x : objs)
        FillObj(x,kObj);
}




/* uint Config::GetNCalo(){ */
/*     if(doPoisson) return (*rng_ncalo)(*rng); */
/*     else return avg_ncalo; */
/* } */
/* Calo Config::GetCalo(){ */
/*     int16_t  pt=Pt(); */
/*     Calo x{pt, 0, int16_t(pt/4), Eta(), Phi(), 0, false}; */
/*     return x; */
/* } */

