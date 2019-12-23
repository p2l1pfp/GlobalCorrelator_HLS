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
    inline int IntRange(float a, float b, float f){return std::round(a+(b-a)*f);}
    inline int16_t Pt() {return IntRange(  10, 100, (*rng_f)(*rng));}
    inline int16_t Eta(){return IntRange(-270, 270, (*rng_f)(*rng));}
    inline int16_t Phi(){return IntRange(-512, 511, (*rng_f)(*rng));}
    inline float Rand() {return (*rng_f)(*rng);}
    inline float Gaus(float mu, float sig){return mu+sig*((*rng_gaus)(*rng));}
    
    void FillAllObjects(std::vector<Calo>  &calos,
                        std::vector<Calo>  &emcalos,
                        std::vector<Track> &tracks, 
                        std::vector<Muon>  &muons);

    
    // gen a poisson number of objects or fixed
    void Poisson();
    // generate events using PF inputs quantities or outputs
    // i.e. pfch -> one track and calo, etc...
    void GenFromPF();

    // fill objects with random inputs
    uint GetNPF(PFType kObj);
    uint GetNObj(InputType kObj);
    template<typename T> void FillObj(T& obj, InputType kObj);
    template<typename T> void FillObjs(std::vector<T>& objs, InputType kObj);
    
    Config();
    ~Config();
 private:
    std::default_random_engine *rng;
    std::uniform_real_distribution<float> *rng_f;
    std::normal_distribution<float> *rng_gaus;
    
    bool _doPoisson=false;
    std::poisson_distribution<int> * rng_nObj[NINTYPE];
    std::poisson_distribution<int> * rng_nPF[NPFTYPE];

    bool _genFromPF=false;
    template<typename T1, typename T2> void Link(T1 &to, T2 &from);
    float _PF_EMprob=0.3; // prob for a hadron to leave energy in EM also
    float _PF_EMshare=0.2; // prob for a hadron to leave energy in EM also
};

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
    rng_f  = new std::uniform_real_distribution<float>(0,1); // utility
    rng_gaus  = new std::normal_distribution<float>(); // defaults to unit
}

Config::~Config(){
    /* if (!_initialized) return; */
    delete rng;
    delete rng_f;
    delete rng_gaus;

    if(!_doPoisson) return;
    for(int i=0;i< NINTYPE;i++)
        delete rng_nObj[i];
    for(int i=0;i< NPFTYPE;i++)
        delete rng_nPF[i];
}

void Config::Poisson(){
    _doPoisson=true;
    for(int i=0;i< NINTYPE;i++)
        rng_nObj[i] = new std::poisson_distribution<int>(avg_nObj[i]);
    for(int i=0;i< NPFTYPE;i++)
        rng_nPF[i] = new std::poisson_distribution<int>(avg_nPF[i]);
}

void Config::GenFromPF(){
    _genFromPF=true;
}

template<typename T1, typename T2>
void Config::Link(T1 &to, T2 &from){
    to.hwEta = from.hwEta;
    to.hwPhi = from.hwPhi;
}

void Config::FillAllObjects(std::vector<Calo>  &calos, 
                            std::vector<Calo>  &emcalos,
                            std::vector<Track> &tracks,
                            std::vector<Muon>  &muons){
    if(!_genFromPF){
        FillObjs(calos, kCa);
        FillObjs(emcalos, kEM);
        FillObjs(tracks, kTk);
        FillObjs(muons, kMu);
        return;
    }
    
    // generate from PF objects
    for(int i=0;i<GetNPF(kPFCh);i++){ // charged
        Track t; FillObj(t, kTk);
        Calo c; Link(c, t);
        if (Rand()>_PF_EMprob){ // em calo contribution
            Calo cEM; Link(cEM, t);
            c.hwPt   = Gaus((1.-_PF_EMshare)*t.hwPt, 0.1*t.hwPt);
            cEM.hwPt = Gaus(_PF_EMshare*t.hwPt, 0.1*t.hwPt);
            cEM.isEM = true;
            emcalos.push_back(cEM); 
        } else {
            c.hwPt = Gaus(t.hwPt, 0.1*t.hwPt);
        }        
        tracks.push_back(t);
        calos.push_back(c); 
    }
    for(int i=0;i<GetNPF(kPFNe);i++){ // neutral hadron
        Calo c; FillObj(c,kCa); 
        if (Rand()>_PF_EMprob){ // em calo contribution
            Calo cEM; Link(cEM, c);
            cEM.hwPt = Gaus(_PF_EMshare*c.hwPt, 0.1*c.hwPt);
            cEM.isEM = true;
            emcalos.push_back(cEM); 
        }
        calos.push_back(c);         
    }
    for(int i=0;i<GetNPF(kPFPh);i++){ // photon
        Calo c; FillObj(c,kCa);
        c.isEM = true;
        emcalos.push_back(c);
    }
    for(int i=0;i<GetNPF(kPFCh);i++){ // electron
        Track t; FillObj(t, kTk);
        Calo c; Link(c, t);
        c.hwPt   = Gaus(t.hwPt, 0.1*t.hwPt);
        c.isEM = true;
        tracks.push_back(t);
        emcalos.push_back(c); 
    }
    for(int i=0;i<GetNPF(kPFMu);i++){ // muon
        Track t; FillObj(t, kTk);
        Muon c; Link(c, t);
        c.hwPt   = Gaus(t.hwPt, 0.1*t.hwPt);
        tracks.push_back(t);
        muons.push_back(c); 
    }

}




uint Config::GetNPF(PFType kObj){
    if(_doPoisson) return (*rng_nPF[kObj])(*rng);
    else return uint(avg_nPF[kObj]);
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

template<typename T>
void Config::FillObjs(std::vector<T>& objs, InputType kObj){
    /* int N = ; */
    /* cout << N << endl; */
    /* for(int i=0;i<N;i++)  */
    objs.resize(GetNObj(kObj));
    for(auto& x : objs)
        FillObj(x,kObj);
}



