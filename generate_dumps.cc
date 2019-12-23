//
// Generate randomized input dump files
// g++ generate_dumps.cc -std=c++11
//
#include <random>

#include "generate_dumps.h"
using namespace std;


// template<typename T> 
// std::vector<T> GetObjs(Config& c, InputType kObj){
//     std::vector<T> x{};
//     int N = c.GetNObj(kObj);
//     while (x.size()<N){x.push_back(c.GetObj(kObj));}
//     return x;
// }
// std::vector<Calo> GetCalos(Config& c){
//     std::vector<Calo> x{};
//     int N = c.GetNCalo();
//     while (x.size()<N){x.push_back(c.GetCalo());}
//     return x;
// }
// std::vector<Calo> GetEMCalos(Config& c){
//     std::vector<Calo> x{};
//     return x;
// }
// std::vector<Track> GetTracks(Config& c){
//     std::vector<Track> x{};
//     return x;
// }
// std::vector<Muon> GetMuons(Config& c){
//     std::vector<Muon> x{};
//     return x;
// }


void WriteEvent(Config& c,FILE *f, uint64_t event){
    // See Event::readFromFile in DiscretePFInputs_IO.h

    uint32_t run = 1;
    uint32_t lumi = 1000;
    fwrite(&run, sizeof(uint32_t), 1, f);
    fwrite(&lumi, sizeof(uint32_t), 1, f);
    fwrite(&event, sizeof(uint64_t), 1, f);

    // add region
    uint32_t one = 1; // 1 reg per event
    fwrite(&one, 4, 1, f); 
    float etaCenter = 0.;
    float etaMin = -1.5;
    float etaMax = 1.5;
    float phiCenter = 0.;
    float phiHalfWidth = 3.14159;
    float etaExtra = 0.;
    float phiExtra = 0.;
    Region r{etaCenter, etaMin, etaMax, phiCenter, phiHalfWidth,
            etaExtra, phiExtra};
    // fill objects according to config, add to region
    std::vector<Calo>  calos  ; 
    std::vector<Calo>  emcalos; 
    std::vector<Track> tracks ; 
    std::vector<Muon>  muons  ; 
    c.FillAllObjects(calos, emcalos, tracks, muons);

    cout << "Event counts: calos "
         << calos.size()  << ", EM "
         << emcalos.size() << ", tracks "
         << tracks.size() << ", muons "
         << muons.size() << std::endl;
    r.calo   = calos  ;
    r.emcalo = emcalos;
    r.track  = tracks ;
    r.muon   = muons  ;
    r.writeToFile(f);

    // vertex / puppi info
    float z0 = 0.;
    fwrite(&z0, sizeof(float), 1, f);
    fwrite(&z0, sizeof(float), 1, f);
    float alpha = 1.;
    fwrite(&alpha, sizeof(float), 1, f);
    fwrite(&alpha, sizeof(float), 1, f);
    fwrite(&alpha, sizeof(float), 1, f);
    fwrite(&alpha, sizeof(float), 1, f);
}

void WriteFile(Config& c, const char* fname){
    //initialize rng service if necessary
    // c.Init();

    // write file
    FILE *f = fopen(fname, "wb");
    for(int evt_index=1;evt_index<=c.num_events;evt_index++){
        WriteEvent(c, f, evt_index);
    }
    fclose(f);

}


int main(int argc, char **argv)
{
    //setup random number generator
    std::default_random_engine generator(1776); // seed
    std::uniform_real_distribution<float> pt_dist(10.,100.); 
    std::poisson_distribution<int> test_dist(10.5); //mean
    float x2 = pt_dist(generator);
 
    // write standard events
    Config c;
    c.num_events=10;
    c.Poisson();
    WriteFile(c,"data/test.dump");

    // write events based on PF objects
    c.GenFromPF();
    WriteFile(c,"data/testPF.dump");

    // l1tpf_int::CaloCluster x{1,2,3};
    // //    x.fill(10.,5.,2.,1,1,false,0);
    // cout << x.hwPt << endl;
    // cout << x.hwEmPt << endl;
    // cout << x.hwPtErr<< endl;
    // cout << x.hwEta  << endl;
    // cout << x.hwPhi  << endl;
    return 0;

}
