#ifndef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_H
#define FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_H

#if defined(__GXX_EXPERIMENTAL_CXX0X__) or defined(CMSSW)
#include <cstdint>
#define FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_MORE
#else
#include <stdint.h>
#endif

// the serialization may be hidden if needed
#define FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_IO

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>

namespace l1tpf_int { 

  struct CaloCluster {
      int16_t  hwPt;   
      int16_t  hwEmPt;   
      int16_t  hwPtErr;   
      int16_t  hwEta;   
      int16_t  hwPhi;   
      uint16_t hwFlags;
      bool     isEM, used;

      // sorting
      bool operator<(const CaloCluster &other) const { return hwPt > other.hwPt; }

#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_MORE
      static constexpr float PT_SCALE = 4.0;     // quantize in units of 0.25 GeV (can be changed)
      static constexpr float ETAPHI_FACTOR = 4;  // size of an ecal crystal in phi in integer units (our choice)
      static constexpr float ETAPHI_SCALE = ETAPHI_FACTOR*(180./M_PI);  // M_PI/180 is the size of an ECal crystal; we make a grid that is 4 times that size
      static constexpr int16_t PHI_WRAP = 360*ETAPHI_FACTOR;            // what is 3.14 in integer

      // filling from floating point
      void fill(float pt, float emPt, float ptErr, float eta, float phi, bool em, unsigned int flags) {
          hwPt  = round(pt  * CaloCluster::PT_SCALE);
          hwEmPt  = round(emPt  * CaloCluster::PT_SCALE);
          hwPtErr = round(ptErr  * CaloCluster::PT_SCALE);
          hwEta = round(eta * CaloCluster::ETAPHI_SCALE);
          hwPhi = int16_t(round(phi * CaloCluster::ETAPHI_SCALE)) % CaloCluster::PHI_WRAP;
          isEM  = em;
          used = false;
          hwFlags = flags;
      }

      float floatPt() const { return float(hwPt) / CaloCluster::PT_SCALE; }
      float floatEmPt() const { return float(hwEmPt) / CaloCluster::PT_SCALE; }
      float floatPtErr() const { return float(hwPtErr) / CaloCluster::PT_SCALE; }
      static float minFloatPt() { return float(1.0) / CaloCluster::PT_SCALE; }
      float floatEta() const { return float(hwEta) / CaloCluster::ETAPHI_SCALE; }
      float floatPhi() const { return float(hwPhi) / CaloCluster::ETAPHI_SCALE; }
      void  setFloatPt(float pt) { hwPt  = round(pt  * CaloCluster::PT_SCALE); }
      void  setFloatEmPt(float emPt) { hwEmPt  = round(emPt  * CaloCluster::PT_SCALE); }
#endif

#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_IO
      void writeToFile(FILE *file) const {
        fwrite(&hwPt, 2, 1, file);
        fwrite(&hwEmPt, 2, 1, file);
        fwrite(&hwPtErr, 2, 1, file);
        fwrite(&hwEta, 2, 1, file);
        fwrite(&hwPhi, 2, 1, file);
        fwrite(&hwFlags, 2, 1, file);
        fwrite(&isEM, 1, 1, file); 
        // used is not written out
      }
      void readFromFile(FILE *file) {
        fread(&hwPt, 2, 1, file);
        fread(&hwEmPt, 2, 1, file);
        fread(&hwPtErr, 2, 1, file);
        fread(&hwEta, 2, 1, file);
        fread(&hwPhi, 2, 1, file);
        fread(&hwFlags, 2, 1, file);
        fread(&isEM, 1, 1, file); 
        used = false;
      }
#endif
  };

  // https://twiki.cern.ch/twiki/bin/view/CMS/L1TriggerPhase2InterfaceSpecifications
  struct InputTrack {
      uint16_t hwInvpt;
      int32_t  hwVtxEta;
      int32_t  hwVtxPhi;
      bool     hwCharge;
      int16_t  hwZ0;
      uint16_t hwChi2, hwStubs;
      uint16_t hwFlags;

#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_MORE
      static constexpr float INVPT_SCALE   = 2E4;    // 1%/pt @ 100 GeV is 2 bits 
      static constexpr float VTX_PHI_SCALE = 1/2.5E-6; // 5 micro rad is 2 bits
      static constexpr float VTX_ETA_SCALE = 1/1E-5;   // no idea, but assume it's somewhat worse than phi
      static constexpr float Z0_SCALE      = 20;     // 1mm is 2 bits
      static constexpr int32_t VTX_ETA_1p3 = 1.3 * InputTrack::VTX_ETA_SCALE;

      // filling from floating point
      void fillInput(float pt, float eta, float phi, int charge, float dz, unsigned int flags) {
          hwInvpt  = round(1/pt  * InputTrack::INVPT_SCALE);
          hwVtxEta = round(eta * InputTrack::VTX_ETA_SCALE);
          hwVtxPhi = round(phi * InputTrack::VTX_PHI_SCALE);
          hwCharge = (charge > 0);
          hwZ0     = round(dz  * InputTrack::Z0_SCALE);
          hwFlags = flags;
      }

      float floatVtxPt() const { return 1/(float(hwInvpt) / InputTrack::INVPT_SCALE); }
      float floatVtxEta() const { return float(hwVtxEta) / InputTrack::VTX_ETA_SCALE; }
      float floatVtxPhi() const { return float(hwVtxPhi) / InputTrack::VTX_PHI_SCALE; }
      float floatDZ()     const { return float(hwZ0) / InputTrack::Z0_SCALE; }
      int intCharge()     const { return hwCharge ? +1 : -1; }
#endif

#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_IO
      void writeToFile(FILE *file) const {
        fwrite(&hwInvpt, 2, 1, file);
        fwrite(&hwVtxEta, 4, 1, file);
        fwrite(&hwVtxPhi, 4, 1, file);
        fwrite(&hwCharge, 1, 1, file); 
        fwrite(&hwZ0, 2, 1, file);
        fwrite(&hwChi2, 2, 1, file);
        fwrite(&hwStubs, 2, 1, file);
        fwrite(&hwFlags, 2, 1, file);
      }
      void readFromFile(FILE *file) {
        fread(&hwInvpt, 2, 1, file);
        fread(&hwVtxEta, 4, 1, file);
        fread(&hwVtxPhi, 4, 1, file);
        fread(&hwCharge, 1, 1, file); 
        fread(&hwZ0, 2, 1, file);
        fread(&hwChi2, 2, 1, file);
        fread(&hwStubs, 2, 1, file);
        fread(&hwFlags, 2, 1, file);
      }
#endif
  };

  struct PropagatedTrack : public InputTrack {
      int16_t  hwPt;
      int16_t  hwPtErr;
      int16_t  hwCaloPtErr;
      int16_t  hwEta; // at calo
      int16_t  hwPhi; // at calo
      bool     muonLink;
      bool     used; // note: this flag is not used in the default PF, but is used in alternative algos
      bool     fromPV;

      // sorting
      bool operator<(const PropagatedTrack &other) const { return hwPt > other.hwPt; }

#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_MORE
      void fillPropagated(float pt, float ptErr, float caloPtErr, float eta, float phi, unsigned int flags) {
          hwPt  = round(pt  * CaloCluster::PT_SCALE);
          hwPtErr = round(ptErr  * CaloCluster::PT_SCALE);
          hwCaloPtErr = round(caloPtErr  * CaloCluster::PT_SCALE);
          hwEta = round(eta * CaloCluster::ETAPHI_SCALE);
          hwPhi = int16_t(round(phi * CaloCluster::ETAPHI_SCALE)) % CaloCluster::PHI_WRAP;
          muonLink = false;
          used = false;
      }

      float floatPt() const { return float(hwPt) / CaloCluster::PT_SCALE; }
      float floatPtErr() const { return float(hwPtErr) / CaloCluster::PT_SCALE; }
      float floatCaloPtErr() const { return float(hwCaloPtErr) / CaloCluster::PT_SCALE; }
      float floatEta() const { return float(hwEta) / CaloCluster::ETAPHI_SCALE; }
      float floatPhi() const { return float(hwPhi) / CaloCluster::ETAPHI_SCALE; }
#endif

#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_IO
      void writeToFile(FILE *file) const {
        InputTrack::writeToFile(file);
        fwrite(&hwPt, 2, 1, file);
        fwrite(&hwPtErr, 2, 1, file);
        fwrite(&hwCaloPtErr, 2, 1, file);
        fwrite(&hwEta, 2, 1, file);
        fwrite(&hwPhi, 2, 1, file);
        // muonLink, used, fromPV are transient
      }
      void readFromFile(FILE *file) {
        InputTrack::readFromFile(file);
        fread(&hwPt, 2, 1, file);
        fread(&hwPtErr, 2, 1, file);
        fread(&hwCaloPtErr, 2, 1, file);
        fread(&hwEta, 2, 1, file);
        fread(&hwPhi, 2, 1, file);
        muonLink = false; used = false; fromPV = false;
      }
#endif


  };

  struct Muon {
      int16_t  hwPt;   
      int16_t  hwEta;   // at calo
      int16_t  hwPhi;   // at calo
      uint16_t hwFlags;
      bool     hwCharge;

      // sorting
      bool operator<(const Muon &other) const { return hwPt > other.hwPt; }

#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_MORE
      void fill(float pt, float eta, float phi, int charge, unsigned int flags) {
          // we assume we use the same discrete ieta, iphi grid for all particles 
          hwPt  = round(pt  * CaloCluster::PT_SCALE);
          hwEta = round(eta * CaloCluster::ETAPHI_SCALE);
          hwPhi = int16_t(round(phi * CaloCluster::ETAPHI_SCALE)) % CaloCluster::PHI_WRAP;
          hwCharge = (charge > 0);
          hwFlags = flags;
      }
      float floatPt() const { return float(hwPt) / CaloCluster::PT_SCALE; }
      float floatEta() const { return float(hwEta) / CaloCluster::ETAPHI_SCALE; }
      float floatPhi() const { return float(hwPhi) / CaloCluster::ETAPHI_SCALE; }
      int intCharge()     const { return hwCharge ? +1 : -1; }
#endif

#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_IO
      void writeToFile(FILE *file) const {
        fwrite(&hwPt, 2, 1, file);
        fwrite(&hwEta, 2, 1, file);
        fwrite(&hwPhi, 2, 1, file);
        fwrite(&hwFlags, 2, 1, file);
        fwrite(&hwCharge, 1, 1, file); 
      }
      void readFromFile(FILE *file) {
        fread(&hwPt, 2, 1, file);
        fread(&hwEta, 2, 1, file);
        fread(&hwPhi, 2, 1, file);
        fread(&hwFlags, 2, 1, file);
        fread(&hwCharge, 1, 1, file); 
      }
#endif
  };

  struct PFParticle {
      int16_t         hwPt;   
      int16_t         hwEta;  // at calo face 
      int16_t         hwPhi;   
      uint8_t         hwId; // CH=0, EL=1, NH=2, GAMMA=3, MU=4 
      int16_t         hwVtxEta;  // propagate back to Vtx for charged particles (if useful?)
      int16_t         hwVtxPhi;   
      uint16_t        hwFlags;
      CaloCluster     cluster;
      PropagatedTrack track;
      bool            chargedPV;
      uint16_t        hwPuppiWeight;
      uint16_t        hwStatus; // for debugging

      // sorting
      bool operator<(const PFParticle &other) const { return hwPt > other.hwPt; }

#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_MORE
      static constexpr float PUPPI_SCALE = 100;

      float floatPt() const { return float(hwPt) / CaloCluster::PT_SCALE; }
      float floatEta() const { return float(hwEta) / CaloCluster::ETAPHI_SCALE; }
      float floatPhi() const { return float(hwPhi) / CaloCluster::ETAPHI_SCALE; }
      float floatVtxEta() const { return (track.hwPt > 0 ? track.floatVtxEta() : float(hwVtxEta) / CaloCluster::ETAPHI_SCALE); }
      float floatVtxPhi() const { return (track.hwPt > 0 ? track.floatVtxPhi() : float(hwVtxPhi) / CaloCluster::ETAPHI_SCALE); }
      float floatDZ() const { return float(track.hwZ0) / InputTrack::Z0_SCALE; }
      float floatPuppiW() const { return float(hwPuppiWeight) / PUPPI_SCALE; }
      int intCharge()     const { return (track.hwPt > 0 ? track.intCharge() : 0); }
      void setPuppiW(float w) {
            hwPuppiWeight = std::round(w * PUPPI_SCALE);
      }
      void  setFloatPt(float pt) { hwPt  = round(pt  * CaloCluster::PT_SCALE); }
#endif
  };


#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_IO
  template<typename T>
  void writeManyToFile(const std::vector<T> & objs, FILE *file) {
    uint32_t number = objs.size(); 
    fwrite(&number, 4, 1, file);
    for (uint32_t i = 0; i < number; ++i) objs[i].writeToFile(file);
  }

  template<typename T>
  void readManyFromFile(std::vector<T> & objs, FILE *file) {
    uint32_t number;
    fread(&number, 4, 1, file);
    objs.resize(number); 
    for (uint32_t i = 0; i < number; ++i) objs[i].readFromFile(file);
  }
#endif

  struct InputRegion {
    float etaCenter, etaMin, etaMax, phiCenter, phiHalfWidth;
    float etaExtra, phiExtra;
    std::vector<CaloCluster>      calo;
    std::vector<CaloCluster>      emcalo;
    std::vector<PropagatedTrack>  track;
    std::vector<Muon>             muon;

    InputRegion() : etaCenter(), etaMin(), etaMax(), phiCenter(), phiHalfWidth(), etaExtra(), phiExtra() {}
    InputRegion(float etacenter, float etamin, float etamax, float phicenter, float phihalfwidth, float etaextra, float phiextra) :
        etaCenter(etacenter), etaMin(etamin), etaMax(etamax), phiCenter(phicenter), phiHalfWidth(phihalfwidth), etaExtra(etaextra), phiExtra(phiextra) {}

#ifdef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPFINPUTS_IO
    void writeToFile(FILE *file) const {
        assert(4 == sizeof(float));
        fwrite(&etaCenter, 4, 1, file);
        fwrite(&etaMin,    4, 1, file);
        fwrite(&etaMax,    4, 1, file);
        fwrite(&phiCenter, 4, 1, file);
        fwrite(&phiHalfWidth, 4, 1, file);
        fwrite(&etaExtra, 4, 1, file);
        fwrite(&phiExtra, 4, 1, file);
        writeManyToFile(calo, file);
        writeManyToFile(emcalo, file);
        writeManyToFile(track, file);
        writeManyToFile(muon, file);
    }
    void readFromFile(FILE *file) {
        assert(4 == sizeof(float));
        fread(&etaCenter, 4, 1, file);
        fread(&etaMin,    4, 1, file);
        fread(&etaMax,    4, 1, file);
        fread(&phiCenter, 4, 1, file);
        fread(&phiHalfWidth, 4, 1, file);
        fread(&etaExtra, 4, 1, file);
        fread(&phiExtra, 4, 1, file);
        readManyFromFile(calo, file);
        readManyFromFile(emcalo, file);
        readManyFromFile(track, file);
        readManyFromFile(muon, file);
    }
#endif
    
  };

} // namespace
#endif
