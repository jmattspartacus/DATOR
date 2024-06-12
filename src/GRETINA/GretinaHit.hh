#ifndef LIBGRET_GRETINAHIT_HH
#define LIBGRET_GRETINAHIT_HH

#include "GRETINA/RawData.hh"
#include "GRETINA/GretinaConf.hh"

namespace Gret {
  class Event;
  
  enum class FirstX {
    kSmallestZ = 1,
    kMeanPos = 2,
    kEnergyWeighted = 3,
    kCrystalCent = 4,
    kMaxEnergy = 5,
    kHybrid = 6
  };
    
  class GretinaHit {
  public:
    int Hole;
    int CrystalID;
    double RawEnergy; //un-calibrated
    double TotalEnergy; //calibrated
    int nInteractions;
    int PAD;
    float chisq;
    float t0;
    int64_t timestamp;
    bool Fix;

    double x;
    double y;
    double z;
    
    double PosX;
    double PosY;
    double PosZ;

    double Theta;
    double Phi;
    
    double Time;

    double Efficiency;

    bool valid;

    GretinaHit() : Fix(false), valid(true) {}
    int Build(const int64_t GEBtimestamp,
              const crys_intpts *data);

  };
}

#endif
