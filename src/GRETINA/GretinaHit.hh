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
    double SegEnergy[36];
    int SegID[36];
    int PAD;
    float chisq;
    float t0;
    int64_t timestamp;
    bool BadIntE;
    bool BadT0;

    double x;
    double y;
    double z;
    
    double PosX;
    double PosY;
    double PosZ;

    double Theta;
    double Phi;

    double MaxIntEn;
    
    long long int Time;

    double Efficiency;

    bool valid;

    GretinaHit() : BadIntE(false), BadT0(false), valid(true) {}
    int Build(const int64_t GEBtimestamp,
              const crys_intpts *data);
    bool operator<(GretinaHit &other) { return TotalEnergy < other.TotalEnergy; }

  };
}

#endif
