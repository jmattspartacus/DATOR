#ifndef LIBGRET_GRETINACONF_HH
#define LIBGRET_GRETINACONF_HH

#include <map>
#include <vector>

namespace Gret {
  enum class FirstX;

  enum class EffType {
    kLogLogPoly = 1,
    kRadWare = 2
  };
                      
  
  class Configuration  {
  public:
    std::map<int, std::pair<double, double> > Calibration;
    std::map<std::pair<int, int>, double[4][4]> CrystalMap;
    double PosXOffset;
    double PosYOffset;
    double PosZOffset;
    double BeamLeftOffset;
    double BeamRightOffset;
    double MeanT0;
    FirstX FX;

    double ClusterAngle;
    float EfficiencyCal[7];
    float AbsEff;
    EffType EfficiencyType;
    bool Addback;
    bool PAD128;
    bool PAD128Fix;
    bool FixValid;
    int CrystalSwap[120];
    int NoisyCrystals[120];
    int AddbackDets[120];
    bool warnings;
  public:
    Configuration() : PAD128(false),
                      PAD128Fix(false),
                      FixValid(false),
                      PosXOffset(0),
                      PosYOffset(0),
                      PosZOffset(0),
                      BeamLeftOffset(0),
                      BeamRightOffset(0),
                      warnings(false) {
      for (int i=0; i<120; ++i) { CrystalSwap[i] = i; }
      for (int i=0; i<120; ++i) { NoisyCrystals[i] = 0; }
      for (int i=0; i<120; ++i) { AddbackDets[i] = 1; }
    }
    Configuration(std::string name, std::string title) : Configuration() {
      //      SetName(name.c_str()); SetTitle(title.c_str());
    }
    ~Configuration() {}
    
    int ReadCrystalMap(std::string fn);
    int ReadCalibration(std::string fn);
    int ReadEfficiencyCal(std::string fn);
    int ReadCrystalSwap(std::string fn);
    int ReadNoisyCrystals(std::string fn);
    int ReadAddbackDets(std::string fn);
    std::vector<double> GetPosition(int hole, int xtl, double x, double y, double z) const;
    float Efficiency(const double &e) const;
    void PrintAngles();

  };
}

#endif
