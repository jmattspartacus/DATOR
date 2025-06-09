#ifndef LIBGRET_GRETINACONF_HH
#define LIBGRET_GRETINACONF_HH

#include <map>
#include <vector>

namespace Gret {
  enum class AddbackType {
    NoAddback,
    Cluster2Xtl,
    ClusterAnyXtl,
    Nearest2Xtl,
    NearestAnyXtl
  };

  enum class FirstX;

  enum class EffType {
    kLogLogPoly = 1,
    kRadWare = 2
  };
                      
  class Array {
  public:
    int dim;
    double *data;
    Array();
    Array(Array& other);
    Array(int d);
    ~Array();
    void Allocate(int d);
    void Set(int i, int j, double val);
    double Get(int i, int j) const;
    void operator=(Array& other);
  };
  
  class Configuration  {
  public:
    //std::map<int, std::pair<double, double> > Calibration;
    double CalGain[124];
    double CalOffset[124];
    
    //std::map<std::pair<int, int>, double[4][4]> CrystalMap;
    Array CrystalMap[124];
    double PosXOffset;
    double PosYOffset;
    double PosZOffset;
    double BeamLeftOffset;
    double BeamRightOffset;
    double BeamLeftZOffset;
    double BeamRightZOffset;
    double MeanT0;
    FirstX FX;

    double ClusterAngle;
    double AddbackTDiff;
    float EfficiencyCal[7];
    float AbsEff;
    EffType EfficiencyType;
    AddbackType Addback; //0 = none, 1 = max 2 crystals cluster addback, 2 = any crystals cluster addback
    bool PAD128;
    bool PAD128Fix;
    bool FixValid;
    int CrystalSwap[124];
    int NoisyCrystals[124];
    int AddbackDets[124];
    double CrystalTheta[124];
    double CrystalPhi[124];
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
                      BeamLeftZOffset(0),
                      BeamRightZOffset(0),
                      AddbackTDiff(400.0),
                      warnings(false) {
      for (int i=0; i<124; ++i) {
        CrystalSwap[i] = i;
        NoisyCrystals[i] = 0; 
        AddbackDets[i] = 1;
        CalGain[i] = 1;
        CalOffset[i] = 0;
      }
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
