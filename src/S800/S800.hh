#ifndef LIBS800_S800_HH
#define LIBS800_S800_HH

#include <iostream>
#include <fstream>
#include <functional>

#include "Reader/Reader.hh"

#include "S800/RawData.hh"
#include "S800/S800Confs.hh"

namespace S800 {
    
  class InvMap {
  public:
    class Block {
    public:
      int number;
      int ncoeffs[INVMAP_MAX_ORDER];
      double coeffs[INVMAP_MAX_ORDER][INVMAP_MAX_COEFFS];
      int exponents[INVMAP_MAX_ORDER][INVMAP_MAX_COEFFS][INVMAP_MAX_EXPONENTS];

      Block();
    };

    Block blocks[INVMAP_NBLOCKS];

    InvMap();
    int ReadInvMap(std::string fn);
    double Eval(const int &block, const int &order, double *x);
    void Print(std::ostream &out, int order);
  };

    
  class Configuration {
  public:
    Configuration() {};

    bool warning = false;

    //these are to align e1dn with e1up as a function of crdc2 x position
    double e1dn_slope = 1.;
    double e1dn_offset = 0.0;

    //these are to correct the obj-E1 TOF
    double obj_x1_slope = 1.0;
    double obj_x1_offset = 0.0;

    double obj_y1_slope = 1.0;
    double obj_y1_offset = 0.0;

    double obj_x2_slope = 1.0;
    double obj_x2_offset = 0.0;

    double obj_y2_slope = 1.0;
    double obj_y2_offset = 0.0;      
      
    double obj_afp_slope = 1.0;
    double obj_afp_offset = 0.0;

    double obj_bfp_slope = 1.0;
    double obj_bfp_offset = 0.0;

    double obj_xfp_slope = 1.0;
    double obj_xfp_offset = 0.0;

    //these are to correct the xfp-E1 TOF

    double xfp_afp_slope = 1.0;
    double xfp_afp_offset = 0.0;

    double xfp_x_slope = 1.0;
    double xfp_x_offset = 0.0;
      
      
    ICConf ionChamber;
    CRDCConf crdc;
    MTDCConf mtdc;
    InvMap invMap;
      
  };

  class Trigger  {
  public:
    unsigned short int trigPattern;

    unsigned short int trigTimes[4];
    int nTrigTimes;

    bool GetS800();
    bool GetCoinc();
    bool GetExt1();
    bool GetExt2();
    bool GetSecondary();

    void Print(std::ostream &out);

  };

  class TimeOfFlight  {
  public:
    int nTimes;
    unsigned short int chans[8];
    unsigned short int times[8];
    void Print(std::ostream &out);

    unsigned short int obj;
    unsigned short int xfp;

  };

    
  class Scintillator  {
  public:
    int nScints;
    unsigned short int chans[3];
    unsigned short int energies[3];
    unsigned short int times[3];
    void Print(std::ostream &out);
      
      
  };

  class IonChamber  {
  public:
    int nChans;
    unsigned short int chans[16];
    unsigned short int energies[16];
    double cal[16];
    void Print(std::ostream &out);

    double dE;
    double dE_corr;
    double dE_corr2;

    bool insideROI;

  };

  class CRDC  {
  public:
    int nHits;
    //unsigned short int padEnergies[S800_CRDC_MAX_HITS];

    std::vector<double> padEnergies;
    std::vector<unsigned short int> samples;
    std::vector<unsigned short int> pads;
    /*
      double padEnergies[S800_CRDC_MAX_HITS]; //the pedestal subtraction + gain matching is applied straight away
      unsigned short int samples[S800_CRDC_MAX_HITS]; 
      unsigned short int pads[S800_CRDC_MAX_HITS];
    */

    unsigned short int nSamples[S800_CRDC_NPADS]; // this will be zeroed every event
    //unsigned short int sampleNums[S800_CRDC_NPADS][S800_CRDC_MAX_SAMPLE];
    //double energies[S800_CRDC_NPADS][S800_CRDC_MAX_SAMPLE];
      
    unsigned short int anodeEnergy;
    unsigned short int anodeTime;
    int label;

    double gainx, gainy;
    double offsetx, offsety;

    double x;
    double y;
    double xcal;
    double ycal;

    int maxVal, hitInd, maxPad;

    CRDC();

    void Reset();
      
    void MaxHit(int &maxVal, int &hitInd, int &maxPad);
      
    double SetX();
    double SetY();
    double SetXCal();
    double SetYCal();
      
    double GetX();
    double GetY();
    double GetXCal();
    double GetYCal();

    void SetCal(double gx, double ox, double gy, double oy) { gainx=gx; gainy=gy; offsetx=ox; offsety=oy; }

    void Print(std::ostream &out);

  };

  class MTDC  {
  public:
    //int nHits[S800_MTDC_MAXCHANS];
    //int hitNum[S800_MTDC_MAXCHANS][S800_MTDC_MAXHITS];
    //int times[S800_MTDC_MAXCHANS][S800_MTDC_MAXHITS];
    int nHits;
    int mults[S800_MTDC_MAXCHANS]; //multiplcity per channel
    std::vector<int> chans;
    std::vector<int> hitNum;
    std::vector<int> times;

    int validTimes[S800_MTDC_MAXCHANS];

    int trigPattern;

    MTDC(); 
    void Reset();
    void Print(std::ostream &out);
    void Validate(MTDCConf *conf);      

  };

  class Event : public DATOR::Processor {
  public:
    int nHits = 0;
    int timeMult = 0;
    int evtNumMult = 0;
    int trigMult = 0;
    int tofMult = 0;
    int scintMult = 0;
    int icMult = 0;
    int crdcMult = 0;
    int crdc1Mult = 0;
    int crdc2Mult = 0;
    int mtdcMult = 0;

    bool valid;
    int validCode; //0 = valid
    //1-6 = invalid multiplicity
    //100-103 = invalid MTDC time (100+x)
    //21 = invalid e1 up/dn time difference
      
    unsigned long long headerTime;
    unsigned long long time;
    unsigned long long eventNumber;

    unsigned long long int nS800Hits;
    unsigned long long int nValidS800;
    unsigned long long int nS800Hits_IC;
    unsigned long long int nValidS800_IC;

    //these are the packets of the data stream
    Trigger trig;
    TimeOfFlight tof;
    Scintillator scint;
    IonChamber ionChamber;
    CRDC crdcs[2];
    MTDC mtdc;

    //these are "adopted" values after processing
    double e1up;
    double e1dn;
    double obj;
    double xfp;
    double rf;
    double obj_corr;
    double xfp_corr;

    static Configuration conf;

    Event() {};
    void ResetCounters() { nS800Hits = 0; nValidS800 = 0; nS800Hits_IC = 0; nValidS800_IC = 0; }
    void Reset();    
    void Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length);
    void ProcessFinal();
    void PrintSummary(std::ostream &out);
      
    int ReadTime(short unsigned int *data);
    int ReadEventNumber(short unsigned int *data);
    int ReadTrigger(short unsigned int *data, int len);
    int ReadTOF(short unsigned int *data, int len);
    int ReadScint(short unsigned int *data, int len);
    int ReadIonChamber(short unsigned int *data, int len);
    int ReadCRDC(short unsigned int *data, int len);
    int ReadMTDC(short unsigned int *data, int len);
    void Print(std::ostream &out);
    void SetConf(Configuration c);
    int GetAFP(double &theta);
    int GetBFP(double &theta);
    int Validate();

  };


  /*
    class Hodoscope {

    };

    class TPPACs {

    };

    class ObjPin {

    };

    class FocalPin {

    };

    class Galotte {

    };

    class LaBr {

    };

  */
      
}

#endif
