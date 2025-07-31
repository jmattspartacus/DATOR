#ifndef LIBORR_ORRUBA_HH
#define LIBORR_ORRUBA_HH

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

#include "Reader/Reader.hh"
#include "ORRUBA/SingleParticle.hh"
#include "ORRUBA/SX3.hh"
#include "ORRUBA/QQQ5.hh"
#include "ORRUBA/BB10.hh"
#include "ORRUBA/TDC.hh"

namespace Orruba { 

  class Kinematics {
  public:
    double m_beam;
    double m_targ;
    double m_eject; //light
    double m_rec;   //heavy

    double m_beam_MeV;
    double m_targ_MeV;
    double m_eject_MeV; //light
    double m_rec_MeV;   //heavy
      
    double beam_en; //kinetic
    double beam_v;
    double com_v;
    double e_in;
    double Q_gs;

    double p_beam;
    double p_eject;
      
    void setup() {
      beam_v = std::sqrt(2.0 * beam_en/m_beam);
      com_v = beam_v * (m_beam)/(m_beam + m_targ);
      e_in = 0.5 * m_beam * std::pow(beam_v - com_v, 2) + 0.5 * m_targ * std::pow(com_v,2);

      //convert masses into energies in MeV
      m_beam_MeV = m_beam * 931.478;
      m_targ_MeV = m_targ * 931.478;
      m_eject_MeV = m_eject * 931.478;
      m_rec_MeV = m_rec * 931.478;

      p_beam = std::sqrt(beam_en*beam_en + 2.0*m_beam*beam_en);
    }
      
    double getExRel(double energy, double theta) {        

      double p_eject = std::sqrt(energy*energy + 2.0*m_eject_MeV*energy);
      double e_beam = beam_en + m_beam_MeV; //total energy
      double e_eject = m_eject_MeV + energy; //total energy 

      double e_rec = e_beam+m_targ_MeV-e_eject;

      p_beam = std::sqrt(beam_en*beam_en + 2.0*m_beam_MeV*beam_en);

      double relQval = m_beam_MeV + m_targ_MeV - m_eject_MeV
        - std::sqrt(m_beam_MeV*m_beam_MeV + m_targ_MeV*m_targ_MeV + m_eject_MeV*m_eject_MeV + 2.0*m_targ_MeV*e_beam -
                    2.0*e_eject*(e_beam + m_targ_MeV) + 2.*p_beam*p_eject*std::cos(theta));

      return (Q_gs-relQval);
    }

    double getExNonRel(double energy, double theta) {
      double vel = std::sqrt(2.0*energy/m_eject);
      double velCoM2 = com_v*com_v + vel*vel - 2.0*com_v*vel*std::cos(theta);
      double e_com = 0.5*(m_eject + m_eject*m_eject/m_rec)*velCoM2;
      double nrQval = e_com - e_in;
      return (Q_gs-nrQval);
    }
  };

  class Configuration {
  public:
    std::vector<DetType> types; //per ADC channel
    std::vector<int> detID;
    std::vector<int> layer;
    std::vector<int> side; //front/back
    std::vector<int> subID;
    std::vector<int> LR;
    std::vector<int> UDS;  //up/down stream
    std::vector<float> pedestal;
    std::vector<float> offset;
    std::vector<float> gain;
    std::vector<float> threshold;

    std::vector<sx3cal> sx3cals;
    std::vector<qqq5cal> qqq5cals;
    std::vector<bb10cal> bb10cals;

    float radii[32];
    float dEthetas[32];
    float Ethetas[32];
    float dEdist = -76.0;
    float Edist = -80.0;

    float r_dE_strip[4];
    float r_dE_strip_bb10[8];
    float r_E_strip[4];
    float r_dE_barrel = 98.8;
    float r_E_barrel = 100.2;
    float sx3width = 40.0;
    float sx3zoffset = 2.0; //active area from z=0 offset

    Configuration();
    Configuration(std::string filename);
    Configuration(std::string name, std::string title, std::string filename);
    void Set(std::string name, std::string title, std::string filename);
    ~Configuration();
    void ReadConfiguration(std::string filename);
    void ReadCalibration(std::string filename);
    void ReadPositionCalibration(std::string filename);
    void ReadRadii(std::string filename);
    void SetThresholds(float thresh);
    void SetThresholds(int chan, float thresh);
    void SetThresholds(DetType type, float thresh);
    void SetThresholds(DetType type, int s, float thresh);

  };
      
  class TrackerHit  {
  public:
    float value;
    float cal;
    float tdc;

    bool valid;
    unsigned short int ID;

    TrackerHit() {}
    TrackerHit(unsigned short int v,
               unsigned short int c) : value(v), cal(c) {}

    int Validate();

  };
      
  class Tracker  {
  public:
    int ID;
    //std::vector<TrackerHit> xhits;
    //std::vector<TrackerHit> yhits;
    TrackerHit xhits[16];
    TrackerHit yhits[16];

    int presentX;
    int presentY;
    int nGoodX;
    int nGoodY;
    bool validX;
    bool validY;
    float cathode_val;
    float cathode;

    //these are in wires
    float x;
    float y;

    float xcal;  //these are in the gretina frame of reference
    float ycal;  //+x = down, +y = BR 

    void SetX();
    void SetY();
      
    float GetX() { return x; };
    float GetY() { return y; };

    float GetXCal() { return xcal; }
    float GetYCal() { return ycal; }

    Tracker() {
      for (int i=0; i<16; ++i) {
        xhits[i].ID = i;
        yhits[i].ID = i;
      }
    }
    void Set(unsigned short int chan, unsigned short int value);
    void Reset() {
      for (int i=0; i<16; ++i) {
        xhits[i].valid = false;
        xhits[i].value = 0;
        xhits[i].tdc = 0;
        yhits[i].valid = false;
        yhits[i].value = 0;
        yhits[i].tdc = 0;
      }
      cathode = 0; ID = -1; x=-1; y=-1; }
    void Print(std::ostream &out);

  };
  
  class Event : public DATOR::Processor {
  public:
    static Configuration conf;
    
    std::vector<SX3> sx3s;
    std::vector<QQQ5> qqq5s;
    std::vector<BB10> bb10s;
    Tracker tracker; //! prevent this from being written to TTree
    TDC tdc;

    //don't like that this is a vector of pointers... means we'll be dynamically allocating each event which is expensive
    std::vector<SingleParticle*> single_parts; //these come from a single SX3 or QQQ5, with front/back together
      
    //std::vector<Particle> parts; //these are matched dE/E

    static float dither;
    static float Dither() { dither += 0.1; if (dither >=1.0) { dither = 0.0; } return dither; }
      
    unsigned long long timestamp;
    std::vector<unsigned short int> chans;
    std::vector<unsigned short int> vals;
    //unsigned short int chans[4096];
    //unsigned short int vals[4096];
    int nhits;

    //per event counters
    int sx3hits;
    int qqq5hits;
    int bb10hits;
    int sx3parts;
    int qqq5parts;
    int bb10parts;

    bool SpuriousMyRIAD;     //! don't need this written to the tree

    //counters that do not get reset
    unsigned long long nQQQ5Hits;
    unsigned long long nSX3Hits;
    unsigned long long nBB10Hits;
    //hits above threshold
    unsigned long long nQQQ5HitsTh;
    unsigned long long nSX3HitsTh;
    unsigned long long nBB10HitsTh;
    //reconstructed particles
    unsigned long long nQQQ5Particles;
    unsigned long long nSX3Particles;
    unsigned long long nBB10Particles;
    unsigned long long nSpuriousMyRIAD;
    //events with at least one hit above threshold
    unsigned long long nQQQ5Evts;
    unsigned long long nSX3Evts;
    unsigned long long nBB10Evts;
    //events with at least one hit but no reconstructed event
    unsigned long long nBadQQQ5Evts;
    unsigned long long nBadSX3Evts;
    unsigned long long nBadBB10Evts;

    Event() { sx3s.reserve(24); qqq5s.reserve(8); bb10s.reserve(16); chans.reserve(4096); vals.reserve(4096); dither = 0.0; 
    nQQQ5Hits=0;
    nSX3Hits=0;
    nBB10Hits=0;
    nQQQ5HitsTh=0;
    nSX3HitsTh=0;
    nBB10HitsTh=0;
    nQQQ5Particles=0;
    nSX3Particles=0;
    nBB10Particles=0;
    nQQQ5Evts=0;
    nSX3Evts=0;
    nBB10Evts=0;
    nBadQQQ5Evts=0;
    nBadSX3Evts=0;
    nBadBB10Evts=0;
    nSpuriousMyRIAD=0;
    }
    
    void Set();

    int AddQQQ5(unsigned short int channel, unsigned short int value);
    int AddSX3(unsigned short int channel, unsigned short int value);
    int AddBB10(unsigned short int channel, unsigned short int value);
    void Reset() { 
      sx3s.clear(); qqq5s.clear(); bb10s.clear(); tracker.Reset(); tdc.Reset(); 
      for (int i=0; i<single_parts.size(); ++i) {
        delete single_parts[i];
      }
      single_parts.clear();
      chans.clear();
      vals.clear();
      nhits = 0;
      sx3hits = 0;
      qqq5hits = 0;
      bb10hits = 0;
      sx3parts = 0;
      qqq5parts = 0;
      bb10parts = 0;
    }
    ~Event() { Reset(); }

    void Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length);
    void ProcessFinal() {};
    
    void SetConf(Configuration c) { conf = c; }
    void PrintSummary(std::ostream &out);

  };

  class Basic : public DATOR::Processor {
  public:
    unsigned long long timestamp;
    int nhits;
    std::vector<unsigned short int> chans;
    std::vector<unsigned short int> vals;

    Basic() { chans.reserve(4096); vals.reserve(4096); }

    void Reset() { chans.clear(); vals.clear(); nhits = 0; }
    void Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length);
    void ProcessFinal() {};
  };
}  

#endif
