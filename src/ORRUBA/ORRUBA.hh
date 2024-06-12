#ifndef LIBORR_ORRUBA_HH
#define LIBORR_ORRUBA_HH

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

#include "Reader/Reader.hh"

namespace Orruba {
  
  class SX3;
  class QQQ5;

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
    
  class sx3cal  {
  public:
    float padoff[4][4];
    float padgain[4][4];

    float stripoff[2][4][4];
    float stripgain[2][4][4];

    float stripPosGain[4];
    float stripPosOffset[4];
    
    sx3cal() {
      for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
          padoff[i][j] = 0.0;
          padgain[i][j] = 1.0;
          for (int k=0; k<1; ++k) {
            stripoff[k][i][j] = 0.0;
            stripgain[k][i][j] = 1.0;
          }
        }
          
        stripPosGain[i] = 1.0;
        stripPosOffset[i] = 0.0;
      }
    }

  };

  class qqq5cal  {
  public:
    float ringoff[4][32];
    float ringgain[4][32];
    
    float secoff[4][32];
    float secgain[4][32];

    qqq5cal() {
      for (int i=0; i<4; ++i) {
        for (int j=0; j<32; ++j) {
          ringoff[i][j] = 0.0;
          ringgain[i][j] = 1.0;

          secoff[i][j] = 0.0;
          secgain[i][j] = 1.0;    
        }
      }
    }

  };

  class Configuration {
  public:
    std::vector<std::string> types; //per ADC channel
    std::vector<int> detID;
    std::vector<int> layer;
    std::vector<int> side; //front/back
    std::vector<int> subID;
    std::vector<int> LR;
    std::vector<int> UDS;  //up/down stream
    std::vector<float> pedestal;
    std::vector<float> offset;
    std::vector<float> gain;

    std::vector<sx3cal> sx3cals;
    std::vector<qqq5cal> qqq5cals;

    float radii[32];
    float dEthetas[32];
    float Ethetas[32];
    float dEdist = -76.0;
    float Edist = -80.0;

    float r_dE_strip[4];
    float r_E_strip[4];
    float r_dE_barrel = 98.8;
    float r_E_barrel = 100.2;
    float sx3width = 40.0;

    Configuration();
    Configuration(std::string filename);
    Configuration(std::string name, std::string title, std::string filename);
    void ReadConfiguration(std::string filename);
    void ReadCalibration(std::string filename);
    void ReadPositionCalibration(std::string filename);
    void ReadRadii(std::string filename);


  };

  class SX3FrontHit  {
  public:
    unsigned short int side; //0=L, 1=R
    unsigned short int ID; //0-3 for strips 1-4
    float value;  //pedestal subtracted
    float cal; //calibrated (basic)

    SX3FrontHit() {}
    SX3FrontHit(unsigned short int s, unsigned short int id, float val, float c) : side(s), ID(id), value(val), cal(c) {}

  };

  class SX3BackHit  {
  public:
    unsigned short int ID; //0-3 for pad 1-4
    float value;  //pedestal subtracted
    float cal; //calibrated (basic)

    SX3BackHit() {}
    SX3BackHit(unsigned short int id, float val, float c) : ID(id), value(val), cal(c) {}

  };

  class SingleParticle  {
  public:
    unsigned int detType; //0=QQQ5, 1=SX3
    unsigned int detID;
    unsigned int frontID; //strip/ring
    unsigned int backID; //pad/sector
    unsigned int layer; //dE/E

    float frontEnergy;
    float backEnergy;
      
    bool valid;

    //coordinate system is gretina system: +z along beam axis, +x towards floor, +y towards beam left
    float r; //radius
    float x,y,z; //cartesian
    float theta; //in spherical coordinates, angle from +z axis
    float phi; //azimuth, angle from +x towards +y

    float r_off;
    float x_off, y_off, z_off;
    float theta_off;
    float phi_off;

    SingleParticle() {}
    SingleParticle(unsigned short int dt, unsigned short int did,
                   unsigned short int fid, unsigned short int bid,
                   unsigned short int lay,
                   float fe, float be, bool val) :
      detType(dt), detID(did), backID(bid), frontID(fid), layer(lay),
      frontEnergy(fe), backEnergy(be), valid(val) {};

    void OffsetBeam(float beamx, float beamy, bool verbose=false);
    
  };
      
  class SX3Particle : public SingleParticle {
  public:
    //some extra SX3 only things
    float stripL;
    float stripR;
    float position; //-1 -> 1;

    SX3Particle() {}
    SX3Particle(unsigned short int dt, unsigned short int did,
                unsigned short int pid, unsigned short int sid,
                unsigned short int lay,
                float pe, float sl, float sr, float pos, bool val) :
      SingleParticle(dt, did, sid, pid, lay, sl + sr, pe, val),
      stripL(sl), stripR(sr), position(pos) {}
    void MakeCoords(SX3 *detector);


  };

  class QQQ5Particle : public SingleParticle {
  public:
      
    void MakeCoords(QQQ5 *detector);
    QQQ5Particle() {}
    QQQ5Particle(unsigned short int dt, unsigned short int did,
                 unsigned short int rid, unsigned short int sid,
                 unsigned short int lay,
                 float re, float se, bool val) :
      SingleParticle(dt, did, rid, sid, lay, re, se, val) {}

  };


  class SX3  {
  public:
    int ID;
    int layer;
    int uds;
    std::vector<SX3FrontHit> frontHits;
    std::vector<SX3BackHit> backHits;
      
    //derived quantitites
    //std::vector<SX3Particle> particles;

    SX3() {}
    SX3(unsigned short int channel, unsigned short int value);
    void AddHit(unsigned short int channel, unsigned short int value);
    void MakeParticles(std::vector<SingleParticle*> &single_parts);
    //void Print(std::ofstream &out);

  };

  class QQQ5Hit  {
  public:
    unsigned short int ID; //ring or sector
    float value;  //pedestal subtracted
    float cal;    //calibrated (basic)

    QQQ5Hit() {}
    QQQ5Hit(unsigned short int id, float val, float c) : ID(id), value(val), cal(c) {}

  };

  class QQQ5   {
  public:
    int ID;
    int layer; //0=dE, 1=E
    int uds; //-1 = upstream, +1 = downstream
    int nGoodFronts;
    int nGoodBacks;
    std::vector<QQQ5Hit> frontHits;
    std::vector<QQQ5Hit> backHits;

    //derived quantitites
    //std::vector<QQQ5Particle> particles; //matched front and back hits

    QQQ5() {}
    QQQ5(unsigned short int channel, unsigned short int value);
    void AddHit(unsigned short int channel, unsigned short int value);
    void MakeParticles(std::vector<SingleParticle*> &single_parts);

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
    Tracker tracker; //! prevent this from being written to TTree

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

    bool SpuriousMyRIAD;     //! don't need this written to the tree

    Event() { sx3s.reserve(24); qqq5s.reserve(8); chans.reserve(4096); vals.reserve(4096); dither = 0.0; }
    void Set();

    void AddQQQ5(unsigned short int channel, unsigned short int value);
    void AddSX3(unsigned short int channel, unsigned short int value);
    void Reset() { sx3s.clear(); qqq5s.clear(); tracker.Reset();
      for (int i=0; i<single_parts.size(); ++i) {
        delete single_parts[i];
      }
      single_parts.clear();
      chans.clear();
      vals.clear();
      nhits = 0;
    }

    void Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length);
    void ProcessFinal() {};
    
    void SetConf(Configuration c) { conf = c; }

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
