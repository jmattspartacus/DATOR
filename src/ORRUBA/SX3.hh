#ifndef LIBORR_SX3_HH
#define LIBORR_SX3_HH

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

#include "Reader/Reader.hh"
#include "ORRUBA/SingleParticle.hh"

namespace Orruba {
  
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
    int AddHit(unsigned short int channel, unsigned short int value);
    int MakeParticles(std::vector<SingleParticle*> &single_parts);
    //void Print(std::ofstream &out);

  };

  
  class SX3Particle : public SingleParticle {
  public:
    //some extra SX3 only things
    float stripL;
    float stripR;
    float position; //-1 -> 1;

    SX3Particle() {}
    SX3Particle(DetType dt, unsigned short int did,
                unsigned short int pid, unsigned short int sid,
                unsigned short int lay,
                float pe, float sl, float sr, float pos, bool val);
    void MakeCoords(SX3 *detector);


  };
}

#endif
