#ifndef LIBORR_QQQ5_HH
#define LIBORR_QQQ5_HH

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

#include "Reader/Reader.hh"
#include "ORRUBA/SingleParticle.hh"

namespace Orruba {

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
    
    ~qqq5cal(){ }
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
    int AddHit(unsigned short int channel, unsigned short int value);
    int MakeParticles(std::vector<SingleParticle*> &single_parts);

  };


  class QQQ5Particle : public SingleParticle {
  public:
      
    void MakeCoords(QQQ5 *detector);
    QQQ5Particle() {}
    QQQ5Particle(DetType dt, unsigned short int did,
                 unsigned short int rid, unsigned short int sid,
                 unsigned short int lay,
                 float re, float se, bool val);

  };
}

#endif
