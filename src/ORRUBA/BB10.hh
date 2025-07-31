#ifndef LIBORR_BB10_HH
#define LIBORR_BB10_HH

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

#include "Reader/Reader.hh"
#include "ORRUBA/SingleParticle.hh"

namespace Orruba {
  
  class bb10cal  {
  public:
    float stripoff[8];
    float stripgain[8];
    
    bb10cal() {
      for (int i=0; i<8; ++i) {
        stripoff[i] = 0.0;
        stripgain[i] = 1.0;    
      }
    }

  };

  class BB10Hit  {
  public:
    unsigned short int ID; //strip
    float value;  //pedestal subtracted
    float cal;    //calibrated (basic)

    BB10Hit() {}
    BB10Hit(unsigned short int id, float val, float c) : ID(id), value(val), cal(c) {}

  };

  class BB10   {
  public:
    int ID;
    int layer; //0=dE, 1=E
    int uds; //-1 = upstream, +1 = downstream
    std::vector<BB10Hit> frontHits;

    BB10() {}
    BB10(unsigned short int channel, unsigned short int value);
    int AddHit(unsigned short int channel, unsigned short int value);
    int MakeParticles(std::vector<SingleParticle*> &single_parts);

  };

  class BB10Particle : public SingleParticle {
  public:
      
    void MakeCoords(BB10 *detector);
    BB10Particle() {}
    BB10Particle(DetType dt, unsigned short int did,
                 unsigned short int sid,
                 unsigned short int lay,
                 float se, bool val);

  };

}


#endif
