#ifndef LIBS800PHYSICS_S800_HH
#define LIBS800PHYSICS_S800_HH

#include <iostream>
#include <fstream>
#include <functional>

#include "Reader/Reader.hh"

#include "S800/RawData.hh"
#include "S800/S800.hh"
#include "S800/S800Confs.hh"

namespace S800 {
  class PhysicsEvent : public DATOR::Processor {
  public:
    unsigned long long int timestamp;
    float ata, bta, dta, yta;
    float pTheta, pPhi;
    bool valid;
    
    void Process(unsigned long long int ts, unsigned short int *data, unsigned short int length);
    void ProcessFinal();
    void Reset();
  };

}

#endif
