#ifndef LIBORR_ORRUBA_AUX_HH
#define LIBORR_ORRUBA_AUX_HH

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

#include <vector>
#include <string>

#include "Reader/Reader.hh"

namespace Orruba {
  
  class Auxiliary : public DATOR::Processor {
  public:
    unsigned long long timestamp;
    int nRows;
    char rows[80][81];
    
    Auxiliary() { nRows = 0; }

    void Reset() { nRows = 0; }
    void Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length);
    void ProcessFinal() {};

    void Print();
  };
}  

#endif
