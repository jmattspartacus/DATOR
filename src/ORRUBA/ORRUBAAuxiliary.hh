#ifndef LIBORR_ORRUBA_AUX_HH
#define LIBORR_ORRUBA_AUX_HH

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <ctime>
#include <iomanip>

#include <vector>
#include <string>

#include "Reader/Reader.hh"

namespace Orruba {
  class Table {
  public:
    int ncols;
    int nrows;

    double *data;
    std::vector<std::string> names;
    
    Table();
    Table(int nr, int nc);
    ~Table();

    void Initialize(int nr, int nc);
    double Get(int row, int col);
    void Set(int row, int col, double val);
    void Print();    
  };
    
  
  class Auxiliary : public DATOR::Processor {
  public:
    unsigned long long timestamp;
    int nRows;
    char rows[80][81];

    std::string header;
    time_t time_since_epoch;
    Table scalers;
    
    Auxiliary() { nRows = 0; }

    void Reset() { nRows = 0; }
    void Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length);
    void ProcessFinal() {};

    void Print();
  };
}  

#endif
