#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include "ORRUBA/ORRUBAAuxiliary.hh"

namespace Orruba {
  void Auxiliary::Process(unsigned long long int ts,
                          unsigned short int *data,
                          unsigned short int length) {
    timestamp = ts;
    nRows = 0;
    
    int charCt = 0;
    for (int i=0; i<length; i=i+80) {
      nRows += 1;
      std::memcpy(&rows[nRows-1][0], &data[80*(nRows-1)], 80);
      rows[nRows-1][80]='\0';
    }
  }

  void Auxiliary::Print() {
    if (nRows == 0) { return; }
    std::cout << std::endl;
    for (int i=0; i<nRows; ++i) {
      if((int)rows[i][0] == 0) { continue; }
      std::cout << rows[i] << std::endl;
    }
  }
  
}
