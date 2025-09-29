#ifndef LIBORR_TDC_HH
#define LIBORR_TDC_HH

#include <vector>
#include <string>

namespace Orruba {
  class TDCcal {
  public: 
    std::vector<float> offsets;
    std::vector<float> gains;
  };

  class TDC {
  public:
    int nChannels;
    std::vector<int> raw_vals;

    TDC();
    TDC(int nChans);
    void SetChan(unsigned short int channel, unsigned short int val);
    void Reset();
  };
};

#endif
