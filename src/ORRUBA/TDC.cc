#include "ORRUBA/ORRUBA.hh"
#include "ORRUBA/TDC.hh"

namespace Orruba {

  TDC::TDC() : TDC(32) {}

  TDC::TDC(int nChans) {
    nChannels = nChans;
    for (int i=0; i<nChannels; ++i) {
      raw_vals.push_back(0);
    }
  }

  void TDC::Reset() {
    std::memset(raw_vals.data(), 0, nChannels);  
  }

  void TDC::SetChan(unsigned short int channel, unsigned short int val) {
    raw_vals[Event::conf.subID[channel-1]] = val;
  }
}

