#ifndef LIBS800_S800CONFS_HH
#define LIBS800_S800CONFS_HH

#include <cmath>
#include <cstring>
#include <sstream>
#include <vector>
#include <functional>

#include "S800/RawData.hh"

namespace S800 {
  class ConfFile { //prototype to handle blocks of tables
  public:
    std::vector< std::function <void(std::stringstream &ss)> > readers;
    std::vector<std::string> headers;

    FILE *file;

    ConfFile() {}
    ConfFile(std::string fn) { Open(fn); }
    int Open(std::string fn);
    void Read(std::string fn);
    void Read();
  };

  class ICConf : public ConfFile {
  public:
    double gain[16];
    double ROI[2];

    ICConf();

    void ReadIC(std::stringstream &ss);
    void ReadROI(std::stringstream &ss);
  };

  class CRDCConf : public ConfFile {
  public:
    double gainx[2];
    double gainy[2];
    double offsetx[2];
    double offsety[2];      

    int pedestals[2][224];
    double gains[2][224];

    void ReadCal(std::stringstream &ss);
    void ReadPeds(std::stringstream &ss);
    void ReadGains(std::stringstream &ss);

    CRDCConf();
    //int Read(std::string fn);
  };

  class MTDCConf : public ConfFile {
  public:
      
    double low[S800_MTDC_MAXGATES][S800_MTDC_MAXCHANS]; //index is based on trigger condition
    double high[S800_MTDC_MAXGATES][S800_MTDC_MAXCHANS];

    //double obj_drifts[100];
    //double xfp_drifts[100];
      
    double e1low = -65536;
    double e1high = 65536;

    //double xfp_offset = 0;

    MTDCConf();

    void ReadRawGates(std::stringstream &ss);
    void ReadE1Gate(std::stringstream &ss);      
  };


}

#endif
