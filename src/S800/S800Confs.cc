#include <cmath>
#include <cstring>
#include <sstream>

#include "S800/S800Confs.hh"

namespace S800 {
  int ConfFile::Open(std::string fn) {
    if ( file = fopen(fn.c_str(), "ra") ) {
      return 1;
    }
    std::cerr << "Warning! Configuration file " << fn << " not found!" << std::endl;
    return 0;
  }

  void ConfFile::Read(std::string fn) {
    if (Open(fn)) {
      Read();
    }
  }

  void ConfFile::Read() {
    std::stringstream ss;
    char cline[2048];

    while(std::fgets(cline, sizeof cline, file)!=NULL) {    
      std::string line(cline);
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      ss.clear();
      ss.str(line);

      std::string header;
      ss >> header;

      for (int i=0; i<headers.size(); ++i) {
        if (!header.compare(headers[i])) {
          while(std::fgets(cline, sizeof cline, file)!=NULL) {    
            std::string line(cline);
            if (line.size() == 0) { continue; }
            if (line[0] == '#') { continue; }
            if (line[0] == ';') { continue; }

            ss.clear();
            ss.str(line);

            ss >> header;
            if (!header.compare("END")) { break; }

            ss.clear();
            ss.str(line);

            readers[i](ss);
          }
        }
      }
    }
  }

  ICConf::ICConf() : ConfFile() {
    for (int i=0; i<16; ++i) { gain[i] = 1.0; }
    headers.push_back("CALIBRATION");
    readers.push_back([this](std::stringstream &ss) { ReadIC(ss); });
    headers.push_back("ROI");
    readers.push_back([this](std::stringstream &ss) { ReadROI(ss); });        
  }

  void ICConf::ReadIC(std::stringstream &ss) {
    int seg;
    double g;
    ss >> seg >> g;
    gain[seg] = g;
  }
  void ICConf::ReadROI(std::stringstream &ss) {
    double low, high;
    ss >> low >> high;
    ROI[0] = low;
    ROI[1] = high;
  }

  CRDCConf::CRDCConf() {
    for (int i=0; i<2; ++i) {
      gainx[i] = 1.0; gainy[i] = 1.0; offsetx[i] = 0.0; offsety[i] = 0.0;
      for (int j=0; j<224; ++j) {
        pedestals[i][j] = 0;
        gains[i][j] = 1.0;
      }
    }

    readers.push_back([this](std::stringstream &ss) { ReadCal(ss); });
    headers.push_back("CALIBRATION");

    readers.push_back([this](std::stringstream &ss) { ReadPeds(ss); });
    headers.push_back("PEDESTALS");

    readers.push_back([this](std::stringstream &ss) { ReadGains(ss); });
    headers.push_back("GAINS");
  }
  
  void CRDCConf::ReadCal(std::stringstream &ss) {
    int indx;
    double gx,ox;
    double gy,oy;
    ss >> indx >> ox >> gx >> oy >> gy;
    gainx[indx] = gx;
    gainy[indx] = gy;
    offsetx[indx] = ox;
    offsety[indx] = oy;
  }

  void CRDCConf::ReadPeds(std::stringstream &ss) {
    int indx;
    int ped0;
    int ped1;
    ss >> indx >> ped0 >> ped1;
    pedestals[0][indx] = ped0;
    pedestals[1][indx] = ped1;
  }

  void CRDCConf::ReadGains(std::stringstream &ss) {
    int indx;
    double gain0, gain1;
    ss >> indx >> gain0 >> gain1;
    gains[0][indx] = gain0;
    gains[1][indx] = gain1;
  }

    
  MTDCConf::MTDCConf() {
    for (int i=0; i<S800_MTDC_MAXCHANS; ++i) {
      for (int j=0; j<S800_MTDC_MAXGATES; ++j) {
        low[j][i] = 0;
        high[j][i] = 65536;
      }
    }

    headers.push_back("RAWGATES");
    readers.push_back([this](std::stringstream &ss) { ReadRawGates(ss); } );

    headers.push_back("E1GATE");
    readers.push_back([this](std::stringstream &ss) { ReadE1Gate(ss); } );
  }

  void MTDCConf::ReadRawGates(std::stringstream &ss) {
    int trig, chan;
    double l1,l2;
    double h1,h2;
    ss >> trig >> chan >> l1 >> h1;
          
    low[trig][chan] = l1;
    high[trig][chan] = h1;
  }

  void MTDCConf::ReadE1Gate(std::stringstream &ss) {
    ss >> e1low >> e1high;
  }    
    
}

