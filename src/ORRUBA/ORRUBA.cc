#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include "ORRUBA/ORRUBA.hh"

namespace Orruba {
  float Event::dither;
  Configuration Event::conf;

  void Event::Process(unsigned long long int ts,
      unsigned short int *data,
      unsigned short int length) {

    timestamp = ts;

    unsigned int w1 = (timestamp & 0xFFFF00000000) >> 32;
    unsigned int w2 = (timestamp & 0xFFFF0000) >> 16;
    unsigned int w3 = (timestamp & 0xFFFF);

    if (w1 == w2 && w1 == w3) {
      SpuriousMyRIAD = true;
      nSpuriousMyRIAD += 1;
    }    

    nhits = length/4;
    chans.clear();
    vals.clear();
    for (int i=0; i<nhits; ++i) {
      if (data[2*i] == 0xFFFF && data[2*i+1] == 0xFFFF) { nhits = i; break; }
      chans.push_back(data[2*i] & 0x7FFF);
      vals.push_back(data[2*i+1]);
    }
    Set();
  }

  void Basic::Process(unsigned long long int ts,
      unsigned short int *data,
      unsigned short int length) {
    timestamp = ts;
    nhits = length/4;
    chans.clear();
    vals.clear();
    for (int i=0; i<nhits; ++i) {
      if (data[2*i] == 0xFFFF && data[2*i+1] == 0xFFFF) { nhits = i; break; }
      chans.push_back(data[2*i] & 0x7FFF);
      vals.push_back(data[2*i+1]);
    }
  }

  Configuration::~Configuration() { 
  }

  void Configuration::ReadConfiguration(std::string filename) {
    std::ifstream file(filename.c_str());
    std::string line;
    while (std::getline(file, line)) {
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      std::stringstream ss(line);

      std::string type;
      int chan, ID, lyr, fb, sub, lr, uds;
      float ped, off, gn;
      // should allow trailing comments or other things to not cause a segfault
      std::string trailing;
      ss >> chan >> type >> ID >> lyr >> fb >> uds >> sub >> lr >> ped >> off >> gn >> trailing;
      if (!type.compare("QQQ5")) { types[chan-1] = DetType::QQQ5; }
      else if (!type.compare("SX3")) { types[chan-1] = DetType::SX3; }
      else if (!type.compare("BB10")) { types[chan-1] = DetType::BB10; }
      else if (!type.compare("Track")) { types[chan-1] = DetType::Track; }
      else if (!type.compare("TDC")) { types[chan-1] = DetType::TDC; }
      else { std::cout << "Warning! Unrecognized detector type for channel " << chan << std::endl; }
      detID[chan-1]=ID;
      layer[chan-1]=lyr;
      side[chan-1]=fb;
      LR[chan-1] = lr;
      if (uds == 0) { //upstream
        UDS[chan-1] = -1;
      } else if (uds == 1) { //downstream
        UDS[chan-1] = +1;
      }
      subID[chan-1]=sub;
      pedestal[chan-1]=ped;
      offset[chan-1]=off;
      gain[chan-1]=gn;
      threshold[chan-1] = 0;
    }    
  }

  void Configuration::ReadCalibration(std::string filename) {
    std::ifstream file(filename.c_str());
    std::string line;
    while (std::getline(file, line)) {
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      std::stringstream ss(line);

      std::string type;
      ss >> type;

      int det, pad, strip, ring, sector;
      float padoff, padgain, stripLoff, stripLgain, stripRoff, stripRgain;
      float ringoff, ringgain, secoff, secgain;
      if (!type.compare("SX3")) {
        ss >> det >> pad >> strip >> padoff >> padgain >> stripLoff >> stripLgain >> stripRoff >> stripRgain;
        sx3cals[det-1].padoff[pad][strip] = padoff;
        sx3cals[det-1].padgain[pad][strip] = padgain;
        sx3cals[det-1].stripoff[0][pad][strip] = stripLoff;
        sx3cals[det-1].stripgain[0][pad][strip] = stripLgain;
        sx3cals[det-1].stripoff[1][pad][strip] = stripRoff;
        sx3cals[det-1].stripgain[1][pad][strip] = stripRgain;

      }
      else if (!type.compare("QQQ5")) {
        ss >> det >> ring >> sector >> ringoff >> ringgain >> secoff >> secgain;
        qqq5cals[det-1].ringoff[sector][ring] = ringoff;
        qqq5cals[det-1].ringgain[sector][ring] = ringgain;
        qqq5cals[det-1].secoff[sector][ring] = secoff;
        qqq5cals[det-1].secgain[sector][ring] = secgain;
      }      
    }
  }

  void Configuration::ReadRadii(std::string filename) {
    std::ifstream file(filename.c_str());
    std::string line;
    while (std::getline(file, line)) {
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      std::stringstream ss(line);

      int ringID;
      float avRadius;
      while (ss >> ringID >> avRadius) {
        radii[ringID] = avRadius;

        dEthetas[ringID] = std::atan(avRadius/dEdist)*180./M_PI;
        Ethetas[ringID] = std::atan(avRadius/Edist)*180./M_PI;

        if (dEthetas[ringID] < 0) { dEthetas[ringID] += 180.;}
        if (Ethetas[ringID] < 0) { Ethetas[ringID] += 180.;}
      }
    }
  }

  void Configuration::ReadPositionCalibration(std::string filename) {
    std::ifstream file(filename.c_str());
    std::string line;
    while (std::getline(file, line)) {
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      std::stringstream ss(line);

      int detID, stripID;
      double left, right;

      ss >> detID >> stripID >> left >> right;
      sx3cals[detID-1].stripPosGain[stripID] = 1./(right-left) * 75.;
      sx3cals[detID-1].stripPosOffset[stripID] = left;

      //pos = gain*(pos - left)
    }

    //calculate barrel dimensions
    r_dE_strip[0] = std::sqrt(std::pow(r_dE_barrel,2) + std::pow(3.*sx3width/8.,2));
    r_dE_strip[1] = std::sqrt(std::pow(r_dE_barrel,2) + std::pow(sx3width/8.,2));
    r_dE_strip[2] = std::sqrt(std::pow(r_dE_barrel,2) + std::pow(sx3width/8.,2));
    r_dE_strip[3] = std::sqrt(std::pow(r_dE_barrel,2) + std::pow(3.*sx3width/8.,2));

    r_dE_strip_bb10[0] = std::sqrt(std::pow(r_dE_barrel,2) + std::pow(7.*sx3width/16.,2));
    r_dE_strip_bb10[1] = std::sqrt(std::pow(r_dE_barrel,2) + std::pow(5.*sx3width/16.,2));
    r_dE_strip_bb10[2] = std::sqrt(std::pow(r_dE_barrel,2) + std::pow(3.*sx3width/16.,2));
    r_dE_strip_bb10[3] = std::sqrt(std::pow(r_dE_barrel,2) + std::pow(1.*sx3width/16.,2));
    r_dE_strip_bb10[4] = std::sqrt(std::pow(r_dE_barrel,2) + std::pow(1.*sx3width/16.,2));
    r_dE_strip_bb10[5] = std::sqrt(std::pow(r_dE_barrel,2) + std::pow(3.*sx3width/16.,2));
    r_dE_strip_bb10[6] = std::sqrt(std::pow(r_dE_barrel,2) + std::pow(5.*sx3width/16.,2));
    r_dE_strip_bb10[7] = std::sqrt(std::pow(r_dE_barrel,2) + std::pow(7.*sx3width/16.,2));

    r_E_strip[0] = std::sqrt(std::pow(r_E_barrel,2) + std::pow(3.*sx3width/8.,2));
    r_E_strip[1] = std::sqrt(std::pow(r_E_barrel,2) + std::pow(sx3width/8.,2));
    r_E_strip[2] = std::sqrt(std::pow(r_E_barrel,2) + std::pow(sx3width/8.,2));
    r_E_strip[3] = std::sqrt(std::pow(r_E_barrel,2) + std::pow(3.*sx3width/8.,2));
  }


  Configuration::Configuration() 
    : types(2048), detID(2048), layer(2048), side(2048), subID(2048), LR(2048), UDS(2048), pedestal(2048), offset(2048), gain(2048), threshold(2048), sx3cals(48), qqq5cals(8) {
    }

  Configuration::Configuration(std::string filename) : Configuration() {
    ReadConfiguration(filename);
  }

  Configuration::Configuration(std::string name, std::string title, std::string filename) : Configuration(filename) {
    //SetName(name.c_str());
    //SetTitle(title.c_str());
  }

  void Configuration::Set(std::string name, std::string title, std::string filename) {
    ReadConfiguration(filename);
  }

  void Configuration::SetThresholds(int chan, float thresh) {
    threshold[chan-1] = thresh;
  }

  void Configuration::SetThresholds(float thresh) {
    for (int i=0; i<threshold.size(); ++i) {
      threshold[i] = thresh;
    }
  }

  void Configuration::SetThresholds(DetType type, float thresh) {
    for (int i=0; i<threshold.size(); ++i) {
      if (types[i] == type) {
        threshold[i] = thresh;
      }
    }
  }

  void Configuration::SetThresholds(DetType type, int s, float thresh) {
    for (int i=0; i<threshold.size(); ++i) {
      if (types[i] == type && side[i] == s) {
        threshold[i] = thresh;
      }
    }
  }

  //  void Event::Set(unsigned short int *channels,
  //                  unsigned short int *values,
  //                  int &nhits) {
  void Event::Set() {
    bool sx3evt = false;
    bool qqq5evt = false;
    bool bb10evt = false;
    int qqq5hitsth = 0;
    for (int i=0; i<nhits; ++i) {
      unsigned short int chan = chans[i];
      unsigned short int val = vals[i];
      if (conf.types[chan-1] == DetType::QQQ5) {
        qqq5hits += 1;
        int retval = AddQQQ5(chan, val);
        if (retval) { qqq5hitsth+= 1; nQQQ5HitsTh += 1; }
        if (!qqq5evt && retval)  { qqq5evt = true; nQQQ5Evts += 1; }
      }
      else if (conf.types[chan-1] == DetType::SX3) {
        sx3hits += 1;
        int retval = AddSX3(chan, val);
        if (retval) { nSX3HitsTh += 1; }
        if (!sx3evt && retval) { sx3evt = true; nSX3Evts += 1; }
      }
      else if (conf.types[chan-1] == DetType::BB10) {
        bb10hits += 1;
        int retval = AddBB10(chan, val);
        if (retval) { nBB10HitsTh += 1; }
        if (!bb10evt && retval) { bb10evt = true; nBB10Evts += 1; }
      }
      else if (conf.types[chan-1] == DetType::Track) {
        tracker.Set(chan, val);
      }
      else if (conf.types[chan-1] == DetType::TDC) {
        tdc.SetChan(chan, val);
      }
    }

    for (int i=0; i<sx3s.size(); ++i) {
      sx3parts += sx3s[i].MakeParticles(single_parts);
    }
    nSX3Particles += sx3parts;

    for (int i=0; i<qqq5s.size(); ++i) {
      qqq5parts += qqq5s[i].MakeParticles(single_parts);
    }
    nQQQ5Particles += qqq5parts;
    if (qqq5parts*2 > qqq5hitsth) { std::cout << "In this event there are " << qqq5parts << " reconstructed QQQ5 particles but only " << qqq5hitsth << " hits above threshold" << std::endl; }

    for (int i=0; i<bb10s.size(); ++i) {
      bb10parts += bb10s[i].MakeParticles(single_parts);
    }
    nBB10Particles += bb10parts;

    if (sx3evt && sx3parts == 0) {
      nBadSX3Evts += 1;
    }
    if (qqq5evt && qqq5parts == 0) {
      nBadQQQ5Evts += 1;
    }
    if (bb10evt && bb10parts == 0) {
      nBadBB10Evts += 1;
    }

    tracker.presentX = 0;
    tracker.presentY = 0;      
    for (int i=0; i<16; ++i) {        
      tracker.presentX += tracker.xhits[i].Validate();
      tracker.presentY += tracker.yhits[i].Validate();
    }
    tracker.SetX();
    tracker.SetY();
    //tracker.Print(std::cout);
  }

  int Event::AddQQQ5(unsigned short int channel, unsigned short int value) {
    int retval = 0;
    for (int i=0; i<qqq5s.size(); ++i) {
      if (conf.detID[channel-1] == qqq5s[i].ID) {
        retval = qqq5s[i].AddHit(channel, value);
        nQQQ5Hits += 1;
        return retval;
      }
    }
    if (((float)value-conf.pedestal[channel-1]) > conf.threshold[channel-1]) {
      qqq5s.emplace_back(channel, value);
      retval = 1;
    }
    nQQQ5Hits += 1;
    return retval;
  }

  int Event::AddSX3(unsigned short int channel, unsigned short int value) {
    int retval = 0;
    for (int i=0; i<sx3s.size(); ++i) {
      if (conf.detID[channel-1] == sx3s[i].ID) {
        retval = sx3s[i].AddHit(channel, value);
        nSX3Hits += 1;
        return retval;
      }
    }
    if (((float)value-conf.pedestal[channel-1]) > conf.threshold[channel-1]) {
      sx3s.emplace_back(channel, value);
      retval = 1;
    }
    nSX3Hits += 1;
    return retval;
  }

  int Event::AddBB10(unsigned short int channel, unsigned short int value) {
    int retval = 0;
    for (int i=0; i<bb10s.size(); ++i) {
      if (conf.detID[channel-1] == bb10s[i].ID) {
        retval = bb10s[i].AddHit(channel, value);
        nBB10Hits += 1;
        return retval;
      }
    }
    if (((float)value - conf.pedestal[channel-1]) > conf.threshold[channel-1]) {
      bb10s.emplace_back(channel, value);
      retval = 1;
    }
    nBB10Hits += 1;
    return retval;
  }

  void Event::PrintSummary(std::ostream &out) {
    out << "--------- ORRUBA Summary --------" << std::endl;
    if (nQQQ5Hits>0) {
      out << "   QQQ5 " << ANSI_COLOR_YELLOW << nQQQ5Hits << ANSI_COLOR_RESET << " -> " << ANSI_COLOR_YELLOW << nQQQ5HitsTh << ANSI_COLOR_RESET << " -> " << ANSI_COLOR_GREEN << nQQQ5Particles << ANSI_COLOR_RESET << " reconstructed particles" << std::endl;
      out << "         " << ANSI_COLOR_GREEN << (double)nQQQ5HitsTh/(double)nQQQ5Hits*100.0 << "%" << ANSI_COLOR_RESET << " hits above threshold" << std::endl;
      out << "         " << ANSI_COLOR_GREEN << (double)nQQQ5HitsTh/(double)nQQQ5Particles << ANSI_COLOR_RESET << " hits above threshold/particle" << std::endl;
      out << "         " << ANSI_COLOR_RED << nBadQQQ5Evts << ANSI_COLOR_RESET << " events with QQQ5 hits but no reconstructed particle (" << ANSI_COLOR_RED << (double)nBadQQQ5Evts/(double)nQQQ5Evts*100.0 << "%" << ANSI_COLOR_RESET << ")" << std::endl;
    }
    if (nSX3Hits>0) {
      out << "   SX3 " << ANSI_COLOR_YELLOW << nSX3Hits << ANSI_COLOR_RESET << " -> " << ANSI_COLOR_YELLOW << nSX3HitsTh << ANSI_COLOR_RESET << " -> " << ANSI_COLOR_GREEN << nSX3Particles << ANSI_COLOR_RESET << " reconstructed particles" << std::endl;
      out << "         " << ANSI_COLOR_GREEN << (double)nSX3HitsTh/(double)nSX3Hits*100.0 << "%" << ANSI_COLOR_RESET << " hits above threshold" << std::endl;
      out << "         " << ANSI_COLOR_GREEN << (double)nSX3HitsTh/(double)nSX3Particles << ANSI_COLOR_RESET << " hits above threshold/particle" << std::endl;
      out << "         " << ANSI_COLOR_RED << nBadSX3Evts << ANSI_COLOR_RESET << " events with SX3 hits but no reconstructed particle (" << ANSI_COLOR_RED << (double)nBadSX3Evts/(double)nSX3Evts*100.0 << "%" << ANSI_COLOR_RESET << ")" << std::endl;
    }
    if (nBB10Hits>0) {
      out << "   BB10 " << ANSI_COLOR_YELLOW << nBB10Hits << ANSI_COLOR_RESET << " -> " << ANSI_COLOR_YELLOW << nBB10HitsTh << ANSI_COLOR_RESET << " -> " << ANSI_COLOR_GREEN << nBB10Particles << ANSI_COLOR_RESET << " reconstructed particles" << std::endl;
      out << "         " << ANSI_COLOR_GREEN << (double)nBB10HitsTh/(double)nBB10Hits*100.0 << "%" << ANSI_COLOR_RESET << " hits above threshold" << std::endl;
      out << "         " << ANSI_COLOR_GREEN << (double)nBB10HitsTh/(double)nBB10Particles << ANSI_COLOR_RESET << " hits above threshold/particle" << std::endl;
      out << "         " << ANSI_COLOR_RED << nBadBB10Evts << ANSI_COLOR_RESET << " events with BB10 hits but no reconstructed particle (" << ANSI_COLOR_RED << (double)nBadBB10Evts/(double)nBB10Evts*100.0 << "%" << ANSI_COLOR_RESET << ")" << std::endl;
    }
    if (nSpuriousMyRIAD>0) {
      out << "   " << ANSI_COLOR_RED << nSpuriousMyRIAD << ANSI_COLOR_RESET << " events " << std::endl;
    }
  }

  void SingleParticle::OffsetBeam(float beamx, float beamy, bool verbose) {

    if (verbose) {
      std::cout << "============= " << beamx << ", " << beamy << " ===========" << std::endl;
      std::cout << "=========== OLD ===========" << std::endl;
      std::cout << x << "   " << y << "  " << z << std::endl;
      std::cout << r << "   " << phi << std::endl;
      std::cout << theta << std::endl;
    }

    x_off = x - beamx;
    y_off = y - beamy;
    z_off = z;

    r_off = std::sqrt(x_off*x_off + y_off*y_off);
    theta_off = std::atan( r_off / z_off );
    if (theta_off < 0) { theta_off += M_PI; }
    if (x_off >= 0) {
      phi_off = std::asin(y_off / r_off);
    }
    else if (x_off < 0) {
      phi_off = -std::asin(y_off / r_off) + M_PI;
    }
    if (phi_off < 0) { phi_off += 2.0*M_PI; }
    //phi = std::acos(x / r);
    /*
       if (phi != std::asin(y / r)) {
       std::cout << "Warning! Phi does not match" << std::endl;
       std::cout << phi << "   " << std::asin(y/r) << std::endl;
       }
       */

    if (verbose) {
      std::cout << "=========== NEW ===========" << std::endl;
      std::cout << x_off << "   " << y_off << "  " << z_off << std::endl;
      std::cout << r_off << "   " << phi_off << std::endl;
      std::cout << theta_off << std::endl;
    }
  }

  int TrackerHit::Validate() {
    int present = 0;
    if (value != 0 || tdc != 0) { present = 1; }
    if (value == 0) { valid = false; return present; }
    if (tdc == 0) { valid = false; return present; }

    if (tdc > 1000 && tdc < 2000) { valid = true; return present; }
    valid = false;
    return present;

  }

  void Tracker::Set(unsigned short channel, unsigned short value) {
    if (Event::conf.side[channel-1] == 0) {
      cathode_val = value - Event::conf.pedestal[channel-1] + Orruba::Event::Dither();
      cathode = cathode_val * Event::conf.gain[channel-1] + Event::conf.offset[channel-1];
      return;
    }
    else if (Event::conf.side[channel-1] == 1) { //energies
      if (Event::conf.layer[channel-1] == 0) {
        float val = value - Event::conf.pedestal[channel-1] + Orruba::Event::Dither();
        float cal = val * Event::conf.gain[channel-1] + Event::conf.offset[channel-1];
        xhits[Event::conf.subID[channel-1]].value =  value;
        xhits[Event::conf.subID[channel-1]].cal =  cal;
        xhits[Event::conf.subID[channel-1]].valid = true;
        return;
      }
      else if (Event::conf.layer[channel-1] == 1) {
        float val = value - Event::conf.pedestal[channel-1] + Orruba::Event::Dither();
        float cal = val * Event::conf.gain[channel-1] + Event::conf.offset[channel-1];
        yhits[Event::conf.subID[channel-1]].value =  value;
        yhits[Event::conf.subID[channel-1]].cal =  cal;
        yhits[Event::conf.subID[channel-1]].valid = true;
        return;
      }
    }
    else if (Event::conf.side[channel-1] == 2) { //times
      if (Event::conf.layer[channel-1] == 0) {  //x
        float val = value - Event::conf.pedestal[channel-1] + Orruba::Event::Dither();
        xhits[Event::conf.subID[channel-1]].tdc =  val;
        return;
      }
      else if (Event::conf.layer[channel-1] == 1) { //y
        float val = value - Event::conf.pedestal[channel-1] + Orruba::Event::Dither();
        yhits[Event::conf.subID[channel-1]].tdc =  val;
        return;
      }
    }
  }    

  void Tracker::SetX() {
    float maxval;
    int maxx = -1;
    float wsum = 0;
    float sum = 0;
    nGoodX = 0;
    float goodEns[2];
    float goodWires[2];
    validX = false;
    for (int i=0; i<16; ++i) {
      if (xhits[i].valid == false) { continue; }
      if (xhits[i].cal > 20) {
        wsum += i*xhits[i].cal;
        sum += xhits[i].cal;
        if (nGoodX<2) {
          goodEns[nGoodX] = xhits[i].cal;
          goodWires[nGoodX] = i;
        }
        nGoodX += 1;
      }
      if (xhits[i].cal > maxval) {
        maxval = xhits[i].cal;
        maxx = i;
      }
    }
    x = maxx;
    if (nGoodX == 1) { validX = true; x = goodWires[0]; }
    else if (nGoodX == 2) {
      if (abs(goodWires[0] - goodWires[1]) == 1) {
        validX = true;
        x = (goodWires[0] + goodWires[1])/2.0;
      }
    }

    //x = wsum/sum;
    //if (x > 0 && x < 16) { validX = true; }
    ycal = -(x*2.0 - 8.5*2.0); //now into gretina frame
  }

  void Tracker::SetY() {
    float maxval;
    int maxy = -1;
    float wsum = 0;
    float sum = 0;
    nGoodY = 0;
    validY = false;
    float goodEns[2];
    float goodWires[2];
    for (int i=0; i<16; ++i) {
      if (yhits[i].valid == false) { continue; }
      if (yhits[i].cal > 20) {
        wsum += i*yhits[i].cal;
        sum += yhits[i].cal;
        if (nGoodY<2) {
          goodEns[nGoodY] = yhits[i].cal;
          goodWires[nGoodY] = i;
        }
        nGoodY += 1;
      }
      if (yhits[i].cal > maxval) {
        maxval = yhits[i].cal;
        maxy = i;
      }
    }
    y = maxy;

    if (nGoodY == 1) { validY = true; y = goodWires[0]; }
    else if (nGoodY == 2) {
      if (abs(goodWires[0] - goodWires[1]) == 1) {
        validY = true;
        y = (goodWires[0] + goodWires[1])/2.0;
      }
    }

    //y = wsum/sum ;
    //if (y > 0 && y < 16) { validY = true; }
    xcal = -(y*2.0 - 8.5*2.0); //now into gretina frame
  }

  void Tracker::Print(std::ostream &out) {
    out << "=================" << std::endl;
    for (int i=0; i<16; ++i) {
      out << "xh: " << i  << "   " << xhits[i].cal << std::endl;
    }
    for (int i=0; i<16; ++i) {
      out << "yh: " << i << "   " << yhits[i].cal << std::endl;
    }
  }
}
