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
      ss >> chan >> type >> ID >> lyr >> fb >> uds >> sub >> lr >> ped >> off >> gn;
      types[chan-1]=type;
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

    r_E_strip[0] = std::sqrt(std::pow(r_E_barrel,2) + std::pow(3.*sx3width/8.,2));
    r_E_strip[1] = std::sqrt(std::pow(r_E_barrel,2) + std::pow(sx3width/8.,2));
    r_E_strip[2] = std::sqrt(std::pow(r_E_barrel,2) + std::pow(sx3width/8.,2));
    r_E_strip[3] = std::sqrt(std::pow(r_E_barrel,2) + std::pow(3.*sx3width/8.,2));
  }
    

  Configuration::Configuration() 
    : types(2048), detID(2048), layer(2048), side(2048), subID(2048), LR(2048), UDS(2048), pedestal(2048), offset(2048), gain(2048), sx3cals(24), qqq5cals(8) {
  }
  
  Configuration::Configuration(std::string filename) : Configuration() {
    ReadConfiguration(filename);
  }

  Configuration::Configuration(std::string name, std::string title, std::string filename) : Configuration(filename) {
    //SetName(name.c_str());
    //SetTitle(title.c_str());
  }

  //  void Event::Set(unsigned short int *channels,
  //                  unsigned short int *values,
  //                  int &nhits) {
  void Event::Set() {
    for (int i=0; i<nhits; ++i) {
      unsigned short int chan = chans[i];
      unsigned short int val = vals[i];
      if (!conf.types[chan-1].compare("QQQ5")) {
        AddQQQ5(chan, val);
      }
      else if (!conf.types[chan-1].compare("SX3")) {
        AddSX3(chan, val);
      }
      else if (!conf.types[chan-1].compare("Track")) {
        tracker.Set(chan, val);
      }
    }

    for (int i=0; i<sx3s.size(); ++i) {
      sx3s[i].MakeParticles(single_parts);
    }

    for (int i=0; i<qqq5s.size(); ++i) {
      qqq5s[i].MakeParticles(single_parts);
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

  void Event::AddQQQ5(unsigned short int channel, unsigned short int value) {
    for (int i=0; i<qqq5s.size(); ++i) {
      if (conf.detID[channel-1] == qqq5s[i].ID) {
        qqq5s[i].AddHit(channel, value);
        return;
      }
    }
    qqq5s.emplace_back(channel, value);
  }

  void Event::AddSX3(unsigned short int channel, unsigned short int value) {
    for (int i=0; i<sx3s.size(); ++i) {
      if (conf.detID[channel-1] == sx3s[i].ID) {
        sx3s[i].AddHit(channel, value);
        return;
      }
    }
    sx3s.emplace_back(channel, value);
  }

  SX3::SX3(unsigned short int channel, unsigned short int value) {
    ID = Event::conf.detID[channel-1];
    layer = Event::conf.layer[channel-1];
    uds = Event::conf.UDS[channel-1];
    float val = (value-Event::conf.pedestal[channel-1] + Orruba::Event::Dither());
    float cal = val*Event::conf.gain[channel-1] + Event::conf.offset[channel-1];
    frontHits.reserve(8);
    backHits.reserve(4);
    if (Event::conf.side[channel-1] == 0) { //front
      frontHits.emplace_back(Event::conf.LR[channel-1], Event::conf.subID[channel-1], val, cal);
    }
    else if (Event::conf.side[channel-1] == 1) { //back
      backHits.emplace_back(Event::conf.subID[channel-1], val, cal);
    }
    else {
      std::cerr << "Warning! Channel " << channel << " with invalid front/back value" << std::endl;
    }
  }
  void SX3::AddHit(unsigned short int channel, unsigned short int value) {
    float val = (value-Event::conf.pedestal[channel-1] + Orruba::Event::Dither());
    float cal = val*Event::conf.gain[channel-1] + Event::conf.offset[channel-1];    
    if (Event::conf.side[channel-1] == 0) { //front
      frontHits.emplace_back(Event::conf.LR[channel-1], Event::conf.subID[channel-1], val, cal);
    }
    else if (Event::conf.side[channel-1] == 1) { //back
      backHits.emplace_back(Event::conf.subID[channel-1], val, cal);
    }
    else {
      std::cerr << "Warning! Channel " << channel << " with invalid front/back value" << std::endl;
    }
  }

  void SX3Particle::MakeCoords(SX3 *detector) {
    //this r is in the cylindrical coordinate system
    if (detector->layer == 0) { r = Event::conf.r_dE_strip[frontID]; }
    else if (detector->layer == 1) { r = Event::conf.r_E_strip[frontID]; }

    //theta in spherical polar
    theta = std::atan(r/position);
    if (theta < 0 ) { theta += 3.14159265358799323846; }      
      
    float deltaphi = std::atan(Event::conf.sx3width/8.*std::abs((int)(3-2*frontID))/Event::conf.r_dE_barrel);
    int sign = 1;
    if (detector->layer == 1) {
      if (frontID < 2) { sign = 1; } //low strip ID = high phi for E
      else if (frontID >= 2) { sign = -1; }
    }
    else if (detector->layer == 0) {
      if (frontID < 2) { sign = -1; } //other way around for dE layer
      else if (frontID >= 2) { sign = 1; }
    }

    //phi of SX3 centre, common between cylindrical + spherical systems
    phi = (detector->ID+6)%12 * 2.*3.14159265/12.; //0-> down in lab frame-> +x, 6->up
    phi += sign*deltaphi;
    if (phi<0) { phi += 2.*3.14159265; }

    z = position;
    x = r*std::cos(phi);
    y = r*std::sin(phi);     
  }

  void SX3::MakeParticles(std::vector<SingleParticle*> &single_parts) {
    int frontPairs = 0;
    double frontSum = 0;
    double pos = 0;
    int stripID = 0;
    float stripL = 0;
    float stripR = 0;

    for (int i=0; i<frontHits.size(); ++i) {
      int subID1 = frontHits[i].ID;
      int end1 = frontHits[i].side;
      for (int j=i+1; j<frontHits.size(); ++j) {
        int subID2 = frontHits[j].ID;
        int end2 = frontHits[j].side;
        if (subID1 == subID2) {
          if (end1 == 1 && end2 == 0) {
            frontPairs += 1;
            stripID = subID1;
            stripL = frontHits[j].value;
            stripR = frontHits[i].value;
            frontSum = frontHits[i].cal + frontHits[j].cal;
          }
          else if (end2 == 1 && end1 == 0) {
            frontPairs += 1;
            stripID = subID1;
            stripL = frontHits[i].value;
            stripR = frontHits[j].value;
            frontSum = frontHits[i].cal + frontHits[j].cal;
          }
        }
      }
    }
    int validBacks = 0;
    float padEnergy = 0;
    int padID = 0;
    for (int i=0; i<backHits.size(); ++i) {
      if (backHits[i].cal > 250.0) {
        validBacks += 1;
        padID = backHits[i].ID;
        padEnergy = backHits[i].value; 
      }
    }
    if (frontPairs == 1 && validBacks > 0) {
      //validBacks == 1) {
      //look through for closest back using naive energy calibration
      float closest = 99999;
      for (int i=0; i<backHits.size(); ++i) {
        if (abs(backHits[i].cal - frontSum) < closest) {
          closest = abs(backHits[i].cal - frontSum);
          padID = backHits[i].ID;
          padEnergy = backHits[i].value;
        }
      }
      
      //now we have matched, do calibrations
      float padgain = Event::conf.sx3cals[ID-1].padgain[padID][stripID];
      float padoff = Event::conf.sx3cals[ID-1].padoff[padID][stripID];

      float stripLgain = Event::conf.sx3cals[ID-1].stripgain[0][padID][stripID];
      float stripLoff = Event::conf.sx3cals[ID-1].stripoff[0][padID][stripID];

      float stripRgain = Event::conf.sx3cals[ID-1].stripgain[1][padID][stripID];
      float stripRoff = Event::conf.sx3cals[ID-1].stripoff[1][padID][stripID];

      padEnergy = padEnergy*padgain + padoff;
      stripL = stripL*stripLgain + stripLoff;
      stripR = stripR*stripRgain + stripRoff;

      frontSum = stripL + stripR;
      //check energy equivalence
      pos = (stripR - stripL)/frontSum;
      pos = (pos - Event::conf.sx3cals[ID-1].stripPosOffset[stripID])*Event::conf.sx3cals[ID-1].stripPosGain[stripID];
      
      bool valid = (padEnergy/frontSum  >= 0.95 && padEnergy/frontSum <= 1.05);

      SX3Particle *part = new SX3Particle(1, ID,
                                          padID, stripID,
                                          layer,
                                          padEnergy, stripL, stripR, pos*uds, valid);
      part->MakeCoords(this);
      single_parts.push_back(part);
    }
    else {   }
  }
    
    
  QQQ5::QQQ5(unsigned short int channel, unsigned short int value) {
    ID = Event::conf.detID[channel-1];
    layer = Event::conf.layer[channel-1];
    uds = Event::conf.UDS[channel-1];
    nGoodFronts = 0;
    nGoodBacks = 0;
    float val = (value-Event::conf.pedestal[channel-1] + Orruba::Event::Dither());
    float cal = val*Event::conf.gain[channel-1] + Event::conf.offset[channel-1];
    /*
      std::cout << "here" << std::endl;
      std::cout << Event::conf.gain[channel-1] << "  " << Event::conf.offset[channel-1] << std::endl;
      std::cout << val << "  " << cal << std::endl;
    */
    frontHits.reserve(32);
    backHits.reserve(4);
    if (Event::conf.side[channel-1] == 0) { //front
      frontHits.emplace_back(Event::conf.subID[channel-1], val, cal);
    }
    else if (Event::conf.side[channel-1] == 1) { //back
      backHits.emplace_back(Event::conf.subID[channel-1], val, cal);
    }
    else {
      std::cerr << "Warning! Channel " << channel << " with invalid front/back value" << std::endl;
    }
  }
  void QQQ5::AddHit(unsigned short int channel, unsigned short int value) {
    float val = (value-Event::conf.pedestal[channel-1] + Orruba::Event::Dither());
    float cal = val*Event::conf.gain[channel-1] + Event::conf.offset[channel-1];
    /*
      std::cout << "here" << std::endl;
      std::cout << channel << "   " << Event::conf.gain[channel-1] << "  " << Event::conf.offset[channel-1] << std::endl;
      std::cout << val << "  " << cal << std::endl;
    */
    if (Event::conf.side[channel-1] == 0) { //front
      frontHits.emplace_back(Event::conf.subID[channel-1], val, cal);
    }
    else if (Event::conf.side[channel-1] == 1) { //back
      backHits.emplace_back(Event::conf.subID[channel-1], val, cal);
    }
    else {
      std::cerr << "Warning! Channel " << channel << " with invalid front/back value" << std::endl;
    }
  }

  void QQQ5::MakeParticles(std::vector<SingleParticle*> &single_parts) {
    nGoodFronts = 0;
    nGoodBacks = 0;
    int frontIDs[2];
    int backIDs[2];
    float frontEns[2];
    float backEns[2];
    float frontThresh;
    float backThresh = 200;
    float verbose = false;
    /*
      for (int i=0; i<frontHits.size(); ++i) {
      if (frontHits[i].cal > 10000 && ID == 1) { verbose = true; }
      }
    */
    if (verbose) { std::cout << "=====================================" << std::endl; }
    if (ID <= 4) {
      frontThresh = 300;
    }
    else {
      frontThresh = 100;
    }
    for (int i=0; i<frontHits.size(); ++i) {
      if (verbose) {std::cout << "fh " << i << " " << frontHits[i].ID << "  " << frontHits[i].cal << std::endl;}
      if (frontHits[i].cal > frontThresh) { //use rough calibration here
        nGoodFronts += 1;
      }
      else { continue; }
      if (nGoodFronts > 2) { continue; }

      frontIDs[nGoodFronts-1] = frontHits[i].ID;
      frontEns[nGoodFronts-1] = frontHits[i].cal;      
    }
    for (int i=0; i<backHits.size(); ++i) {
      if (verbose) {std::cout << "bh " << i << " " << backHits[i].ID << "  " << backHits[i].cal << std::endl;}
      if (backHits[i].cal > backThresh) { //use rough calibration here
        nGoodBacks += 1;
      }
      else { continue; }
      if (nGoodBacks > 2) { continue; }

      backIDs[nGoodBacks-1] = backHits[i].ID;
      backEns[nGoodBacks-1] = backHits[i].cal;      
    }

    if ( nGoodFronts > 2 ) { return; }

      
    
    float frontEn = 0;      
    int frontID = 0;
      
    if (nGoodFronts == 1) {
      frontEn = frontEns[0];
      frontID = frontIDs[0];
    }
    else if (nGoodFronts == 2) {
      if (abs(frontIDs[0]-frontIDs[1]) == 1) {
        frontEn = frontEns[0] + frontEns[1];
        if (frontEns[0] > frontEns[1]) { frontID = frontIDs[0]; }
        else { frontID = frontIDs[1]; }
      }
      else {        
        //std::cout << "F" << frontIDs[0] << " : " << frontEns[0] << std::endl;
        //std::cout << "F" << frontIDs[1] << " : " << frontEns[1] << std::endl;
        return;
      }      
    }

    float backEn = 0;
    int backID = 0;

    if (nGoodBacks == 1) {
      backEn = backEns[0];
      backID = backIDs[0];
    }
    /*
      else if (nGoodBacks == 2) {
      if (abs(backIDs[0]-backIDs[1]) == 1) {
      backEn = backEns[0] + backEns[1];
      if (backEns[0] > backEns[1]) { backID = backIDs[0]; }
      else { backID = backIDs[1]; }
      }
      else {
      //std::cout << "B" << backIDs[0] << " : " << backEns[0] << std::endl;
      //std::cout << "B" << backIDs[1] << " : " << backEns[1] << std::endl;
      return;
      }      
      }
    */
    else {
      //search for back with closest energy to fired front
      float closest = 99999;
      for (int i=0; i<backHits.size(); ++i) {
        if (abs(backHits[i].cal - frontEn) < closest) {
          closest = abs(backHits[i].cal - frontEn);
          backID = backHits[i].ID;
          backEn = backHits[i].cal;
        }
      }
    }

    bool valid = (frontEn/backEn >=0.95 && frontEn/backEn<=1.05);
    QQQ5Particle *part = new QQQ5Particle(0, ID,
                                          frontID, backID,
                                          layer,
                                          frontEn, backEn,
                                          valid);
    part->MakeCoords(this);
    single_parts.push_back(part);

  }

  void QQQ5Particle::MakeCoords(QQQ5 *detector) {      
    if (detector->layer == 0) { theta = Event::conf.dEthetas[frontID]*M_PI/180.0; }
    else if (detector->layer == 1) { theta = Event::conf.Ethetas[frontID]*M_PI/180.0; }

    float sector_phi = 2.0*M_PI/(16.); //width
    float dphi = (backID)*sector_phi + sector_phi/2.;

    float sector_off = 0.0; //pointing down, +x

    //detector 0 -> 180, 1-> 270, 2-> 0, 3-> 90
    //this is ABCD clockwise from beam down, starting at top right
    //so A starts at +x, down, phi = 0
    //B starts at +y BL, phi = 90
    //C starts at -x up, phi = 180
    //D starts at -y BR, phi = 270
    sector_off = M_PI/2. * (float)((((detector->ID-1)%4) + 2)%4);

    phi = sector_off + dphi;
    r = Event::conf.radii[frontID]; //cylindrical r, not spherical polar
      
    if (detector->layer == 0) { z = Event::conf.dEdist; }
    else if (detector->layer == 1) { z = Event::conf.Edist; }

    x = r*std::cos(phi);
    y = r*std::sin(phi);        
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
