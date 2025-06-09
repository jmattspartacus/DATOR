#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include "ORRUBA/ORRUBA.hh"
#include "ORRUBA/SX3.hh"

namespace Orruba {

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
  int SX3::AddHit(unsigned short int channel, unsigned short int value) {
    float val = (value-Event::conf.pedestal[channel-1] + Orruba::Event::Dither());
    float cal = val*Event::conf.gain[channel-1] + Event::conf.offset[channel-1];    
    if (val < Event::conf.threshold[channel-1]) { return 0; }
    if (Event::conf.side[channel-1] == 0) { //front
      frontHits.emplace_back(Event::conf.LR[channel-1], Event::conf.subID[channel-1], val, cal);
    }
    else if (Event::conf.side[channel-1] == 1) { //back
      backHits.emplace_back(Event::conf.subID[channel-1], val, cal);
    }
    else {
      std::cerr << "Warning! Channel " << channel << " with invalid front/back value" << std::endl;
      return 0;
    }
    return 1;
  }

  SX3Particle::SX3Particle(DetType dt, unsigned short int did,
                unsigned short int pid, unsigned short int sid,
                unsigned short int lay,
                float pe, float sl, float sr, float pos, bool val) :
      SingleParticle(dt, did, sid, pid, lay, sl + sr, pe, val),
      stripL(sl), stripR(sr), position(pos) {}
  
  void SX3Particle::MakeCoords(SX3 *detector) {
    //this r is in the cylindrical coordinate system
    float perp_r = 0;
    if (detector->layer == 0) { r = Event::conf.r_dE_strip[frontID]; perp_r = Event::conf.r_dE_barrel; }
    else if (detector->layer == 1) { r = Event::conf.r_E_strip[frontID]; perp_r = Event::conf.r_E_barrel; }

    //theta in spherical polar
    theta = std::atan(r/position);
    if (theta < 0 ) { theta += 3.14159265358799323846; }      
      
    float deltaphi = std::atan(Event::conf.sx3width/8.*std::abs((int)(3-2*frontID))/perp_r);
    int sign = 1;
    if (detector->layer == 1) {
      if (frontID < 2) { sign = 1; } //low strip ID = high phi for E
      else if (frontID >= 2) { sign = -1; }
    }
    else if (detector->layer == 0) {
      if (frontID < 2) { sign = -1; } //other way around for dE layer
      else if (frontID >= 2) { sign = 1; }
    }
    if (detector->uds > 0) { //downstream, reverse strip ordering as detectors are flipped
      sign = -sign;
    }

    //phi of SX3 centre, common between cylindrical + spherical systems
    phi = (detector->ID+6)%12 * 2.*3.14159265/12.; //0-> down in lab frame-> +x, 6->up
    phi += sign*deltaphi;
    if (phi<0) { phi += 2.*3.14159265; }

    z = position;
    x = r*std::cos(phi);
    y = r*std::sin(phi);     
  }

  int SX3::MakeParticles(std::vector<SingleParticle*> &single_parts) {
    int frontPairs = 0;
    double frontSum = 0;
    double pos = 0;
    int stripID = 0;
    float stripL = 0;
    float stripR = 0;
    float stripLcal = 0;
    float stripRcal = 0;

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
            stripLcal = frontHits[j].cal;
            stripRcal = frontHits[i].cal;
            frontSum = frontHits[i].cal + frontHits[j].cal;
          }
          else if (end2 == 1 && end1 == 0) {
            frontPairs += 1;
            stripID = subID1;
            stripL = frontHits[i].value;
            stripR = frontHits[j].value;
            stripLcal = frontHits[i].cal;
            stripRcal = frontHits[j].cal;
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
      frontSum = stripL*stripLgain + stripLoff + stripR*stripRgain + stripRoff;
      //check energy equivalence
      pos = (stripRcal - stripLcal)/(stripRcal + stripLcal);
      pos = (pos - Event::conf.sx3cals[ID-1].stripPosOffset[stripID])*Event::conf.sx3cals[ID-1].stripPosGain[stripID] + Event::conf.sx3zoffset;
      if (frontSum==0) { return 0; }
      
      bool valid = (padEnergy/frontSum  >= 0.95 && padEnergy/frontSum <= 1.05);

      SX3Particle *part = new SX3Particle(DetType::SX3, ID,
                                          padID, stripID,
                                          layer,
                                          padEnergy, stripL, stripR, pos*uds, valid);
      part->MakeCoords(this);
      single_parts.push_back(part);
      return 1;
    }
    else { return 0; }
  }
}


