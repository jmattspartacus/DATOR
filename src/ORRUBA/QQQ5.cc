#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include "ORRUBA/ORRUBA.hh"
#include "ORRUBA/QQQ5.hh"

namespace Orruba {
  
  QQQ5::QQQ5(unsigned short int channel, unsigned short int value) {
    ID = Event::conf.detID[channel-1];
    if (ID < 1) { std::cerr << "Warning! QQQ5 ID = " << ID << " for channel " << channel << ", check orruba.conf" << std::endl; }
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
  int QQQ5::AddHit(unsigned short int channel, unsigned short int value) {
    float val = (value-Event::conf.pedestal[channel-1] + Orruba::Event::Dither());
    float cal = val*Event::conf.gain[channel-1] + Event::conf.offset[channel-1];
    if (val < Event::conf.threshold[channel-1]) { return 0; }
    if (Event::conf.side[channel-1] == 0) { //front
      frontHits.emplace_back(Event::conf.subID[channel-1], val, cal);
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

  int QQQ5::MakeParticles(std::vector<SingleParticle*> &single_parts) {
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
    if (verbose) { std::cout << "============== " << ID << " =======================" << std::endl; }
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

    if ( nGoodFronts > 2 ) { return 0; }

    if ( nGoodFronts == 0) { return 0; }
    if ( nGoodBacks == 0) { return 0; }
    
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
        return 0;
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
    else if (nGoodBacks > 1) {
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
    QQQ5Particle *part = new QQQ5Particle(DetType::QQQ5, ID,
                                          frontID, backID,
                                          layer,
                                          frontEn, backEn,
                                          valid);
    part->MakeCoords(this);
    single_parts.push_back(part);
    return 1;
  }

  QQQ5Particle::QQQ5Particle(DetType dt, unsigned short int did,
                             unsigned short int rid, unsigned short int sid,
                             unsigned short int lay,
                             float re, float se, bool val) :
    SingleParticle(dt, did, rid, sid, lay, re, se, val) {}


  void QQQ5Particle::MakeCoords(QQQ5 *detector) {      
    if (detector->layer == 0) { theta = Event::conf.dEthetas[frontID]*M_PI/180.0; }
    else if (detector->layer == 1) { theta = Event::conf.Ethetas[frontID]*M_PI/180.0; }

    float sector_phi = 2.0*M_PI/(16.); //width
    float dphi = (backID)*sector_phi + sector_phi/2.;

    float sector_off = 0.0; //pointing down, +x

    if (detector->uds <= 0) {
      sector_off = M_PI/2. * (float)((((detector->ID-1)%4) + 2)%4);
    }
    else {
      sector_off = M_PI/2. * (float)((((detector->ID-1)%4) + 2)%4) - M_PI/4.;
    }

    phi = sector_off + dphi;
    r = Event::conf.radii[frontID]; //cylindrical r, not spherical polar
      
    if (detector->layer == 0) { z = Event::conf.dEdist; }
    else if (detector->layer == 1) { z = Event::conf.Edist; }
    if (detector->uds <= 0) { z = -z; theta = M_PI-theta; }

    x = r*std::cos(phi);
    y = r*std::sin(phi);        
  }

}
