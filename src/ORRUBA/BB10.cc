#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include "ORRUBA/ORRUBA.hh"
#include "ORRUBA/BB10.hh"

namespace Orruba {

  BB10::BB10(unsigned short int channel, unsigned short int value) {
    ID = Event::conf.detID[channel-1];
    layer = Event::conf.layer[channel-1];
    uds = Event::conf.UDS[channel-1];
    float val = (value-Event::conf.pedestal[channel-1] + Orruba::Event::Dither());
    float cal = val*Event::conf.gain[channel-1] + Event::conf.offset[channel-1];
    frontHits.reserve(10);
    frontHits.emplace_back(Event::conf.subID[channel-1], val, cal);
  }
  int BB10::AddHit(unsigned short int channel, unsigned short int value) {
    float val = (value-Event::conf.pedestal[channel-1] + Orruba::Event::Dither());
    float cal = val*Event::conf.gain[channel-1] + Event::conf.offset[channel-1];    
    if (val < Event::conf.threshold[channel-1]) { return 0; }
    frontHits.emplace_back(Event::conf.subID[channel-1], val, cal);
    return 1;
  }

  BB10Particle::BB10Particle(DetType dt, unsigned short int did,
                 unsigned short int sid,
                 unsigned short int lay,
                 float se, bool val) :
    SingleParticle(dt, did, sid, 0, lay, se, 0, val) {}
  
  void BB10Particle::MakeCoords(BB10 *detector) {
    //this r is in the cylindrical coordinate system

    float perp_r = 0;   
    if (detector->layer == 0) { r = Event::conf.r_dE_strip_bb10[frontID]; perp_r = Event::conf.r_dE_barrel; }
    else if (detector->layer == 1) { r = Event::conf.r_dE_strip_bb10[frontID]; perp_r = Event::conf.r_E_barrel; }

    //theta in spherical polar
    float position = 75./2.*detector->uds; //center of detector, this is pretty meaningless
    theta = std::atan(r/position);
    if (theta < 0 ) { theta += 3.14159265358799323846; }      
      
    float deltaphi = std::atan(Event::conf.sx3width/16.*std::abs((int)(7-2*frontID))/perp_r);
    int sign = 1;
    if (detector->layer == 1) {
      if (frontID < 4) { sign = -1; } //low strip ID = high phi for E
      else if (frontID >= 4) { sign = 1; }
    }
    else if (detector->layer == 0) {
      if (frontID < 4) { sign = 1; } //other way around for dE layer
      else if (frontID >= 4) { sign = -1; }
    }
    if (detector->uds > 0) { //downstream, reverse strip since detectors are flipped
      sign = -sign;
    }
    
    //phi of BB10 centre, common between cylindrical + spherical systems
    phi = (detector->ID+6)%12 * 2.*3.14159265/12.; //0-> down in lab frame-> +x, 6->up
    phi += sign*deltaphi;
    if (phi<0) { phi += 2.*3.14159265; }

    z = position;
    x = r*std::cos(phi);
    y = r*std::sin(phi);     
    
  }

  int BB10::MakeParticles(std::vector<SingleParticle*> &single_parts) {
    int nparts = 0;
    for (int i=0; i<frontHits.size(); ++i) {      
      BB10Particle *part = new BB10Particle(DetType::BB10, ID,
                                            frontHits[i].ID,
                                            layer,
                                            frontHits[i].cal, true);
      part->MakeCoords(this);
      single_parts.push_back(part);
      nparts += 1;
    }
    return nparts;
  }
}
