#include <cmath>

#include "GRETINA/Gretina.hh"
#include "GRETINA/GretinaHit.hh"

namespace Gret {
  //Configuration Event::conf;
  int GretinaHit::Build(const int64_t GEBtimestamp,
                        const crys_intpts *data) {

    if (data->crystal_id < 4) { std::cerr << "bad Gretina hit with CrystalID="<<data->crystal_id<<"! ignoring" << std::endl; valid = false; return -1; }
    // IMPORTANT: Please note!
    //
    // cryst_intpts is a struct of fixed length, corresponding to the MAXIMUM number of interaction points (16).
    // However, we want to be able to read so-called "cropped" type-1 data, where the (non-existent) interaction point 
    // zeros are cropped from the end of each event
    //
    // To do this with minimal distruption, we remove the check on the payload size in Gretina.cc, and cast the
    // pointer to the (potentially incomplete) type-1 payload to (*cryst_intpts)
    //
    // This means that we can continue processing AS LONG AS we do not read the data->intpts array beyond element 
    // data->intpts[data->num-1]. Beyond this, behaviour is undefined and you should expect buffer overflow 
    // errors. User beware!
    //
    //
    valid = true;
    BadIntE = false;
    BadT0 = false;
    CrystalID = Event::conf.CrystalSwap[data->crystal_id];
    if (true) {
      if (CrystalID != data->crystal_id ) {
        std::cout << "Crystal " << data->crystal_id << " treated like " << CrystalID << std::endl;
      }
    }
    Hole = CrystalID/4;  //holes are one-indexed, CrystalID starts at 4
    PAD = data->pad;
    chisq = data->chisq;
    t0 = data->t0;
    timestamp = GEBtimestamp; //units 10 ns


    if (Event::conf.PAD128) {
      if (PAD != 0 && PAD != 128) { valid = false; }
    }
    else if (Event::conf.PAD128Fix) {
      if (PAD != 0 && PAD != 128) { valid = false; }      
      else if (PAD == 128 && Event::conf.NoisyCrystals[CrystalID] == 0) { valid = false; }
    }
    else {
      if (PAD != 0) { valid = false; }
    }

    // if (valid == false) {
    //   std::cout << "Bad PAD = " << PAD << ", crystal ID = " << CrystalID << std::endl;
    // }
    
    double offset = Event::conf.CalOffset[CrystalID];
    double gain = Event::conf.CalGain[CrystalID];

    RawEnergy = data->tot_e;
    TotalEnergy = (data->tot_e)*gain + offset;

    if (data->t0 != data->t0) {
      if (true) {
        std::cerr << "Warning! crystal t0 is NaN: " << data->t0 << std::endl;
      }
      BadT0 = true;
      Time = 10*data->timestamp + (long long int)Event::conf.MeanT0;
      //Time = data->timestamp;
      
    }
    else {
      Time = 10*data->timestamp + (long long int)data->t0;
      //Time = data->timestamp;
    }

    double minz = 1e6;
    double sum_energy = 0;
    MaxIntEn = 0;
    x = 0;
    y = 0;
    z = 0;


    //std::cout << TotalEnergy << "   " << RawEnergy << std::endl;
    nInteractions = data->num;
    if (data->num >= max_intpts) {
      if (Event::conf.warnings) {std::cerr << "Bad number of interaction points (" << data->num << "). PAD = " << data->pad << std::endl;}
      valid = false;        
    }
    else if (data->num == 0) {
      if (Event::conf.warnings) {std::cerr << "Bad number of interaction points (" << data->num << "). PAD = " << data->pad << std::endl;}
      valid = false;        
    }
    else {
      for (int i=0; i<data->num; ++i) {
        double int_x = data->intpts[i].x;
        double int_y = data->intpts[i].y;
        double int_z = data->intpts[i].z;
        double int_e = data->intpts[i].e;

	SegEnergy[i] = data->intpts[i].seg_ener;
	SegID[i] = data->intpts[i].seg;

        //std::cout << i << "   " << data->intpts[i].seg << "   " << data->intpts[i].seg_ener << "   " << data->intpts[i].e  << "   " << int_x << "   " << int_y << "   " << int_z << std::endl;

        if (BadT0) {
          int_x = 0;
          int_y = 0;
          int_z = 0;
          int_e = 1;
        }
      
        if (int_e <= 0) {
          if (Event::conf.warnings) {
            std::cerr << "Interaction (" << i << "/" << data->num << ") energy < 0! " << std::endl;
            std::cerr << "   Interaction energy = " << int_e << std::endl;
            std::cerr << "   Segment energy = " << data->intpts[i].seg_ener << std::endl;
            std::cerr << "   PAD = " << data->pad << std::endl;
          }
          int_x = 0;
          int_y = 0;
          int_z = 0;
          int_e = 1;          
          BadIntE = true;
        }      
      
        switch(Event::conf.FX) {        
        case FirstX::kSmallestZ: 
          //find smallest z
          if (int_z < minz) {
            minz = int_z;
            x = int_x;
            y = int_y;
            z = int_z;
          }
          break;
        case FirstX::kMeanPos: 
          x += int_x/(double)data->num;
          y += int_y/(double)data->num;
          z += int_z/(double)data->num;
          break;
        case FirstX::kEnergyWeighted: 
          x += int_x * int_e;
          y += int_y * int_e;
          z +=  int_z * int_e;
          sum_energy += int_e;
          break;
        case FirstX::kCrystalCent: 
          x = 0;
          y = 0;
          z = 0;
          break;
        case FirstX::kMaxEnergy: 
          if (int_e > MaxIntEn) {
            MaxIntEn = int_e;
            x = int_x;
            y = int_y;
            z = int_z;            
          }
          //std::cout << x << "  " << y << "  " << z << std::endl;
          break;
        case FirstX::kHybrid: 
          if (data->tot_e > 515) { //arbitrary, should this be a parameter?
            if (int_e > MaxIntEn) {
              MaxIntEn = int_e;
              x = int_x;
              y = int_y;
              z = int_z;            
            } 
          }
          else {
            x += int_x * int_e;
            y += int_y * int_e;
            z +=  int_z * int_e;
            sum_energy += int_e;                
          }
          break;  
        }
      }
        
      if (sum_energy > 0.0) {
        x /= sum_energy;
        y /= sum_energy;
        z /= sum_energy;
      }
    }
        
    
    if (x != x || y != y || z != z) {
      std::cout << "tote = " << data->tot_e << std::endl;
      std::cout << "x/y/z = " << x << " " << y << " " << z << std::endl;
    }

    auto pos = Event::conf.GetPosition(Hole, CrystalID%4, x,y,z);    
    PosX = pos[0];
    PosY = pos[1];
    PosZ = pos[2];

    double x2y2 = std::sqrt(PosX*PosX + PosY*PosY);
    double x2y2z2 = std::sqrt(PosX*PosX + PosY*PosY + PosZ*PosZ);
    double costhet = PosZ/x2y2z2;
    if (costhet >= -1.0 and costhet <= 1.0) {
      Theta = std::acos(costhet);
      //std::cout << "Crystal = " << CrystalID << "   Theta = " << Theta*180.0/3.1415926535 << std::endl;
    }
    else {
      std::cout << "Cos(theta) was out of range: " << costhet << std::endl;
      std::cout << PosX << "   " << PosY << "   " << PosZ << std::endl;
      std::cout << "Crystal " << CrystalID << "   PAD = " << data->pad << std::endl;
      std::cout << "FX = " << (int)(Event::conf.FX) << std::endl;
      Theta = 0;
    }
    double cosphi = PosX/x2y2;
    if (cosphi >= -1.0 and cosphi <= 1.0) {
      Phi = std::acos(cosphi);
      if (PosY < 0) {
        Phi = 2.0*3.1415926536 - Phi;
      }
    }
    else {
      std::cout << "Cos(phi) was out of range: " << cosphi << std::endl;
      Phi = 0;
    }

    // if (142 <= Theta*180.0/3.1415926535 && Theta <=143*180.0/3.1415926535) {
    //   if (193 <= Phi*180.0/3.1415926535 && Phi*180.0/3.1415926535 <= 194) {
    //     std::cout << Hole << "  " << CrystalID%4 << std::endl;
    //     std::cout << x << "   " << y << "   " << z << std::endl;
    //     std::cout << PosX << "   " << PosY << "   " << PosZ << std::endl;
    //     if (Fix) { std::cout << "Fix" << std::endl; }
    //   }
    // }

    if (!(Event::conf.FixValid)) {
      if (BadT0 || BadIntE) { valid = false; }
    }

    if (RawEnergy < 50) { valid = false; } //energy threshold

    return 0;
  }
}
