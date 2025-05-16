#ifndef LIBORR_SINGLEPARTICLE_HH
#define LIBORR_SINGLEPARTICLE_HH

namespace Orruba {
  enum class DetType {
    NoType,
    QQQ5,
    SX3,
    BB10,
    Track,
    TDC
  };

class SingleParticle  {
  public:
    DetType detType; //0=QQQ5, 1=SX3, 2=BB10
    unsigned int detID;
    unsigned int frontID; //strip/ring
    unsigned int backID; //pad/sector
    unsigned int layer; //dE/E

    float frontEnergy;
    float backEnergy;
      
    bool valid;

    //coordinate system is gretina system: +z along beam axis, +x towards floor, +y towards beam left
    float r; //radius
    float x,y,z; //cartesian
    float theta; //in spherical coordinates, angle from +z axis
    float phi; //azimuth, angle from +x towards +y

    float r_off;
    float x_off, y_off, z_off;
    float theta_off;
    float phi_off;

    SingleParticle() {}
    SingleParticle(DetType dt, unsigned short int did,
                   unsigned short int fid, unsigned short int bid,
                   unsigned short int lay,
                   float fe, float be, bool val) :
      detType(dt), detID(did), frontID(fid), backID(bid), layer(lay),
      frontEnergy(fe), backEnergy(be), valid(val) {};

    void OffsetBeam(float beamx, float beamy, bool verbose=false);
    
  };
}

#endif
