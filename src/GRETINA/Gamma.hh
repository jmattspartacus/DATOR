#ifndef LIBGRET_GAMMA_HH
#define LIBGRET_GAMMA_HH

#include "GretinaHit.hh"

namespace Gret {

  class Gamma {
  public:
    double Energy;

    std::vector<int> hitInds;
    int nHits;
    
    int FirstInt;
    int ID;
    
    double Theta;
    double Phi;

    long long int Time;

    float Efficiency;

    bool Fix;
  public:
    Gamma() {};
    Gamma(double en, int nh, int *inds); 
    Gamma(double en, int nh, int *inds, int fi, double t, double p, int id, long long int time, bool fx, float eff); 
    Gamma(GretinaHit &hit, int indx);
    void Set(double en, int nh, int *inds, int fi, double t, double p, int id, long long int time, bool fx, float eff); 
    void AddHit(GretinaHit &hit, int indx);

  };
}

#endif
