#ifndef LIBGRET_GAMMA_HH
#define LIBGRET_GAMMA_HH

namespace Gret {
  class Gamma {
  public:
    double Energy;

    int hitInds[2];
    int nHits;
    
    int FirstInt;
    int ID;
    
    double Theta;
    double Phi;

    double Time;

    float Efficiency;

    bool Fix;
  public:
    Gamma() {};
    Gamma(double en, int nh, int *inds) :
      Energy(en), FirstInt(0), Theta(0), Phi(0), ID(0), Time(0), nHits(nh) { for (int i=0; i<nh; ++i) {hitInds[i] = inds[i];} }
    Gamma(double en, int nh, int *inds, int fi, double t, double p, int id, double time, bool fx, float eff) :
      Energy(en), FirstInt(fi), Theta(t), Phi(p), ID(id), Time(time), Fix(fx), Efficiency(eff), nHits(nh) {for (int i=0; i<nh; ++i) {hitInds[i] = inds[i];} }
    void Set(double en, int nh, int *inds, int fi, double t, double p, int id, double time, bool fx, float eff) {
      Energy = en;
      FirstInt = fi;
      Theta = t;
      Phi = p;
      ID = id;
      Time = time;
      Fix = fx;
      Efficiency = eff;
      nHits = nh;
      for (int i=0; i<nh; ++i) {hitInds[i] = inds[i];}
    }

  };
}

#endif
