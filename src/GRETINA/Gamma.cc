#include "GRETINA/Gamma.hh"
#include "GRETINA/Gretina.hh"

namespace Gret {

  Gamma::Gamma(double en, int nh, int *inds) :
      Energy(en), FirstInt(0), Theta(0), Phi(0), ID(0), Time(0), nHits(nh) { 
	      hitInds.reserve(nh); 
	      for (int i=0; i<nh; ++i) {hitInds.push_back(inds[i]);} 
      }


  Gamma::Gamma(double en, int nh, int *inds, int fi, double t, double p, int id, long long int time, bool fx, float eff) :
      Energy(en), FirstInt(fi), Theta(t), Phi(p), ID(id), Time(time), Fix(fx), Efficiency(eff), nHits(nh) 
	{ 
		hitInds.reserve(nh); 
		for (int i=0; i<nh; ++i) { hitInds.push_back(inds[i]);} 
	}

  Gamma::Gamma(GretinaHit &hit, int indx) {
	Energy = hit.TotalEnergy;
        FirstInt = indx;
	Theta = hit.Theta;
 	Phi = hit.Phi;
	ID = hit.CrystalID;
	Time = hit.Time;
	Fix = hit.BadIntE || hit.BadT0;
	Efficiency = Event::conf.Efficiency(hit.TotalEnergy);
	nHits = 1;
	hitInds.push_back(indx);	
    }
    void Gamma::Set(double en, int nh, int *inds, int fi, double t, double p, int id, long long int time, bool fx, float eff) {
      Energy = en;
      FirstInt = fi;
      Theta = t;
      Phi = p;
      ID = id;
      Time = time;
      Fix = fx;
      Efficiency = eff;
      nHits = nh;
      hitInds.clear();
      for (int i=0; i<nh; ++i) {hitInds.push_back(inds[i]);}
    }
    void Gamma::AddHit(GretinaHit &hit, int indx) {
	    if (hit.TotalEnergy > Energy) { std::cout << "Warning! Hits added to gamma out of energy order" << std::endl; }
	    Energy += hit.TotalEnergy;
	    Fix = Fix || (hit.BadIntE || hit.BadT0);
	    Efficiency = Event::conf.Efficiency(Energy);
	    nHits += 1;
	    hitInds.push_back(indx);

    }
}
