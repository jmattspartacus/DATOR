#ifndef LIBGRET_GRETINA_HH
#define LIBGRET_GRETINA_HH

#include <vector>
#include <map>
#include <iostream>

#include "Reader/Reader.hh"

#include "GRETINA/GretinaHit.hh"
#include "GRETINA/GretinaConf.hh"
#include "GRETINA/Gamma.hh"

#define GRETINA_MAXHITS 64

namespace Gret {
  class Event : public DATOR::Processor {
  public:
    int nhits;
    int ngammas;
    std::vector<GretinaHit> hits;
    std::vector<Gamma> gammas;
    //GretinaHit hits[GRETS800_GRETINA_MAXHITS];
    //Gamma gammas[GRETS800_GRETINA_MAXHITS];
    static Configuration conf;
    
    //counters which are not reset
    unsigned long long int nGretina;
    unsigned long long int nValidGretina;
    unsigned long long int nTotalGammas;
    unsigned long long int nEvents;
    
  public:
    Event() : nhits(0), ngammas(0) { hits.reserve(GRETINA_MAXHITS); gammas.reserve(GRETINA_MAXHITS); } ;
    void SetConf(Configuration c) { conf = c; }

    int AddHit(const int64_t GEBtimestamp,
               const crys_intpts *data);
    int BuildGammas();

    void ResetCounters() { nGretina = 0; nValidGretina = 0; nTotalGammas = 0; nEvents = 0; }
    void Reset() { nhits = 0; ngammas = 0; hits.clear(); gammas.clear(); }
    void Process(unsigned long long timestamp, unsigned short int *data, unsigned short int length);
    void ProcessFinal() { nEvents += 1; BuildGammas(); }
    void PrintSummary(std::ostream &out);

    ~Event() { }

  };
}

#endif
