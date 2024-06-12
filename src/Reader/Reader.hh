#ifndef DATOR_READER_HH
#define DATOR_READER_HH

#include <stdint.h>
#include <stdio.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <chrono>
#include <zlib.h>

#include "RawData.hh"

namespace DATOR {
  class Processor {
  public:
    virtual void Reset() = 0;
    virtual void Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length) = 0;
    virtual void ProcessFinal() = 0;
    virtual void PrintSummary(std::ostream &out) {};
  };
  
  class Reader {
  public:
    FILE *DataFile;
    gzFile gDataFile;

    FILE *PrunedFile;
    gzFile gPrunedFile;
    
    int64_t old_timestamp;
    int coinc_window;
    static double Dither;

    int subevt;
    GEBHeader headers[MAX_GEB_SUBEVTS];
    unsigned int sub[MAX_GEB_SUBEVTS][MAX_GEB_PAYLOAD];

    int fileIndx = -1;
    std::vector<std::string> runPaths;
    std::vector<std::string> prunedPaths;
    std::vector<int> runNos;
    bool fCompressed;
    bool PrunedOutput;
    std::string prunedPathRoot;
    std::string prunedName;
    bool eof;
    
    unsigned long long int nOutOfOrder;
    unsigned long long int nEvents;
    unsigned long long int nGEBTypes[MAX_GEB_TYPE];

    long long eventmult;

    size_t fileSize;

    int ts_mode; // 0 - GEB timestamps reset every run, 1 - GEB timestamps never reset
    
    double starttime; //start of current run, in minutes
    double walltime; //timestamp of current event, in minutes
    double run_wt_offset; //in minutes, for last file

    std::chrono::system_clock::time_point start_time;
    std::chrono::system_clock::time_point stop_time;

    std::chrono::system_clock::time_point time_last;
    std::chrono::system_clock::time_point time_now;
    
    bool warning;

    std::vector<Processor*> processors[MAX_GEB_TYPE]; //these process each GEB type
    
    
  public:    
    Reader() : DataFile(0),
	       gDataFile(0),
               fileIndx(-1),
               nOutOfOrder(0),
               nEvents(0),
               old_timestamp(-1),
               coinc_window(171),
               ts_mode(0),
               walltime(0),
               run_wt_offset(0),
               warning(false),
               eof(false),
               fCompressed(false),
               PrunedOutput(false),
               prunedName("GlobalPruned"),
               prunedPathRoot(".")
    {}
    
    Reader(std::string fn) : Reader() { DataFile = fopen(fn.c_str(), "r"); gDataFile = gzdopen(fileno(DataFile), "r"); }
    static double GetDither();
    int LoadPaths(std::string fn);
    int NextFile();
    void AddProcessor(int type, Processor* proc);
    void PrintUpdate(std::ostream &out);
    void Reset();
    int Init();
    void Start();
    void Stop();
    int Read();
    int GetRunNo();
    double GetWallTime() { return walltime; }
    double GetGlobalWallTime() { if (ts_mode == 1) { return walltime; } else if ( ts_mode == 0 ) { return run_wt_offset + walltime; } };
    double GetRunWallTime() { if (ts_mode == 1) { return (walltime - run_wt_offset); } else if ( ts_mode == 0 ) { return walltime; } };
    void SetRunNo(int rn) { runNos.clear(); runNos.push_back(rn); }
    std::string GetPath();
    int currIndex;
    int Next();
    int Write();
    int PrintSummary(std::ostream &out);
  };
}

#endif
