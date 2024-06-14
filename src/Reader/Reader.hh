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
  /*!  Prototype for Processor objects. All processors must inherit from this class, and define their own implementation of the Reset(), Process(), and ProcessFinal() functions.
  */
  class Processor {
  public:
    /*! Reset anything that should be zeroed/reset at the beginning of each physics event */
    virtual void Reset() = 0;
    /*! Do the processing of the raw GEB payload data
      \param timestamp GEB timestamp
      \param data GEB payload
      \param length Length of GEB payload in bytes
    */
    virtual void Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length) = 0;
    /*! Do final processing at the end of the physics event, after all subevents have been loaded. For example, addback for GRETINA
     */
    virtual void ProcessFinal() = 0;
    /*! Print summary of statistics/diagnostics, typically at the end of a run or file. */
    virtual void PrintSummary(std::ostream &out) {};
  };

  /*! Control object for I/O, time correlation (event building), and iterating through the data file(s). Processors are loaded into this class for use.

    Constructing a DATOR::Reader object is easy:

    \code{.cpp}
    DATOR::Reader reader;
    \endcode

    ## Basic useage
    Warnings are enabled/disabled with
    
    \code{.cpp}
    reader.warnings = true;  //or false
    \endcode

    Input files are loaded in with
    \code{.cpp}
    reader.LoadPaths("path/to/files");
    \endcode

    The type of input file is automatically determined based on its extension. A path ending in ".txt" will be interpreted as an ASCII file containing a list of (binary) data files to sort, each identified by a run number: 

    \code
    1   /path/to/run1/Global.dat
    2   /path/to/run1/Global.dat
    5   /path/to/run5/Global.dat
    \endcode

    The function ```Reader::NextFile()``` will switch to the next file in the list, closing the previous one and opening the next. The run number of the current file can be accessed with 

    \code{.cpp}
    reader.GetRunNo()
    \endcode

    which can be useful for creating run-by-run histograms, or loading run-by-run calibrations. Once a file has been loaded, iteration through the file proceeds with

    \code{.cpp}
    reader.Next()
    \endcode

    which will read (and process) sub-events from the file until a sufficiently large timestamp difference between successive events is found. The threshold for this timestamp difference is the coincidence window, set with

    \code{.cpp}
    reader.coinc_window = 3000; //ns
    \endcode
    
    A value of 3 us is fairly typical. The ```Reader::Next()``` function will return 1 if an event is successfully read, and 0 if the end of file is reached. Thus

    \code{.cpp}
    while (reader.Next()) { 
      //fill histograms here
    }
    \endcode
    
    is a convenient way to sort through a file. During the ```Reader::Next()``` function, the sub events are processed according to the DATOR::Processor objects loaded. These should be created and loaded before iterating through the file:

    \code{.cpp}
    Gret::Event gret; //GRETINA processor
    Orruba::Event orr; //ORRUBA processor

    reader.AddProcessor(1, &gret);
    reader.AddProcessor(19, &orr);
    \endcode

    Where the first argument to ```Reader::AddProcessor``` is the GEB type to be associated with that processor. 
    
    ### Single file input
    Single files can also be passed directly to the reader with Reader::LoadPaths. If the path ends fith ".dat", this will be interpreted directly as the path to the binary data file, and loaded (with run number zero). 

    ### Compressed input
    If the binary path (either passed directly or contained within the *.txt file) ends in ".dat.gz", it will be interpeted as a compressed binary file. The zlib library will be used to read this directly. If the input file is compressed, the (pruned) output file will also be written to a compressed file.

    ### Pruned output
    The DATOR::Reader class also provides support for writing out a subset of the raw data to a new set of files. To enable this feature, set

    \code{.cpp}
    reader.PrunedOutput = true;
    \endcode

    The output file is specified by a directory, and a file name.

    \code{.cpp}
    reader.prunedPathRoot = "/directory/for/pruned/output/";
    reader.prunedName = "GlobalPruned.dat";
    \endcode

    If multiple runs/files are loaded using a *.txt file, the reader will create a new pruned output file for every input file, demarkated with the suffix "_[run #].dat". The example *.txt file above would create
    
    \code
    GlobalPruned_1.dat
    GlobalPruned_2.dat
    GlobalPruned_5.dat
    \endcode

    To write to the pruned output, call 

    \code
    reader.Write()
    \endcode
    
    during the event loop. This will write the current event to the pruned output file.

    ### Wall time and ts_mode

    The wall time in minutes is stored for each event, based on the 10 ns GEB timestamp. However, the wall time may be accessed by the user in two ways, depending on the desired case. 

    \code{.cpp}
    reader.GetGlobalWallTime();
    \endcode

    Will get the time in minutes since the beginning of the experiment. This is useful for plotting, for example the rates for the entire experiment all together. However, one may also desire to have the wall time since the start of the current run:

    \code{.cpp}
    reader.GetRunWallTime();
    \endcode

    The calculation of these quantities depends on whether the GEB timestamps are reset every run or not. If they are not, the former calculation is trivial, while the latter requires an offset given by the first event in each run. If the timestamps are reset every run, the latter (Run) wall time is trivial, while the former requires a cumulative offset obtained by adding the last timestamp from each of the former runs. The two modes can be selected with

    \code{.cpp}
    reader.ts_mode; // 0 - GEB timestamps reset every run, 1 - GEB timestamps never reset
    \endcode

    Obviously in the case where GEB timestamps are reset, if runs are sorted out of order or with missing periods of time, the so-called "global wall times" do not in reality correspond to a true wall time. Nonetheless, they can be useful as a global time that distinguishes between runs.

    ### Updates and sort timing
    The ```Reader::PrintUpdate``` function may be used to print an update of the sort progress. This will print the total event sorted so far, the percentage completion (based on file size), and an instantaneous rate estimate of events per second. The timing 

    ### Final summary
    The ```Reader::PrintSummary``` function may be used to print a final summary of the sorting, usually once for each file. It will show the total event sorted, the total sort time, and various other diagnostics relating to the GEB data file. It is at this point that the ```Processor::PrintSummary``` functions are called.

    For accurate timing information, the ```Reader::Start()``` and ```Reader::Stop``` functions should be called immediately before and after the event loop.    
   */
  class Reader {
  public:
    FILE *DataFile;  /**< Input data file pointer */
    gzFile gDataFile; /**< Input data file (compressed) */

    FILE *PrunedFile; /**< Output pruned data file */
    gzFile gPrunedFile; /**< Output pruned data file (compressed) */
    
    int64_t old_timestamp; /**< Timestamp for the previous subevent */
    int coinc_window; /**< Rolling coincidence window in nanoseconds */
    static double Dither; /**< Dither for general use */

    int subevt; /**< Sub-event counter for within current event */
    GEBHeader headers[MAX_GEB_SUBEVTS]; /**< Buffer of GEB headers in current event */
    unsigned int sub[MAX_GEB_SUBEVTS][MAX_GEB_PAYLOAD]; /**< Buffer of GEB payloads in current event */

    int fileIndx = -1; /**< Index for current file */
    std::vector<std::string> runPaths; /**< List of paths for all the loaded files */
    std::vector<std::string> prunedPaths; /**< List of paths for the pruned output */
    std::vector<int> runNos; /**< List of run numbers */
    bool fCompressed; /**< Flag for compression, set for each file */
    bool PrunedOutput; /**< Flag for whether pruned output is enabled */
    std::string prunedPathRoot; /**< Path for the folder where pruned outputs will go */
    std::string prunedName; /**< Name of pruned output file: this will be modified with the run number index */
    bool eof; /**< Flag for when the end of file has been reached */
    
    unsigned long long int nOutOfOrder; /**< Counter for out of order events */
    unsigned long long int nEvents; /**< Counter for total number of physics events */
    unsigned long long int nEventsLast; /**< Counter for number of physics events at the previous printed update (used to calculate rates) */
    unsigned long long int nGEBTypes[MAX_GEB_TYPE]; /**< Counter for the number of sub-events of each GEB type encountered */

    size_t fileSize; /**< File size, used for progress */

    int ts_mode; /**< Timestamp mode: 

                    - (0) GEB timestamps reset every run, 
                    - (1) GEB timestamps never reset
                 */
    
    double starttime; /**< Start of current run, in minutes, set by first sub-event */
    double walltime; /**< Timestamp of current event, in minutes */
    double run_wt_offset; /**< Offset in minutes for the end of the previous file. If ts_mode = 0 this cumulatively increases every file. If ts_mode = 1, this is set to the first timestamp in the current file and used to calculate the run wall time. */

    std::chrono::system_clock::time_point start_time; /**< Time sorting started, set by Start() */
    std::chrono::system_clock::time_point stop_time; /**< Time sorting stopped, set by Stop() */

    std::chrono::system_clock::time_point time_last; /**< Time of previous printed update */
    std::chrono::system_clock::time_point time_now; /**< Current time, used for printing update */
    
    bool warning; /**< Warning flag */

    std::vector<Processor*> processors[MAX_GEB_TYPE]; /**< List of processors, indexed by GEB type */
    
    
  public:    
    Reader() : DataFile(0),
               gDataFile(0),
               PrunedFile(0),
               gPrunedFile(0),
               old_timestamp(-1),
               coinc_window(3000),
               subevt(0),               
               fileIndx(-1),
               fCompressed(false),
               PrunedOutput(false),
               prunedPathRoot("."),
               prunedName("GlobalPruned"),
               eof(false),
               nOutOfOrder(0),
               nEvents(0),
               nEventsLast(0),
               fileSize(0),
               ts_mode(0),
               starttime(0),
               walltime(0),
               run_wt_offset(0),
               warning(false)
      
    {}
    
    Reader(std::string fn) : Reader() { DataFile = fopen(fn.c_str(), "r"); gDataFile = gzdopen(fileno(DataFile), "r"); }
    static double GetDither();
    int LoadPaths(std::string fn);
    int NextFile();
    void AddProcessor(int type, Processor* proc);
    void PrintUpdate(std::ostream &out);
    void Reset();
    void Start();
    void Stop();
    int GetRunNo();
    /*! Get wall time */
    double GetWallTime() { return walltime; }
    /*! Get global wall time - since beginning of experiment */
    double GetGlobalWallTime() { if (ts_mode == 1) { return walltime; } else if ( ts_mode == 0 ) { return run_wt_offset + walltime; } };
    /*! Get run wall time - since beginning of run */
    double GetRunWallTime() { if (ts_mode == 1) { return (walltime - run_wt_offset); } else if ( ts_mode == 0 ) { return walltime; } };
    /*! Sets the run number */
    void SetRunNo(int rn) { runNos.clear(); runNos.push_back(rn); }
    std::string GetPath();
    int Next();
    int Write();
    int PrintSummary(std::ostream &out);
  };
}

#endif
