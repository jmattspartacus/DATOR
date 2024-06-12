/*!
  \file BasicSort.cc 

  BasicSort is a simple example of a sort program using the DATOR framework. It utilizes the Reader for reading the GlobalMerged.dat data file(s), and three processors: one for Type-1 GEB events, and two for Type-19 GEB events.

  Simple output histograms are created, as well as a pruned Global.dat containing events only with high multipliticy events.

  Use this as a starting point for your own sort code.

To compile, copy this file and the Makefile in this directory to a "working directory" associated with your experimental analysis. 


## Understanding the Makefile

    INSTALLDIR=$(HOME)/.local/

This points the makefile at the directory where DATOR is installed.

    CC = g++

This is the C++ compiler program. On linux systems it will be ```g++```, on macOS it may be ```clang++``` instead.

    CFLAGS = `root-config --cflags` -O3 -g -fPIC

These are the default flags for compiling. The first section ```root-config --cflags``` gets the flags which the local ROOT was compiled with. The next three are various options that shouldn't need to be touched.

    LIBS = -Wl,--no-as-needed -lz

These are libraries to include. ```-lz``` refers to zlib, a compression library used to read compressed Globat.dat (i.e. Global.dat.gz) files.

    ROOTLIBS = `root-config --libs --glibs` -Wl,--no-as-needed -lMathMore

These are ROOT libraries to include, again, fetched using the ```root-config``` program.
    
    BasicSort : BasicSort.cc libReader libGRETINA libORRUBA
    	$(CC) $(CFLAGS) -I$(INSTALLDIR)/include/DATOR/ -o BasicDatorSort BasicSort.cc $(ROOTLIBS) $(LIBS) -L$(INSTALLDIR)/lib -lReader -lGRETINA -lORRUBA

This is the actual make command which puts them all together. The 

    -I$(INSTALLDIR)/include/DATOR/

gives the location of the DATOR header files, while 

    -L$(INSTALLDIR)/lib

points to the location of the DATOR library files (*.so).

    -lReader -lGRETINA -lORRUBA

Tells the compiler which libraries to actually link to: i.e. which ones are needed. In this case, we use the Reader (this will always be necessary), and then processor objects defined in the GRETINA and ORRUBA libraries.

## Understanding the sort code

\image html DATOR_sort.png
\image latex DATOR_sort.png

### Initialization
The first important section is

~~~~~~~~~~~~~~~{.cpp}
  DATOR::Reader reader;
  reader.warning = false;
  reader.coinc_window=3000;  //ns
  reader.PrunedOutput = true;
  reader.ts_mode = 0; //0 - GEB timestamps reset each run, 1 - GEB timestamps do not reset
~~~~~~~~~~~~~~~

where we initialize a DATOR::Reader object "reader", and set a few settings in it. We turn warnings off, set the coincidence window to 3000 ns, instruct it to write a the pruned Global.dat with default output location and naming, and set the timestamp mode to "0", indicating that the GEB timestamps will reset every run.

~~~~~~~~~~~~~~~{.cpp}
  std::cout << "Reading input files from " << ANSI_COLOR_BLUE << argv[1] << ANSI_COLOR_RESET << std::endl;
  reader.LoadPaths(argv[1]);
    
  if (argc > 3) {
    reader.SetRunNo(std::atoi(argv[3]));
  }
~~~~~~~~~~~~~~~

We then read the input file(s) from the first command line argument: If this is a *.txt file it will expect a list of files indexed by their run number: i.e.
   
    1   /path/to/run1/Global.dat
    2   /path/to/run2/Global.dat
    10   /path/to/run10/Global.dat

These do not need to be sequential or start from 1. Otherwise, if a single *.dat or *.dat.gz is passed as the first argument, the reader will sort only that file. 

What then follows are three sections where three different processors are initialized and added to the reader. The first is the type-1 GRETINA processor:

~~~~~~~~~~~~~~~{.cpp}
  Gret::Event gret;  //GRETINA processor
  reader.AddProcessor(1, &gret);
~~~~~~~~~~~~~~~

And then there are two type-19 ORRUBA processors

~~~~~~~~~~~~~~~{.cpp}
  Orruba::Event orr;   //ORRUBA processor
  reader.AddProcessor(19, &orr);

  Orruba::Basic orrb;  //ORRUBA "basic" processor
  reader.AddProcessor(19, &orrb);
~~~~~~~~~~~~~~~

The other lines in this section are initializing various parameters for the processing: details of the geometry, loading map and configuration files, etc. etc.

Next is a section where all the histograms are initialized

~~~~~~~~~~~~~~~{.cpp}
  TH1F *histogram = new TH1F("myHistogram", "Title",...);
~~~~~~~~~~~~~~~

### File loop

The first main loop is the "file loop", which goes through each file loaded into the reader. If there is only a single file in the *.txt, or a *.dat or *.dat.gz file is passed, this loop only occurs once.

~~~~~~~~~~~~~~~{.cpp}
  while (reader.NextFile()) {
~~~~~~~~~~~~~~~
The DATOR::Reader::NextFile() function handles the opening and closing of the file, and resetting various counters and diagnostics. At this point it would be a good time to do anything related to initializing anything specific to a particular file or run. For example, creating a histogram that is distinct for each run, or loading in a run-based calibration. 

### Event loop

Next is the main event loop. This achieved through the DATOR::Reader::Next() function, which does the time correlation and calls the ```Processor::Process``` and ```Processor::FinalProcess``` methods for each of the processors previously loaded into the reader. 

~~~~~~~~~~~~~~~{.cpp}
  while (reader.Next()) {
~~~~~~~~~~~~~~~

Inside this event all the main things happen: looping through the number of hits in each detector array, filling histograms accordingly, looking for time differences between the different arrays, applying cuts etc. In this example various GRETINA histograms are filled, including a gamma-gamma histogram, a GRETINA-ORRUBA time difference histogram is populated, and a basic diagnostic ORRUBA histogram is filled. In addition, the equivalence between the "Basic ORRUBA" and "ORRUBA" processors is checked: if they do not read in the same data, then an error is printed.

Finally, the DATOR::Reader::Write() function is called

~~~~~~~~~~~~~~~{.cpp}
  if (nGams > 3 && orr.nhits > 5) { reader.Write(); }
~~~~~~~~~~~~~~~

Which writes the entire physics event to the pruned Global.dat, but only in the event of more than three gamma rays in GRETINA, and more the 5 channels fired in ORRUBA. These are not particularly meaningful cuts, but serve to demsonstrate how the pruned Global.dat can work.

#### Things to try

Inside the event loop is where most additions will be made. To get a feel for how to extend this basic sort, you can try a few additions:
- print out some data from the ORRUBA (orr) or GRETINA (gret) processor objects. For example, 

~~~~~~~~~~~~~~~{.cpp}
  std::cout << gret.nhits << std::endl;
~~~~~~~~~~~~~~~

will print the number of GRETINA crystals hit in each event (this will spam the screen enormously, be prepared to Ctrl-C out of the sort program immediately). You can see the list of members of the Gret::Event class in Gretina.hh, and the corresponding Orruba::Event class in ORRUBA.hh. 
- create and fill a new histogram with data from either the ```orr``` or ```gret``` objects.
- check the ORRUBA-GRETINA time difference histogram, and make sure there's a clear prompt peak.
- create and fill a new histogram only when ORRUBA and GRETINA are present in the event, and when the time difference between them is within the prompt peak identified in the previous step.
- write your own processor for either ORRUBA or GRETINA (see MyProcessor.cc).


### Final things

At the end of the event loop, a summary of statistics is printed to screen with DATOR::Reader::PrintSummary(). This is where each ```Processor::PrintSummary``` function is called also. 

~~~~~~~~~~~~~~~{.cpp}
  reader.PrintSummary(std::cout);
~~~~~~~~~~~~~~~

This could also be written to a log file instead of ```std::cout```. After the file loop exits, the histograms are written to disk and the TFile object is closed

~~~~~~~~~~~~~~~{.cpp}
  file->Write();
  file->Close();
~~~~~~~~~~~~~~~

 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <tuple>
#include <signal.h>
#include <math.h>

#include "TFile.h"
#include "TH1.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include "TCutG.h"
#include "TObject.h"

#include "GRETINA/Gretina.hh"
#include "ORRUBA/ORRUBA.hh"

#include "Reader/Reader.hh"
#include "Reader/RawData.hh"

int quit;

void handler(int sig) {
  std::cout << "Ctrl-C caught, exiting cleanly" << std::endl;
  quit = 1;
}

int main(int argc, const char **argv) {
  if (argc < 2) {
    std::cout << "Usage: ./BasicDatorSort InputFileList OutputFileName [runno]" << std::endl;
    exit(1);
  }
  
  DATOR::Reader reader;
  reader.warning = false;
  reader.coinc_window=3000;  //ns
  reader.PrunedOutput = true;
  reader.ts_mode = 0; //0 - GEB timestamps reset each run, 1 - GEB timestamps do not reset

  int no;

  std::cout << "Reading input files from " << ANSI_COLOR_BLUE << argv[1] << ANSI_COLOR_RESET << std::endl;
  reader.LoadPaths(argv[1]);

  if (argc > 3) {
    reader.SetRunNo(std::atoi(argv[3]));
  }

  std::cout << std::endl;


  ///////////////////////
  // GRETINA PROCESSOR //
  ///////////////////////

  //configurations, map files, geometry, options etc.
  Gret::Configuration gretconf("GretinaConf", "GRETINA Configuration");
  gretconf.ReadCrystalMap("crmat.dat");
  gretconf.FX = Gret::FirstX::kHybrid;
  gretconf.Addback = true;
  //gretconf.warnings = true;
  gretconf.ClusterAngle = 15.0;
  gretconf.BeamLeftOffset = 20.8; //mm
  gretconf.BeamRightOffset = 23.7; //mm
  gretconf.PosXOffset = 0.0;
  gretconf.PosYOffset = 0.0;
  gretconf.PosZOffset = 0.0;
  //gretconf.PrintAngles();

  Gret::Event gret;
  gret.SetConf(gretconf);

  reader.AddProcessor(1, &gret);

  //////////////////////
  // ORRUBA PROCESSOR //
  //////////////////////

  //configurations, map files, geometry, options, etc.
  Orruba::Configuration orrconf("ORRUBAConf", "ORRUBA Configuration", "orruba.conf"); //this contains basic per-channel calibration
  orrconf.ReadCalibration("orruba.cal"); //this contains more complex calibrations for cross-talk effects
  orrconf.dEdist = -85.8;
  orrconf.Edist = -89.8;
  orrconf.ReadRadii("qqq5_radii.dat"); //radii of rings in barrel detectors
  orrconf.ReadPositionCalibration("sx3pos.cal"); //left and right edges of each SX3 in (r-l)/(r+l)

  Orruba::Event orr;
  orr.SetConf(orrconf);  

  reader.AddProcessor(19, &orr);

  ////////////////////////////
  // ORRUBA BASIC PROCESSOR //
  ////////////////////////////

  //this does nothing but load the (channel,value) pairs into vectors
  Orruba::Basic orrb;
  reader.AddProcessor(19, &orrb);
  
  std::cout << "Writing output histograms to " << ANSI_COLOR_YELLOW << argv[2] << ANSI_COLOR_RESET << std::endl;
  TFile *file = new TFile(argv[2], "RECREATE");

    //basic time differences
  TH1F *gret_orr_tdiff = new TH1F("gret_orr_tdiff", "Gretina/ORRUA time difference; #Delta t_{Gretina-ORRUBA} (ns); Counts", 1024, -5120, 5120);
  
  TH1F *delta_t1 = new TH1F("delta_t1", "All Sequential Time Difference; Time (ns); Events/10 ns", 4096, 0, 81920);
  TH1F *delta_t2 = new TH1F("delta_t2", "In-Event Time Difference; Time (ns); Events/20 ns", 4096, 0, 81920);

  //gretina
  TH2I *gret_rates = new TH2I("gret_rates", "Gretina Rates; Time (s); Crystal", 2*3600,0, 2*3600, 120, 0, 120);
  TH2F *e_cal = new TH2F("e_cal", "Energy;Energy(keV);Crystal ID", 8192, 0, 8192, 120,0,120);
  TH2F *e_raw = new TH2F("e_raw", "Energy (Raw);Energy(keV);Crystal ID", 8192, 0, 8192, 120,0,120);
  TH2F *gret_gamma_id = new TH2F("gret_gamma_id", "Gretina Gammas vs ID; Energy (keV); ID", 8192, 0, 8192, 120, 0, 120);
  TH1F *gret_gamma = new TH1F("gret_gamma", "Gretina Gammas; Energy (keV); Counts/1.0 keV", 8192, 0, 8192);  
  TH1F *gret_pad = new TH1F("gret_pad", "Gretina PAD (error) code; PAD; Counts", 256, 0, 256);
  TH2F *hit = new TH2F("hit", "Hit Pattern;Crystal ID;Index", 120, 0, 120, 200,0,200);
  TH2F *ID_ID = new TH2F("id_id", "ID vs ID; Crystal ID; Crystal ID", 120, 0, 120, 120, 0, 120);
  TH2F *gret_theta_phi = new TH2F("gret_theta_phi", "Theta vs Phi; Theta (deg); Phi (deg);", 360, 0, 360, 360, 0, 360);
  TH1F *gret_x = new TH1F("gret_x", "Gretina x (crystal coordinates); x (mm)", 400, -200, 200);
  TH1F *gret_y = new TH1F("gret_y", "Gretina y (crystal coordinates); y (mm)", 400, -200, 200);
  TH1F *gret_z = new TH1F("gret_z", "Gretina z (crystal coordinates); z (mm)", 400, -200, 200);
  
  TH2F *gret_dt_gg = new TH2F("gret_dt_gg", "Gamma vs dt; Energy (keV); Time (ns)", 4096, 0, 4096, 4096, -2048, 2048);
  TH2F *gret_gg = new TH2F("gret_gg", "Gamma vs Gamma (prompt); Energy (keV); Energy (keV)", 4096, 0, 4096, 4096, 0, 4096);
  TH2F *gret_gg_np = new TH2F("gret_gg_np", "Gamma vs Gamma (non-prompt); Energy (keV); Energy (keV)", 4096, 0, 4096, 4096, 0, 4096);

  //orruba 
  TH2I *raw_chans = new TH2I("raw_chans", "Raw ORRUBA channels; ADC Value; Channel", 4096, 0, 4096, 2048, 0.5, 2048.5);
  TH2I *ped_chans = new TH2I("ped_chans", "ORRUBA channels with pedestal subtraction; ADC Value; Channel", 4096, 0, 4096, 2048, 0.5, 2048.5);
  TH2I *cal_chans = new TH2I("cal_chans", "ORRUBA channels calibrated; Energy (keV); Channel", 4096, 0, 16384, 2048, 0.5, 2048.5);
  
  double ggprompt[2] = {-300,300}; //ns
  double ggnonprompt_1[2] = {-500,-300}; //ns
  double ggnonprompt_2[2] = {300,500}; //ns
  
  quit = 0;
  signal(SIGINT, handler);

  unsigned long long int evtCtr;
  int nFiles = 0;
  while (reader.NextFile()) {
    if (quit == 1) { break; }
    std::cout << "Reading from input file " << ANSI_COLOR_BLUE << reader.GetPath() << ANSI_COLOR_RESET << std::endl;
    std::cout << "Cumulative wall time = " << ANSI_COLOR_BLUE << reader.GetGlobalWallTime() << ANSI_COLOR_RESET << " mins" << std::endl;
    int run_no = reader.GetRunNo();

    nFiles += 1;
    reader.Start();
    evtCtr = 0;

    while (reader.Next()) {
        
      if (quit == 1) { break; }
      evtCtr += 1;
      if (evtCtr % 7137 == 0) {
        reader.PrintUpdate(std::cout);
      }

      //GRETINA things
      for (int i=0; i<gret.nhits; ++i) {
        auto grethit = gret.hits[i];

        gret_pad->Fill(grethit.PAD);
        gret_x->Fill(grethit.x);
        gret_y->Fill(grethit.y);
        gret_z->Fill(grethit.z);
        e_raw->Fill(grethit.RawEnergy, grethit.CrystalID);
        e_cal->Fill(grethit.TotalEnergy, grethit.CrystalID);

        for (int j=i+1; j<gret.nhits; ++j) {
          auto grethit2 = gret.hits[j];
          ID_ID->Fill(grethit.CrystalID, grethit2.CrystalID);
        }

        gret_rates->Fill((grethit.timestamp)*10.0/1e9-reader.starttime*60.0, grethit.CrystalID);
      }

      int nGams = gret.ngammas;
      int nHits = gret.nhits;
      for (int i=0; i<nGams; ++i) {
        auto gam_i = gret.gammas[i];

        gret_theta_phi->Fill(gam_i.Theta*180.0/M_PI, gam_i.Phi*180.0/M_PI);
        gret_gamma->Fill(gam_i.Energy);
        gret_gamma_id->Fill(gam_i.Energy, gam_i.ID);

        //GAMMA-GAMMA
        for (int j=i+1; j<nGams; ++j) {          
          auto gam_j = gret.gammas[j];
          double dt = 10.0*(gam_i.Time - gam_j.Time)/2.0;  //why the divide by 2? units ns

          if (gam_i.Energy > 80 && gam_i.Energy < 2000 && gam_j.Energy > 80 && gam_j.Energy < 2000) {
            gret_dt_gg->Fill(gam_i.Energy, dt);
            gret_dt_gg->Fill(gam_i.Energy, -dt);
            gret_dt_gg->Fill(gam_j.Energy, dt);
            gret_dt_gg->Fill(gam_j.Energy, -dt);            
          }

          if (ggprompt[0] < dt && dt < ggprompt[1]) {
            gret_gg->Fill(gam_i.Energy, gam_j.Energy);
            gret_gg->Fill(gam_j.Energy, gam_i.Energy);
          }

          if ((ggnonprompt_1[0] < dt && dt < ggnonprompt_1[1]) ||
              (ggnonprompt_2[0] < dt && dt < ggnonprompt_2[1])) {
            gret_gg_np->Fill(gam_i.Energy, gam_j.Energy);
            gret_gg_np->Fill(gam_j.Energy, gam_i.Energy);
          }
        }
        
        if (orr.timestamp) {
          double dt = 10.0*((long long)gam_i.Time - (long long)orr.timestamp)/2.0;  //why the divide by 2? units ns
          gret_orr_tdiff->Fill(dt);          
        }
      } //gamma loop

      //check equivalence between basic and full ORRUBA processing
      if (orr.nhits != orrb.nhits) { std::cerr << "Basic and full ORRUBA processing are not equivalent!" << std::endl; std::cerr << "nhits = " << orr.nhits << "  " << orrb.nhits << std::endl; }
      for (int i=0; i<orr.nhits; ++i) {
        if (orr.chans[i] != orrb.chans[i]) { std::cerr << "Basic and full ORRUBA processing are not equivalent!" << std::endl;
          std::cerr << "chan " << i << " = " << orr.chans[i] << "  " << orrb.chans[i] << std::endl; }
        if (orr.vals[i] != orrb.vals[i]) { std::cerr << "Basic and full ORRUBA processing are not equivalent!" << std::endl;
        std::cerr << "val " << i << " = " << orr.vals[i] << "  " << orrb.vals[i] << std::endl;}
      }

      //basic ORRUBA histograms
      for (int i=0; i<orr.nhits; ++i) {
        int chan = orr.chans[i];
        float val = orr.vals[i];
        raw_chans->Fill(val, chan);
        val -= orrconf.pedestal[chan-1];
        ped_chans->Fill(val, chan);
        double cal = (val + Orruba::Event::Dither())*orrconf.gain[chan-1] + orrconf.offset[chan-1];
        cal_chans->Fill(cal, chan);
      }

      if (nGams > 3 && orr.nhits > 5) { reader.Write(); }
    } //event loop
    
    reader.Stop();
    std::cout << std::endl;

    reader.PrintSummary(std::cout);
  } //file loop

  std::cout << "Sorted " << nFiles << " files" << std::endl;
  
  std::cout << "Writing histograms..." << std::flush;

  file->Write();
  file->Close();
  std::cout << " ...done" << std::endl;

  return 0;
}


