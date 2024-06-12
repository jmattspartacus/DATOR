#ifndef LIBS800_RAWDATA_HH
#define LIBS800_RAWDATA_HH

#include <stdio.h>
#include <iostream>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define S800_CRDC_NPADS 224
#define S800_CRDC_MAX_SAMPLE 32
#define S800_CRDC_MAX_HITS 8192

#define S800_MTDC_MAXCHANS 32
#define S800_MTDC_MAXHITS 32
#define S800_MTDC_MAXGATES 5

#define INVMAP_MAX_ORDER 10
#define INVMAP_MAX_COEFFS 32
#define INVMAP_MAX_EXPONENTS 6
#define INVMAP_NBLOCKS 4

#define S800_CRDC_SEPARATION 1073 //mm

namespace S800 {
    enum class HeaderType {
                         kGretinaSignalDecomp = 1,
                         kGretinaTracked = 3,
                         kS800 = 5,
                         kS800Aux = 6,
                         kS800Physics = 9,
                         kS800AuxTS = 10
  };

  enum class S800Tag {
                      kTimeStamp = 0x5803,
                      kEventNumber= 0x5804,
                      kTrigger = 0x5801,
                      kTimeOfFlight = 0x5802,
                      kScintillator = 0x5810,
                      kIonChamber = 0x5820,
                      kICSubEnergy = 0x5821,
                      kCRDC = 0x5840,
                      kCRDCSubRaw = 0x5841,
                      kCRDCSubAnode = 0x5845,
                      kHodoscope = 0x58B0,
                      kHodoHitPat = 0x0002,
                      kTPPACs = 0x5870,
                      kTPPACSubRaw = 0x5871,
                      kObjPin = 0x58A0,
                      kFocalPin = 0x5805,
                      kGalotte = 0x58D0,
                      kLaBr = 0x58E0,
                      kMTDC = 0x58F0
  };
  
  struct GEBHeader {
    int32_t type;
    int32_t length; /* length of payload following the header, in bytes */
    int64_t timestamp;
  };

  static const int max_det = 30;
  static const int max_det_xtl = 4;
  static const int max_intpts = 16;
  static const int max_xtl = max_det*max_det_xtl;
  static const int nppacs = 20;

  struct crys_intpts {
    int type; /* as of June 2012: abcd5678 */
    int crystal_id;
    int num; /* # of interaction points from decomp, or # of nets on decomp error */
    float tot_e; /* CC energy for the central contact selected for use in decomp (calibrated, and for 10 MeV CC channels, includes DNL correction. */
    int core_e[4]; /* 4 raw core energies from FPGA filter (uncalibrated) */
    long long int timestamp; //(10 ns)
    long long int trig_time;
    float t0; //(1 ns)
    float cfd;
    float chisq;
    float norm_chisq;
    float baseline;
    float prestep; /* avg trace value before step (baseline) */
    float poststep; /* avg trace value after step (flat-top) */
    int pad; /* non-0 with a decomp error, value gives error type */ 
    /*    pad = 1   a null pointer was passed to dl_decomp()
          = 2   total energy below threshold
          = 3   no net charge segments in evt
          = 4   too many net charge segments
          = 5   chi^2 is bad following decomp (in this case
          crys_intpts is non-zero but post-processing step is not applied)
          = 6   bad build, i.e. <40 segment+CC channels found

          pad|= 128  PileUp, i.e. pileup flag or deltaT1<6usec

          e.g. pad = 128  pileup+Good
          = 133  pileup+BadChisq                */
    struct {
      float x,y,z,e; /* here e refers to the fraction */
      int seg; /* segment number hit */
      float seg_ener; /* energy of the hit segment */
    } intpts[max_intpts];
  };

  struct Mask {
    unsigned int BitMask;
    unsigned int BitShift;
    unsigned int operator()(const unsigned int &data) const {
      return ((data & BitMask) >> BitShift);
    }
  };
  
  struct tracked_gammas {
    int ngam;
    int pad;
    struct {
      float esum;              /* gamma ray energy */
      int ndet;                /* number of interactions */
      float fom;               /* figure of merit */
      int   tracked;           /* 1==if tracked */
      long long int timestamp; /* timestap of first interaction point */
      float x0, y0, z0, e0;    /* first interaction point */
      float x1, y1, z1, e1;    /* second interaction point */
    } gam[max_det];
  };
}


#endif
