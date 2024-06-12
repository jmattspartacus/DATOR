#include <cmath>
#include <cstring>
#include <sstream>
#include <iomanip>

#include "Reader/Reader.hh"
#include "S800/S800.hh"
#include "S800/S800Confs.hh"

namespace S800 {
  Configuration Event::conf;
  
  void Event::SetConf(Configuration c) {
    conf = c;
    crdcs[0].SetCal(conf.crdc.gainx[0], conf.crdc.offsetx[0], conf.crdc.gainy[0], conf.crdc.offsety[0]);
    crdcs[1].SetCal(conf.crdc.gainx[1], conf.crdc.offsetx[1], conf.crdc.gainy[1], conf.crdc.offsety[1]);
  }
    
  void Event::Reset() {
    nHits = 0;
    timeMult = 0;
    headerTime = 0;
    evtNumMult = 0;
    trigMult = 0;
    tofMult = 0;
    scintMult = 0;
    icMult = 0;
    ionChamber.insideROI = false;
    crdcMult = 0;
    crdc1Mult = 0;
    crdc2Mult = 0;
    mtdcMult = 0;
    obj=0;
    xfp=0;
    obj_corr=0;
    xfp_corr=0;
  }  

  void Event::Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length) {
    nS800Hits += 1;
    nHits += 1;
    headerTime = timestamp;
        
    int bytecount = 0;
    bytecount +=10;

    while (bytecount < length) {
      short unsigned int s800datalength = data[bytecount/2];
      short unsigned int s800datatag = data[bytecount/2+1];

      bytecount += 4;
          
      S800Tag tag(static_cast<S800Tag>(s800datatag));
      switch (tag) {
      case S800Tag::kTimeStamp: {
        ReadTime(&data[bytecount/2]);
        break;
      }
      case S800Tag::kEventNumber: {
        ReadEventNumber(&data[bytecount/2]);
        break;
      }
      case S800Tag::kTrigger: {
        ReadTrigger(&data[bytecount/2], s800datalength-2);            
        break;
      }
      case S800Tag::kTimeOfFlight: {
        ReadTOF(&data[bytecount/2], s800datalength-2);
        break;
      }
      case S800Tag::kScintillator: {
        ReadScint(&data[bytecount/2], s800datalength-2);
        break;
      }
      case S800Tag::kIonChamber: {
        ReadIonChamber(&data[bytecount/2], s800datalength-2);
        break;
      }
      case S800Tag::kCRDC: {
        int label = ReadCRDC(&data[bytecount/2], s800datalength-2);
        //crdcs[label].Print(std::cout);
        break;
      }
      case S800Tag::kHodoscope: {
        if (s800datalength > 2 && conf.warning) {
          std::cerr << "Warning! Hodoscope processing not implemented" << std::endl;
        }
        break;
      }
      case S800Tag::kTPPACs: {
        if (s800datalength > 2 && conf.warning) {
          std::cerr << "Warning! Hodoscope processing not implemented" << std::endl;
        }
        break;
      }
      case S800Tag::kObjPin: {
        if (s800datalength > 2 && conf.warning) {
          std::cerr << "Warning! Object Pin processing not implemented" << std::endl;
        }
        break;
      }
      case S800Tag::kFocalPin:{
        if (s800datalength > 2 && conf.warning) {
          std::cerr << "Warning! Focal Pin processing not implemented" << std::endl;
        }
        break;
      }
      case S800Tag::kGalotte:{
        if (s800datalength > 2 && conf.warning) {
          std::cerr << "Warning! Galotte processing not implemented" << std::endl;
        }
        break;
      }
      case S800Tag::kLaBr:{
        if (s800datalength > 2 && conf.warning) {
          std::cerr << "Warning! LaBr processing not implemented" << std::endl;
        }
        break;
      }
      case S800Tag::kMTDC: {
        ReadMTDC(&data[bytecount/2], s800datalength-2);
        //mtdc.Print(std::cout);
        break;
      }
      default:
        std::cout << "Unknown S800 data packet " << s800datatag << std::endl;
      }
      bytecount += (s800datalength-2)*2;
      //std::cout << "    tag = " << s800datatag << "   length = " << s800datalength << std::endl;                   
    }
  }

  void Event::PrintSummary(std::ostream &out) {
    out << "--------- S800 Summary --------" << std::endl;
    out << "    " << ANSI_COLOR_YELLOW << nValidS800 << ANSI_COLOR_RESET << "/" << ANSI_COLOR_GREEN << nS800Hits << ANSI_COLOR_RESET
        << " valid S800 (" << ANSI_COLOR_YELLOW << std::setprecision(4) << (float)nValidS800/(float)nS800Hits*100.0 << "%" << ANSI_COLOR_RESET << ")     ";
    out << "   " << ANSI_COLOR_YELLOW << nValidS800_IC << ANSI_COLOR_RESET << "/" << ANSI_COLOR_GREEN << nS800Hits_IC << ANSI_COLOR_RESET
        << " inside IC gate (" << ANSI_COLOR_YELLOW << std::setprecision(4) << (float)nValidS800_IC/(float)nS800Hits_IC*100.0 << "%" << ANSI_COLOR_RESET << ")" << std::endl;
  }

  void Event::ProcessFinal() {
    Validate();
  }
    
  int Event::ReadTime(short unsigned int *data) {      
    time = data[0] + ( ((unsigned long long int)data[1])<<16) + (((unsigned long long int)data[2])<<32) + (((unsigned long long int)data[3])<<48);
    timeMult += 1;
    return 0;
  }
  int Event::ReadEventNumber(short unsigned int *data) {      
    eventNumber = data[0] + ( ((unsigned long long int)data[1])<<16) + (((unsigned long long int)data[2])<<32);
    evtNumMult +=1;
    return 0;
  }
  int Event::ReadTrigger(short unsigned int *data, int len) {
    trigMult +=1;
    trig.nTrigTimes = 0;
    if (len == 0) { return -1; }
      
    trig.trigPattern = data[0];

    trig.nTrigTimes = len-1;
    for (int i=1; i<len; ++i) {
      trig.trigTimes[i-1] = data[i];
    }
    return 0;
  }
  int Event::ReadTOF(short unsigned int *data, int len) {
    tofMult +=1;
    tof.nTimes = len;
    if (len==0) { return -1; }
      
    for (int i=0; i<tof.nTimes; ++i) {
      tof.times[i] = (data[i] & 0x0FFF);
      tof.chans[i] = (data[i] & 0xF000)>>12;
    }
    return 0;
  }
  int Event::ReadScint(short unsigned int *data, int len) {
    scintMult +=1;
    scint.nScints = len/2;
    if (len==0) { return -1; }

    for (int i=0; i<scint.nScints; ++i) {
      int chan1 = (data[2*i] & 0xF000) >> 12;
      int chan2 = (data[2*i+1] & 0xF000) >> 12;
      if (chan1 != chan2) { std::cerr << "Warning! Scintillator channels not equal!" << std::endl; }
      scint.energies[i] = (data[2*i] & 0x0FFF);
      scint.times[i] = (data[2*i+1] & 0X0FFF);
      scint.chans[i] = chan1;
    }
    return 0;
  }
  int Event::ReadIonChamber(short unsigned int *data, int len) {
    icMult +=1 ;
    ionChamber.nChans = 0;
    if (len==0) { return -1; }

    int sublen = data[0];
    int subtag = data[1];
    if (static_cast<S800Tag>(data[1]) != S800Tag::kICSubEnergy) { std::cerr << "Severe error! Ion Chamber packet does not make sense!" << std::endl; exit(1); }

    ionChamber.nChans = sublen-2;
    ionChamber.dE = 0;
    for (int i=0; i<ionChamber.nChans; ++i) {
      ionChamber.chans[i] = (data[2+i] & 0xF000) >> 12;
      ionChamber.energies[i] = (data[2+i] & 0x0FFF);
      ionChamber.cal[i] = conf.ionChamber.gain[ionChamber.chans[i]] * ionChamber.energies[i];
      ionChamber.dE += ionChamber.cal[i]/16.0;
    }
    if (ionChamber.dE >= conf.ionChamber.ROI[0] &&
        ionChamber.dE <= conf.ionChamber.ROI[1]) {
      ionChamber.insideROI = true;
      nS800Hits_IC += 1;
    }

    if (len != sublen) { std::cerr << "Warning! Ion Chamber len and sublen do not match: " << len << " vs " << sublen << std::endl; }

    return 0;
  }

  int Event::ReadCRDC(unsigned short int *data, int len) {
    crdcMult +=1;
    if (len==0) { return -1; }
    unsigned short int label = data[0];
    if ( (label != 0)  && (label != 1) ) {std::cerr << "Warning! Invalid CRDC label " << label << std::endl; return -1; }
    crdcs[label].Reset();

    if (label == 0) { crdc1Mult += 1;}
    if (label == 1) { crdc2Mult += 1;}

    crdcs[label].label=label;
      
    unsigned short int rawLen = data[1];
    unsigned short int rawTag = data[2];

    if (static_cast<S800Tag>(rawTag) != S800Tag::kCRDCSubRaw) { std::cerr << "Warning! CRDC Sub-packet tag incorrect!" << std::endl; return -1; }
    if (data[3] != 0) { std::cerr << "Warning! Global threshold is not zero" << std::endl; return -1; }

    unsigned int sample;
    unsigned int channel;
    unsigned int connector;
    unsigned int energy;
    for (int i=4; i<rawLen+1; ++i) {
      if (i==rawLen) { continue; } //temporary to get exact agreement with GRUTINIZER code
      unsigned short int ctrl = (data[i] & 0x8000) >> 15;

      if (ctrl == 1) {
        channel = (data[i] & 0x003F);
        sample = (data[i] & 0x7FC0) >> 6;
      }
      else {
        energy = (data[i] & 0x03FF);
        connector = (data[i] & 0x0C00) >> 10;
        if (sample >= S800_CRDC_MAX_SAMPLE) {
          std::cout << "CRDC " << label+1 << " : pad " << channel + connector*64 << ", sample = " << sample << ", energy = " << energy << std::endl;
          continue;
        }
        if (crdcs[label].nHits >= S800_CRDC_MAX_HITS) { std::cerr << "Warning! More than " << S800_CRDC_MAX_HITS << " hits for CRDC " << label << ", ignoring data. Consider increasing S800_CRDC_MAX_HITS" << std::endl; crdcs[label].Print(std::cerr); continue; }

        unsigned int pad = channel + connector*64;
        //double cal_en = std::max(0.0, ((energy - conf.crdc.pedestals[label][channel+connector*64])+Reader::GetDither())*conf.crdc.gains[label][channel+connector*64]);
        double cal_en = std::max(0.0, ((energy - conf.crdc.pedestals[label][channel+connector*64]))*conf.crdc.gains[label][channel+connector*64]);
        /*
          crdcs[label].pads[crdcs[label].nHits] = pad;
          crdcs[label].samples[crdcs[label].nHits] = sample;
          crdcs[label].padEnergies[crdcs[label].nHits] = cal_en;
          crdcs[label].nHits += 1;
        */
        crdcs[label].pads.push_back(pad);
        crdcs[label].samples.push_back(sample);
        crdcs[label].padEnergies.push_back(cal_en);
        crdcs[label].nHits += 1;

        if (pad >= S800_CRDC_NPADS) { std::cerr << "Warning! Pad number >= " << S800_CRDC_NPADS << " encountered, " << pad << std::endl; }
        //crdcs[label].sampleNums[pad][crdcs[label].nSamples[pad]] = sample;
        //crdcs[label].energies[pad][crdcs[label].nSamples[pad]] = cal_en;
        crdcs[label].nSamples[pad] += 1;          
      }
    }

    unsigned short int anodeLen = data[rawLen+1];
    unsigned short int anodeTag = data[rawLen+2];
    if (anodeLen != 4) { std::cerr << "Warning! CRDC Anode sub-packet length incorrect! " << anodeLen << std::endl; return -1; }
    if (static_cast<S800Tag>(anodeTag) != S800Tag::kCRDCSubAnode) { std::cerr << "Warning! CRDC Anode sub-packet tag incorrect!" << std::endl; return -1; }
    crdcs[label].anodeEnergy = data[rawLen+3];
    crdcs[label].anodeTime = data[rawLen+4];

    if (len != rawLen+5) { std::cerr << "Warning! CRDC len and rawLen do not match: " << len << " vs " << rawLen+5 << std::endl; }

    if (crdcs[label].nHits == 0) {
      crdcMult -=1;
      if (label == 0) { crdc1Mult -= 1;}
      if (label == 1) { crdc2Mult -= 1;}
    }
      
    return label;
  }

  int Event::ReadMTDC(short unsigned int *data, int len) {
    mtdc.Reset();
    mtdcMult += 1;      

    if (trigMult <= 0) { std::cerr << "Warning! MTDC packet without trigger information!" << std::endl; mtdc.trigPattern = -1; }
    mtdc.trigPattern = trig.trigPattern;

    if (len==0) { return -1; }

    for (int i=0; i<len/2; ++i) {
      int chan = (data[2*i] & 0x00FF);
      if (chan >= S800_MTDC_MAXCHANS) { std::cerr << "Warning! MTDC Channel " << chan << " > " << S800_MTDC_MAXCHANS <<", consider increasing S800_MTDC_MAXCHANS" << std::endl; continue; }
      mtdc.nHits += 1;

      mtdc.chans.push_back(chan);
      mtdc.times.push_back(data[2*i+1]);
      mtdc.hitNum.push_back((data[2*i] & 0xFF00) >> 8);
      mtdc.mults[chan] += 1;
        
      /*int n = mtdc.nHits[chan];
        mtdc.hitNum[chan][n] = (data[2*i] & 0xFF00) >> 8;
        mtdc.times[chan][n] = data[2*i+1];
        mtdc.nHits[chan]+=1;
      */
    }

    mtdc.Validate(&(conf.mtdc));

    return 0;
  }

  int Event::Validate() {
    valid = true;
    validCode = 0;
    //check if everything is good
    if (nHits != 1) { valid = false; validCode = 1; return 0; }
    if (trigMult != 1) { valid = false; validCode = 2; return 0; }
    if (mtdcMult != 1) { valid = false; validCode = 3; return 0; }

    //we require valid Mesytec TDC values for E1 up/down and obj and xfp
    if (mtdc.validTimes[0] == 0 ) { validCode = 100; valid = false; return 0; }
    if (mtdc.validTimes[1] == 0 ) { validCode = 101; valid = false; return 0; }
    if (mtdc.validTimes[2] == 0 ) { validCode = 102; valid = false; return 0; }
    if (mtdc.validTimes[3] == 0 ) { validCode = 103; valid = false; return 0; }

    if (mtdc.validTimes[0] < 0 ) { validCode = 150; valid = false; return 0; }
    if (mtdc.validTimes[1] < 0 ) { validCode = 151; valid = false; return 0; }
    if (mtdc.validTimes[2] < 0 ) { validCode = 152; valid = false; return 0; }
    if (mtdc.validTimes[3] < 0 ) { validCode = 153; valid = false; return 0; }
      
    //std::cout << "invalid times" << std::endl;
    //std::cout << mtdc.validTimes[0] << "   " << mtdc.validTimes[1] << "   " << mtdc.validTimes[2] << "  " << mtdc.validTimes[3] << std::endl;
    //std::cout << mtdc.nHits[0] << "   " << mtdc.nHits[1] << "   " << mtdc.nHits[2] << "  " << mtdc.nHits[3] << std::endl;

    if (icMult != 1) { valid = false; validCode = 4; return 0; }
    if (crdc1Mult != 1) { valid = false; validCode = 5; return 0; }
    if (crdc2Mult != 1) { valid = false; validCode = 6; return 0; }

    crdcs[0].SetX();
    crdcs[0].SetY();
    crdcs[1].SetX();
    crdcs[1].SetY();

    crdcs[0].SetXCal();
    crdcs[0].SetYCal();
    crdcs[1].SetXCal();
    crdcs[1].SetYCal();

    //this function also calculates adopted obj and xfp times
    if ((mtdc.validTimes[0] - mtdc.validTimes[1]) < 860 ) { valid = false; validCode = 21; return 0; }
    if ((mtdc.validTimes[0] - mtdc.validTimes[1]) > 1050) { valid = false; validCode = 21; return 0; }

    //x-dependent adjustment to E1 down time
    //e1dn = mtdc.validTimes[1] + crdcs[1].GetXCal()*0.22762508 + 926.7319;
    e1dn = mtdc.validTimes[1] + crdcs[1].GetXCal()*conf.e1dn_slope + conf.e1dn_offset;
    e1up = mtdc.validTimes[0];
      
    double tav = (mtdc.validTimes[0] + e1dn)/2.;
    //double tav = mtdc.validTimes[0];
    //std::cout << mtdc.validTimes[0] << "   " << tav << std::endl;
    obj = mtdc.validTimes[3] - tav;
    xfp = mtdc.validTimes[2] - tav;

    nValidS800 += 1;
    if (ionChamber.insideROI == true) {
      nValidS800_IC += 1;
    }

    return 1;
  }
    
  bool Trigger::GetS800() { return (trigPattern & 0x0001); }
  bool Trigger::GetCoinc() { return (trigPattern & 0x0002) >> 1; }
  bool Trigger::GetExt1() { return (trigPattern & 0x0004) >> 2; }
  bool Trigger::GetExt2() { return (trigPattern & 0x0008) >> 3; }
  bool Trigger::GetSecondary() { return (trigPattern & 0x0010) >> 4; }

  void Trigger::Print(std::ostream &out) {
    out << "Trigger " << trigPattern << std::endl;
    out << "Times: " << std::endl;
    for (int i=0; i<nTrigTimes; ++i) {
      out << "   " << trigTimes[i] << std::endl;
    }
  }

  void TimeOfFlight::Print(std::ostream &out) {
    std::cout << nTimes << " TOFs " << std::endl;
    for (int i=0; i<nTimes; ++i) {
      out << "    " << chans[i] << "   " << times[i] << std::endl;
    }
  }

  void Scintillator::Print(std::ostream &out) {
    out << nScints << " scintillators" << std::endl;
    for (int i=0; i<nScints; ++i) {
      out << chans[i] << "   " << energies[i] << "   " << times[i] << std::endl;
    }
  }

  void IonChamber::Print(std::ostream &out) {
    out << nChans << " ion chamber channels" << std::endl;
    for (int i=0; i<nChans; ++i) {
      out << chans[i] << "   " << energies[i] << std::endl;
    }
  }
    
  void CRDC::Print(std::ostream &out) {
    std::cout << nHits << " data from CRDC" << std::endl;
    std::cout << "   anode energy = " << anodeEnergy << ", anode time = " << anodeTime << std::endl;
    for (int i=0; i<nHits; ++i) {
      std::cout << i << "   sample : " << samples[i] << "    pad : " << pads[i] << "   energy : " << padEnergies[i] << std::endl;
    }
  }

  CRDC::CRDC() {
    nHits = 0;
    padEnergies.reserve(S800_CRDC_MAX_HITS);
    samples.reserve(S800_CRDC_MAX_HITS);
    pads.reserve(S800_CRDC_MAX_HITS);
  }
      
  void CRDC::Reset() {
    nHits = 0;
    std::memset(&nSamples[0], 0, S800_CRDC_NPADS*sizeof(unsigned short int));
    padEnergies.clear();
    samples.clear();
    pads.clear();
  }

  void CRDC::MaxHit(int &maxVal, int &hitInd, int &maxPad) {
    maxVal = -1;
    /*
      for (int i=0; i<nHits; ++i) {
      if (pads[i] == 197) { continue; }
      if (padEnergies[i] > maxVal) { maxVal = padEnergies[i]; hitInd = i; maxPad = pads[i]; }
      }
    */
    /*for (int i=0; i<S800_CRDC_NPADS; ++i) {
      if (nSamples[i] < 2) {
      continue;
      }
      for (int j=0; j<nSamples[i]; ++j) {
      if (energies[i][j] > maxVal) { maxVal = energies[i][j]; hitInd = j; maxPad = i; }
      }
      }
    */

    for (int i=0; i<nHits; ++i) {
      if (nSamples[pads[i]] < 2) { continue; }
      if (padEnergies[i] > maxVal) { maxVal = padEnergies[i]; hitInd = i; maxPad = pads[i]; }
    }
    return;
  }

  double CRDC::SetX() {
    //int maxVal, hitInd, maxPad;      
    MaxHit(maxVal, hitInd, maxPad);

    int width = 14;
    int num_pads = 224;

    int lowpad = std::max(0, maxPad - width/2);
    int highpad = std::min(num_pads-1, maxPad + width/2);

    double sum = 0;
    double wsum = 0;

    int nbad = 0;
    for (int i=0; i<nHits; ++i) {
      if (nSamples[pads[i]] < 2) { 
        if (nSamples[pads[i]] > 0) {
          nbad += 1;
        }
        continue;
      }      
          
      if ((pads[i] < lowpad) || (pads[i] > highpad)) { continue; }

      sum += padEnergies[i];
      wsum += pads[i]*padEnergies[i];
    }

    /*
      if (nbad) {
      std::cout << "bad pads! " << std::endl;
      Print(std::cout);
      }
    */
      
    /*
      for (int i=0; i<nHits; ++i) {
      //       if (i==0 && nHits > 1) {
      //   if (pads[i] != pads[i+1]) {
      //     std::cout << "bad pad! " << i << std::endl;
      //     std::cout << pads[i] << "   " << pads[i+1] << std::endl;
      //     Print(std::cout);
      //     continue;
      //   }
      // }
      // else if (i==nHits-1 && nHits > 2) {
      //   if (pads[i] != pads[i-1]) {
      //     std::cout << "bad pad! " << i << std::endl;
      //     std::cout << pads[i] << "   " << pads[i-1] << std::endl;
      //     Print(std::cout);
      //     continue;
      //   }
      // }
      // else if (nHits > 2) {
      //   if (pads[i] != pads[i-1] && pads[i] != pads[i+1]) {
      //     std::cout << "bad pad! " << i << std::endl;
      //     std::cout << pads[i] << "   " << pads[i+1] << "   " << << std::endl;
      //     Print(std::cout);
      //     continue;
      //   }
      // }
        
      Int pad = pads[i];
      if ((pad < lowpad) || (pad > highpad)) { continue; }

      sum += padEnergies[i];
      wsum += pad*padEnergies[i];
      }
    */

      

    double mean = wsum/sum + 0.5;
    x = mean;

    return mean;
  }

  double CRDC::SetY() {
    y = anodeTime;//+Reader::GetDither();
    return anodeTime;
  }

  double CRDC::SetXCal() {
    xcal = x*gainx + offsetx;
    return x*gainx + offsetx;
  }
    
  double CRDC::SetYCal() {
    ycal = y*gainy + offsety;
    return y*gainy + offsety;
  }

  double CRDC::GetX() { return x; }
  double CRDC::GetY() { return y; }
  double CRDC::GetXCal() { return xcal; }
  double CRDC::GetYCal() { return ycal; }

  int Event::GetAFP(double &theta) {
    if ((crdc1Mult != 1) || (crdc2Mult != 1)) { return 0; }

    //std::cout << crdcs[0].GetXCal() << "   " << crdcs[1].GetXCal() << std::endl;

    theta = std::atan2((crdcs[1].GetXCal() - crdcs[0].GetXCal()),S800_CRDC_SEPARATION);
    return 1;
  }

  int Event::GetBFP(double &theta) {
    if ((crdc1Mult != 1) || (crdc2Mult != 1)) { return 0; }

    //std::cout << crdcs[0].GetXCal() << "   " << crdcs[1].GetXCal() << std::endl;

    theta = std::atan2((crdcs[1].GetYCal() - crdcs[0].GetYCal()),S800_CRDC_SEPARATION);
    return 1;
  }
   
  void Event::Print(std::ostream &out) {
    std::cout << "S800 Multiplicities: " << std::endl;
    std::cout << "   Time stamp    : " << timeMult << std::endl;
    std::cout << "   Event Number  : " << evtNumMult << std::endl;
    std::cout << "   Trigger       : " << trigMult << std::endl;
    std::cout << "   TOF           : " << tofMult << std::endl;
    std::cout << "   Scintillator  : " << scintMult << std::endl;
    std::cout << "   Ion Chamber   : " << icMult << std::endl;
    std::cout << "   CRDC (Total)  : " << crdcMult << std::endl;
    std::cout << "   CRDC 1        : " << crdc1Mult << std::endl;
    std::cout << "   CRDC 2        : " << crdc2Mult << std::endl;
    std::cout << "   Mesytec TDC   : " << mtdcMult << std::endl;      
  }  

  S800::MTDC::MTDC() {
    chans.reserve(S800_MTDC_MAXCHANS*S800_MTDC_MAXHITS);
    times.reserve(S800_MTDC_MAXCHANS*S800_MTDC_MAXHITS);
    hitNum.reserve(S800_MTDC_MAXCHANS*S800_MTDC_MAXHITS);
  }
    
  void S800::MTDC::Reset() {
    std::memset(&mults[0], 0, sizeof(mults[0])*S800_MTDC_MAXCHANS);
    nHits = 0;
    chans.clear();
    hitNum.clear();
    times.clear();
    //for (int i=0; i<S800_MTDC_MAXCHANS; ++i) {
    //  nHits[i] = 0;
    //}
  }

  void S800::MTDC::Print(std::ostream &out) {
    /*
      for (int i=0; i<S800_MTDC_MAXCHANS; ++i) {
      if (nHits[i] == 0) { continue; }
      std::cout << "Mesytec TDC channel " << i << std::endl;
      for (int j=0; j<nHits[i]; ++j) {
      std::cout << "   " << j << "   " << hitNum[i][j] << "   " << times[i][j] << std::endl;
      }
      }
    */
    out << "Mesytec TDC " << std::endl;
    for (int i=0; i<nHits; ++i) {
      out << i << "  " << chans[i] << "  " << hitNum[i] << "  " << times[i] << std::endl;
    }
      
  }

  void S800::MTDC::Validate(MTDCConf *conf) {
    std::memset(&validTimes[0], 0, sizeof(validTimes[0])*S800_MTDC_MAXCHANS);
    if (trigPattern <= 0) { std::cerr << "Warning! No trigger pattern in MTDC; cannot validate!" << std::endl;
      return;
    }

    /*
      for (int i=0; i<S800_MTDC_MAXCHANS; ++i) {
      //validTimes[i] = 0;
      for (int j=0; j<nHits[i]; ++j) {
      if (((times[i][j] > conf.low[trigPattern][i]) && (times[i][j] < conf.high[trigPattern][i]))) {
      if (validTimes[i] == 0) { validTimes[i] = times[i][j]; }
      else { validTimes[i] = -1; break; } //more than one TDC hit inside gates          
      }
      }
      }
    */
    for (int i=0; i<nHits; ++i) {
      int chan = chans[i];
      if (((times[i] > conf->low[trigPattern][chan]) && (times[i] < conf->high[trigPattern][chan]))) {
        if (validTimes[chan] == 0) { validTimes[chan] = times[i]; }
        else { validTimes[chan] = -1; break; } //more than one TDC hit inside gates          
      }
    }
  }

  S800::InvMap::Block::Block() {
    number = 0;
    for (int i=0; i<INVMAP_MAX_ORDER; ++i) {
      ncoeffs[i] = 0;
      for (int j=0; j<INVMAP_MAX_COEFFS; ++j) {
        coeffs[i][j] = 0;
        for (int k=0; k<INVMAP_MAX_EXPONENTS; ++k) {
          exponents[i][j][k] =0;
        }
      }
    }
  }

  S800::InvMap::InvMap() {
    for (int i=0; i<4; ++i) {
      blocks[i].number = i;
    }
  }

  int S800::InvMap::ReadInvMap(std::string fn) {
    FILE *file = fopen(fn.c_str(), "ra");
    std::stringstream ss;
    char cline[2048];

    int block = -1;
    std::fgets(cline, sizeof cline, file);
    while (std::fgets(cline, sizeof cline, file)!=NULL) {
      std::string line(cline);
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      if (line.find("-----------")!=std::string::npos) {
        continue;
      }
              
      if (line.find("COEFFICIENT")!=std::string::npos) {
        block += 1;
        continue;
      }
          
      ss.clear();
      ss.str(line);
      int indx;
      double coeff;
      int order;
      int e1,e2,e3,e4,e5,e6;

      ss >> indx >> coeff >> order >> e1 >> e2 >> e3 >> e4 >> e5 >> e6;

      int ind = blocks[block].ncoeffs[order];
      blocks[block].coeffs[order][ind] = coeff;
      blocks[block].exponents[order][ind][0] = e1;
      blocks[block].exponents[order][ind][1] = e2;
      blocks[block].exponents[order][ind][2] = e3;
      blocks[block].exponents[order][ind][3] = e4;
      blocks[block].exponents[order][ind][4] = e5;
      blocks[block].exponents[order][ind][5] = e6;
      blocks[block].ncoeffs[order] += 1; 
    }
    return 0;
  }

  void S800::InvMap::Print(std::ostream &out, int order) {
    for (int i=0; i<4; ++i) {
      out << "======= Block " << i << " ==============" << std::endl;
      for (int j=1; j<=order; ++j) {
        out << "   Order " << j << std::endl;
        for (int n=0; n<blocks[i].ncoeffs[j]; ++n) {
          out << "      " << n << "   " << blocks[i].coeffs[j][n] << "    " << blocks[i].exponents[j][n][0]<< " " << blocks[i].exponents[j][n][1] << "  " << blocks[i].exponents[j][n][2] << " " << blocks[i].exponents[j][n][3] << std::endl;
        }
      }
    }
  }
  
  double S800::InvMap::Eval(const int &block,
                            const int &order,
                            double *x) {

    double sum = 0;
    if (order < 1) { std::cerr << "Invalid order for InvMap!" << std::endl; return 0.0; }
    for (int o=1; o<=order; ++o) {
      for (int i=0; i<blocks[block].ncoeffs[o]; ++i) {
        double mult = 1.0;
        for (int j=0; j<4; ++j) {
          if (blocks[block].exponents[o][i][j]>0) {
            mult *= std::pow(x[j], blocks[block].exponents[o][i][j]);
          }
          //do we care about exponents 5 and 6?
        }
        sum += blocks[block].coeffs[o][i] * mult;          
      }
    }
    return sum;      
  }
}
