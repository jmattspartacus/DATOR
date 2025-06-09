#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <zlib.h>

#include "Reader.hh"

namespace DATOR {
  double Reader::Dither;

  /*! Loads files into the reader object
    \param fn File to load paths from. This can be either a *.txt file containing two (whitespace-separated) columns: run number and path
    ~~~~~~~~~~~~
    1    /path/to/run1/Global.dat
    2    /path/to/run2/Global.dat
    5    /path/to/run5/Global.dat
    ~~~~~~~~~~~~
    or a single *.dat (uncompressed) or *.dat.gz (compressed) data file
   */
  int Reader::LoadPaths(std::string fn) {
    if (!fn.substr(fn.size()-7).compare(".dat.gz")) {
      fCompressed = true;
      runPaths.clear();
      runNos.clear();
      prunedPaths.clear();

      runPaths.push_back(fn);
      runNos.push_back(0);

      prunedPaths.push_back(prunedPathRoot+"/"+prunedName+".dat.gz");
      return 0;
    }
    else if (!fn.substr(fn.size()-4).compare(".dat")) {
      fCompressed = false;
      runPaths.clear();
      runNos.clear();
      prunedPaths.clear();

      runPaths.push_back(fn);
      runNos.push_back(0);
      prunedPaths.push_back(prunedPathRoot+"/"+prunedName+".dat");
      return 0;
    }
    else if (!fn.substr(fn.size()-4).compare(".txt")) {
      std::ifstream infile(fn);
      if (!infile.is_open()) { std::cerr << fn << " not found!" << std::endl; return -1; }

      int indx;
      std::string path;
      // using a trailing string allows to have end of line comments
      std::string trailing;

      runPaths.clear();
      runNos.clear();
      std::string tline;
      while (std::getline(infile, tline)) {
        if(tline.size() < 2) { continue; }
        if(tline[0] == "#"[0]){
          continue;
        }
        std::stringstream st(tline);
        st >> indx >> path >> trailing;
        // uncomment if you want explicit output of the run numbers!
        // std::cout << "Run num: " << indx << " Path: " << path << " Trailing: " << trailing << std::endl;
        
        if (runNos.size()>0) {
          if (runNos[runNos.size()-1] == indx) {
            std::cerr << "Warning! Same run number encountered twice?" << std::endl;
          }
        }
        runPaths.push_back(path);
        runNos.push_back(indx);
        if (!path.substr(path.size()-7).compare(".dat.gz")) {
          prunedPaths.push_back(prunedPathRoot+"/"+prunedName+"_"+std::to_string(indx)+".dat.gz");
        }
        else if (!path.substr(path.size()-4).compare(".dat")) {
          prunedPaths.push_back(prunedPathRoot+"/"+prunedName+"_"+std::to_string(indx)+".dat");
        }                   
      }
      return 0;
    }
    else {
      std::cerr << "Unrecognized file extension!" << std::endl;
      exit(1);
    }
  }    

  /*! Switches to next file in the list of loaded files. Closes the old file and opens a new one, resetting and incrementing various counters and diagnostics. Also reads the first header of the new file in preparation for a Reader::Next() loop. 
 * 

    \return Error code. <br>
    0 - New file couldn't be opened or opened file is empty. <br>
    1 - No error, file opened and first header read successfully.
    
   */
  int Reader::NextFile() {
    fileIndx += 1;
    if (fCompressed) {
      if (fileIndx >= (int)runPaths.size()) {
        if (DataFile != 0) { gzclose(gDataFile); }
        if (PrunedFile != 0) { gzclose(gPrunedFile); }
        return 0;
      }
    }
    else {
      if (fileIndx >= (int)runPaths.size()) {
        if (DataFile != 0) { fclose(DataFile); }
        if (PrunedFile != 0) { fclose(PrunedFile); }
        return 0;
      }
    }

    if (fCompressed) {
      if (DataFile != 0) { gzclose(gDataFile); }
      if (PrunedFile != 0) { gzclose(gPrunedFile); }
    }
    
    DataFile = fopen(runPaths[fileIndx].c_str(), "r");
    if (!DataFile) { std::cerr << "File " << runPaths[fileIndx] << " does not exist!" << std::endl; return 0; }
    if (PrunedOutput) {
      if (!prunedPaths[fileIndx].compare(runPaths[fileIndx])) { std::cerr << "Pruned file (output) is the same as input file!" << std::endl; exit(1); }
      PrunedFile = fopen((prunedPaths[fileIndx]).c_str(), "w");
    }
    else {
      PrunedFile = NULL;
    }

    if (!runPaths[fileIndx].substr(runPaths[fileIndx].size()-7).compare(".dat.gz")) {
      fCompressed = true;
    }
    else if (!runPaths[fileIndx].substr(runPaths[fileIndx].size()-4).compare(".dat")) {
      fCompressed = false;
    }
    else {
      std::cout << "Warning! Unrecognized file extension" << std::endl;
    }
    
    if (fCompressed) {
      gDataFile = gzdopen(fileno(DataFile), "r");
      if (PrunedOutput) {
        gPrunedFile = gzdopen(fileno(PrunedFile), "w");
      }
    }
    
    if (ts_mode == 0) {
      run_wt_offset += walltime;
    }    
    
    Reset();

    std::fseek(DataFile, 0L, SEEK_END);
    fileSize = std::ftell(DataFile);
    std::rewind(DataFile);

    time_last = std::chrono::system_clock::now();
    eof = false;

    //read first header
    if (fCompressed) {      
      if (gzread(gDataFile, &headers[0], sizeof(headers[0])) == 0) {
        return 0;
      }
    }
    else {
      if (fread(&headers[0], sizeof(headers[0]), 1, DataFile) == 0) {
        return 0;
      }
    }
    
    return 1;
  }

  /*! Prints sort progress

    \param out Stream to print to
  */
  void Reader::PrintUpdate(std::ostream &out) {
    
    double perc_complete;
    if (fCompressed) {
      perc_complete = 100.0*(float)((float)gzoffset(gDataFile))/fileSize;
    }
    else {
      perc_complete = 100.0*(float)((float)ftell(DataFile))/fileSize;
    }
    auto time_now = std::chrono::system_clock::now();
    double evt_diff = nEvents - nEventsLast;
    double diff = double(std::chrono::duration_cast <std::chrono::microseconds> (time_now - time_last).count());
    std::cout << "\r" << ANSI_COLOR_GREEN << nEvents << ANSI_COLOR_RESET << " events sorted, " << std::fixed << std::setprecision(1) << ANSI_COLOR_GREEN << perc_complete << "%" << ANSI_COLOR_RESET ;
    out.unsetf(std::ios_base::floatfield);
    out << std::setprecision(6) << " [" << ANSI_COLOR_YELLOW << evt_diff/diff*1000000 << ANSI_COLOR_RESET << " evts/s]" << std::flush;
    time_last = time_now;
    nEventsLast = nEvents;
  }

  /*! Resets counters, usually in preparation for a new file */
  void Reader::Reset() {
    nOutOfOrder = 0;
    nEvents = 0;
    for (int i=0; i<MAX_GEB_TYPE; ++i) {
      nGEBTypes[i] = 0;
    }
    old_timestamp = -1;
  }

  /*! Gets current run number 
    \return run number
   */
  int Reader::GetRunNo() { return runNos[fileIndx]; }
  /*! Gets current data file path
    \return path
   */
  std::string Reader::GetPath() { return runPaths[fileIndx]; }    

  /*! Gets current Dither value and increments the dither.
    \return dither
   */
  double Reader::GetDither() { 
    Reader::Dither += 0.02;
    if (Reader::Dither >= 0.9999) { Reader::Dither = 0.0; } 
    return Reader::Dither;
  }

  /*! Starts the clock for sorting */
  void Reader::Start() {
    start_time = std::chrono::system_clock::now();
  }

  /*! Stops the clock for sorting */
  void Reader::Stop() {
    stop_time = std::chrono::system_clock::now();
  }

  /*! Adds a processor
    \param type GEB type to associate the processor with
    \param proc Pointer to the processor object
   */
  void Reader::AddProcessor(int type, Processor *proc) {
    processors[type].push_back(proc);
  }

  /*! Main event loop function. Builds a physics event from succesive sub-events until the sequential time difference is greater than coinc_window.
    \return Error code. <br>
    0 - End of file reached, no new data in the event. <br>
    1 - No error, event built successfully.
  */
  int Reader::Next() {
    if (eof) { return 0; }
    if (nEvents > 0) {
      //copy the last header that was read into the first one for this event
      memcpy(&headers[0], &headers[subevt], sizeof(headers[subevt]));
    }
    
    subevt = 0;
    old_timestamp = -1;

    for (int geb=0; geb<MAX_GEB_TYPE; ++geb) {
      for (size_t i=0; i<processors[geb].size(); ++i) {        
        processors[geb][i]->Reset();
      }
    }
    
    nEvents += 1;
    int retval = 1;
    while (1) { //subevent loop
      if (nEvents == 1 && subevt == 0) {
        starttime = headers[subevt].timestamp * 10.0/(1e9*60.0);
        if (ts_mode == 1) { run_wt_offset = starttime; }
      }
      walltime = headers[subevt].timestamp * 10.0/(1e9*60.0);
        
      if (old_timestamp != -1) { 
        long long int dt = ((long long int)headers[subevt].timestamp - (long long int)old_timestamp)*10;
        if (dt < 0) {
          if (warning) {
            printf("SEVERE ERROR: dt = %lld, type = %i\n", dt, headers[subevt].type);
          }
          nOutOfOrder += 1;
        }
        if (dt > coinc_window) {
          break;           
        }  
      }

      if (false) {
        printf("subevt = %i\n", subevt);
        printf("type = %d\n", headers[subevt].type);
        printf("length = %d\n", headers[subevt].length);
        printf("timestamp = %ld\n", headers[subevt].timestamp);
      }
      
      //HeaderType htype(static_cast<HeaderType>(headers[subevt].type));

      //read the payload
      if (headers[subevt].length > MAX_GEB_PAYLOAD) { std::cerr << "Severe error! GEB payload size " << headers[subevt].length << " > " << MAX_GEB_PAYLOAD << ". Increase MAX_GEB_PAYLOAD" << std::endl; exit(1); };
      if (headers[subevt].type > MAX_GEB_TYPE) { std::cerr << "Severe error! GEB type " << headers[subevt].type << " encountered. Incrase MAX_GEB_TYPE" << std::endl; exit(1); }
      if (fCompressed) {
        if (gzread(gDataFile, sub[subevt], headers[subevt].length) == 0) {
          break;
        }
      }
      else {
        if (fread(sub[subevt], headers[subevt].length, 1, DataFile) == 0) {
          break;
        }
      }

      nGEBTypes[headers[subevt].type]++;
      //process the payload
      for (size_t i=0; i<processors[headers[subevt].type].size(); ++i) {
        if (processors[headers[subevt].type][i] != NULL) {
          processors[headers[subevt].type][i]->Process(headers[subevt].timestamp, (unsigned short int*)(&sub[subevt][0]), headers[subevt].length);
        }
        else if (warning) {
          std::cout << "Other event type = " << headers[subevt].type << ", length=" << headers[subevt].length << std::endl;          
        }        
      }
      
      old_timestamp = headers[subevt].timestamp;
      ++subevt;
      if (fCompressed) {
        if (gzread(gDataFile, &headers[subevt], sizeof(headers[subevt])) == 0) {
          eof = true;
          break;
        }
      }
      else {
        if (fread(&headers[subevt], sizeof(headers[subevt]), 1, DataFile) == 0) {
          eof = true;
          break;
        }
      }
      if (subevt >= MAX_GEB_SUBEVTS-1) { std::cerr << "Severe error! Number of sub-events > " << MAX_GEB_SUBEVTS << ". Increase MAX_GEB_SUBEVTS" << std::endl; break; }
    }  //end subevent loop

    if (subevt == 0) { eof = true; return 0; }  //no subevents, end of file

    for (int geb=0; geb<MAX_GEB_TYPE; ++geb) {
      for (size_t i=0; i<processors[geb].size(); ++i) {        
        processors[geb][i]->ProcessFinal(); //end-of-event processing
      }
    }
    
    return retval;
  }

  /*! Write to pruned output file
    \return Number of bytes written
  */
  int Reader::Write() {
    if (!PrunedOutput) { return -1; }
    int retval = 0;
    for (int i=0; i<subevt; ++i) {
      if (fCompressed) {
        gzwrite(gPrunedFile, &headers[i], sizeof(headers[i]));
        gzwrite(gPrunedFile, &(sub[i][0]), headers[i].length);
      }
      else {
        fwrite(&headers[i], sizeof(headers[i]), 1, PrunedFile);
        fwrite(&(sub[i][0]), headers[i].length, 1, PrunedFile);
      }
      retval += sizeof(headers[i]) + headers[i].length;
    }
    return retval;
  }

  /*! Print summary of sorting, usually at the end of each file. This function calls Processor::PrintSummary for each processor loaded.
   */
  int Reader::PrintSummary(std::ostream &out) {
    double duration = (double)(std::chrono::duration_cast <std::chrono::microseconds> (stop_time - start_time).count());
    
    out << "================= " << ANSI_COLOR_GREEN << nEvents << ANSI_COLOR_RESET << " events sorted in " << ANSI_COLOR_GREEN << (int)duration/1000000 << " s " << ANSI_COLOR_RESET << "================= " << std::endl;
    out << "   " << ANSI_COLOR_RED << nOutOfOrder << ANSI_COLOR_RESET << " out of time-order" << std::endl;
    
    for (int geb=0; geb<MAX_GEB_TYPE; ++geb) {
      if (nGEBTypes[geb] > 0) {
        out << "   Type " << geb << ": " << ANSI_COLOR_GREEN << nGEBTypes[geb] << ANSI_COLOR_RESET << " sub-events" << std::endl;
      }
    }
    
    out << "   Run was " << ANSI_COLOR_GREEN << walltime-starttime << ANSI_COLOR_RESET <<" min long from timestamps" << std::endl;
    out << std::setprecision(6);

    for (int geb=0; geb<MAX_GEB_TYPE; ++geb) {
      for (size_t i=0; i<processors[geb].size(); ++i) {
        processors[geb][i]->PrintSummary(out);
      }
    }

    /*
    out << "   " << ANSI_COLOR_YELLOW << nValidGretina << ANSI_COLOR_RESET << "/" << ANSI_COLOR_GREEN << nGretinaHits << ANSI_COLOR_RESET
        << " valid Gretina (" << ANSI_COLOR_YELLOW << std::setprecision(4) << validGret << "%" << ANSI_COLOR_RESET << ")"<< std::endl;
    out << "   " << ANSI_COLOR_YELLOW <<  nORRUBAHits-nSpuriousMyRIAD << ANSI_COLOR_RESET << "/" << ANSI_COLOR_GREEN << nORRUBAHits << ANSI_COLOR_RESET
        << " valid ORRRUBA (" << ANSI_COLOR_YELLOW << std::setprecision(4) << validORRUBA << "%" << ANSI_COLOR_RESET << ")  ";
    out << ANSI_COLOR_RED << nSpuriousMyRIAD << ANSI_COLOR_RESET << " bad MyRIAD timestamps" << std::endl;
    */


    return 0;
  }


  void BasicProcessor::Reset() { fired = false; } ;
  void BasicProcessor::Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length){
    time = timestamp;
    fired = true;
  }
  void BasicProcessor::ProcessFinal() { };
  void BasicProcessor::PrintSummary(std::ostream &out) {};
}
