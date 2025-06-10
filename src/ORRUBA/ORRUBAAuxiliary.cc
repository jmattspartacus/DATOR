#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include "ORRUBA/ORRUBAAuxiliary.hh"

namespace Orruba {
  Table::Table() {
    nrows = 0;
    ncols = 0;
  }

  Table::Table(int nr, int nc) {
    nrows = nr;
    ncols = nc;
    data = new double[nr*nc];    
  }

  Table::~Table() {
    if (nrows > 0 && ncols > 0) {
      delete data;
    }
  }

  void Table::Initialize(int nr, int nc) {    
    if (nrows > 0 && ncols > 0) {
      delete data;
    }
    nrows = nr;
    ncols = nc;
    names.reserve(nr);
    data = new double[nr*nc];
  }

  double Table::Get(int row, int col) {
    return data[row*ncols + col];
  }

  void Table::Set(int row, int col, double val) {
    if (col >= ncols || row >= nrows) { std::cerr << "Warning! Attempt to set value outside of range in ORRUBA Scalers! Make sure the table is initialized correctly" << std::endl; return; } 
    data[row*ncols + col] = val;
  }

  void Table::Print() {
    for (int r=0; r<nrows; ++r) {
      printf("%20s  ", names[r].c_str());
      for (int c=0; c<ncols; ++c) {
        if (c == 0 ) {
          printf("%8.2e   ", Get(r,c));
        }
        else {
          printf("%8.2f   ", Get(r,c));
        }
      }
      printf("\n");
    }
  }
  
  void Auxiliary::Process(unsigned long long int ts,
                          unsigned short int *data,
                          unsigned short int length) {
    timestamp = ts;
    nRows = 0;
    
    int charCt = 0;
    for (int i=0; i<length; i=i+80) {
      nRows += 1;
      std::memcpy(&rows[nRows-1][0], &data[40*(nRows-1)], 80);
      rows[nRows-1][80]='\0';
    }

    header = rows[0];
    std::stringstream ss;
    ss.str(header);
    std::string dump;

    ss >> dump;
    ss >> dump;

    std::string date;
    std::string time;
    ss >> date;
    ss >> time;
    
    std::tm t;
    //make sure to let std::mktime determine if DST is in effect
    //note that this might cause some weirdness if the analysis is taking place in a different timezone (?)
    t.tm_isdst = -1; 
    /*
    if (sscanf(time.c_str(), "%2d:%2d:%2d",
           &t.tm_hour,
           &t.tm_min,
               &t.tm_sec) !=3) { std::cout << "ERROR in sscanf!" << std::endl; exit(1); };

    char month[3];
    int year;
    if (sscanf(date.c_str(), "%2d-%3s-%2d",
               &t.tm_mday,
               &month[0],
               &year
               ) != 3) { std::cout << "ERROR in sscanf!" << std::endl; exit(1); };

    std::vector<std::string> months = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
    for (int i=0; i<12; ++i) {
      if (!months[i].compare(month)) { t.tm_mon = i; break; }
    }
    if (year < 69) { year += 100; }
    t.tm_year = year; //since 1900;    

    */

    date = date+" "+time;
    ss.clear();
    ss.str(date);    
    ss >> std::get_time(&t, "%d-%b-%y%t%H:%M:%S");
    if (ss.fail()) { std::cerr << "Parsing of date failed!" << std::endl; }    
    time_since_epoch = std::mktime(&t);
                        
    if (scalers.nrows == 0 && scalers.ncols == 0) { return; }
    
    scalers.names.clear();

    for (int i=1; i<nRows; ++i) {
      ss.clear();
      std::string line(rows[i]);
      ss.str(line);
      std::string name;
      if (!(ss >> name)) { break; };
      scalers.names.push_back(name);
      std::string val;
      int col = 0;
      while (ss >> val) {
        int pos = val.find("D");
        val.replace(pos, 1, "E");
        
        scalers.Set(i-1, col, std::atof(val.c_str()));
        col++;
      }
    }
  }

  void Auxiliary::Print() {
    if (nRows == 0) { return; }
    std::cout << std::endl;
    if (scalers.nrows == 0 && scalers.ncols == 0) {
      for (int i=0; i<nRows; ++i) {
        if((int)rows[i][0] == 0) { continue; }
        if(rows[i][0] == ' ') { continue; }
        std::cout << rows[i] << std::endl;
      }
    }
    else {
      std::cout << header << std::endl;
      scalers.Print();
    }
  }
}
