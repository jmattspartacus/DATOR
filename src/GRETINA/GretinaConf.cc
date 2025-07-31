#include <cmath>
#include <sstream>
#include <fstream>

#include "GRETINA/Gretina.hh"
#include "GRETINA/GretinaHit.hh"
#include "GRETINA/GretinaConf.hh"

namespace Gret {
  Array::Array() { dim = -1; }

  Array::Array(int d) {
    dim = d;
    data = new double[d*d];
  }

  Array::Array(Array& other) {
    Allocate(other.dim);
    for (int i=0; i<dim; ++i) {
      for (int j=0; j<dim; ++j) {
        data[i*dim + j] = other.Get(i,j);
      }
    }
  }
  
  Array::~Array() {
    if (dim > 0) {
      delete data;
    }
  }
 
  void Array::Allocate(int d) {
    dim = d;
    data = new double[d*d];
  }

  void Array::Set(int i, int j, double val) {
    data[i*dim + j] = val;
  }

  double Array::Get(int i, int j) const {
    return data[i*dim + j];
  }
  
  void Array::operator=(Array &other) {
    if (dim > 0) {
      delete data;
    }
    Allocate(other.dim);
    for (int i=0; i<dim; ++i) {
      for (int j=0; j<dim; ++j) {
        data[i*dim + j] = other.Get(i,j);
      }
    }
  }
  
  int Configuration::ReadCrystalMap(std::string fn) {
    std::ifstream file(fn.c_str());
    std::string line;
    int prog = 0;
    int hole = 0;
    int xtl = 0;
    while (getline(file, line)) {
      //std::cout << line << std::endl;
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      std::stringstream ss(line);
      if (prog == 0) {
        ss >> hole >> xtl;
        CrystalMap[4*(hole-1)+xtl].Allocate(4);
        //std::cout << hole << "   " << xtl << std::endl;
        prog = 1;
      }
      else if (prog > 0 ) {
        double a, b, c, d;
        ss >> a >> b >> c >> d;
        /*
        CrystalMap[{hole, xtl}][prog - 1][0] = a;
        CrystalMap[{hole, xtl}][prog - 1][1] = b;
        CrystalMap[{hole, xtl}][prog - 1][2] = c;
        CrystalMap[{hole, xtl}][prog - 1][3] = d;
        */
        CrystalMap[4*(hole-1)+xtl].Set(prog - 1, 0, a);
        CrystalMap[4*(hole-1)+xtl].Set(prog - 1, 1, b);
        CrystalMap[4*(hole-1)+xtl].Set(prog - 1, 2, c);
        CrystalMap[4*(hole-1)+xtl].Set(prog - 1, 3, d);
        ++prog;
      }
      if (prog == 5) {
        /*
          std::cout << hole << "   " << xtl << "   " << CrystalMap[{hole,xtl}][0][0] << "   "
          << CrystalMap[{hole,xtl}][0][1] << "   "
          << CrystalMap[{hole,xtl}][0][2] << "   "
          << CrystalMap[{hole,xtl}][0][3] <<  std::endl;
          std::cout << hole << "   " << xtl << "   " << CrystalMap[{hole,xtl}][1][0] << "   "
          << CrystalMap[{hole,xtl}][1][1] << "   "
          << CrystalMap[{hole,xtl}][1][2] << "   "
          << CrystalMap[{hole,xtl}][1][3] <<  std::endl;
          std::cout << hole << "   " << xtl << "   " << CrystalMap[{hole,xtl}][2][0] << "   "
          << CrystalMap[{hole,xtl}][2][1] << "   "
          << CrystalMap[{hole,xtl}][2][2] << "   "
          << CrystalMap[{hole,xtl}][2][3] <<  std::endl;
          std::cout << hole << "   " << xtl << "   " << CrystalMap[{hole,xtl}][3][0] << "   "
          << CrystalMap[{hole,xtl}][3][1] << "   "
          << CrystalMap[{hole,xtl}][3][2] << "   "
          << CrystalMap[{hole,xtl}][3][3] <<  std::endl;
        */
        prog = 0;
      }
    }

    //pre-calculate crystal theta/phi coordinates for nearest-neighbour addback
    for (int i=4; i<124; ++i) {
      auto pos = GetPosition(i/4, i%4, 0.,0.,0.);
      double PosX = pos[0];
      double PosY = pos[1];
      double PosZ = pos[2];

      double x2y2 = std::sqrt(PosX*PosX + PosY*PosY);
      double x2y2z2 = std::sqrt(PosX*PosX + PosY*PosY + PosZ*PosZ);
      double costhet = PosZ/x2y2z2;
      if (costhet >= -1.0 and costhet <= 1.0) {
        CrystalTheta[i] = std::acos(costhet);
        //std::cout << "Crystal = " << CrystalID << "   Theta = " << Theta*180.0/3.1415926535 << std::endl;
      }
      else {
        std::cout << "Cos(theta) was out of range: " << costhet << std::endl;
        std::cout << PosX << "   " << PosY << "   " << PosZ << std::endl;
        CrystalTheta[i] = 0;
      }
      double cosphi = PosX/x2y2;
      if (cosphi >= -1.0 and cosphi <= 1.0) {
        CrystalPhi[i] = std::acos(cosphi);
        if (PosY < 0) {
          CrystalPhi[i] = 2.0*3.1415926536 - CrystalPhi[i];
        }
      }
      else {
        std::cout << "Cos(phi) was out of range: " << cosphi << std::endl;
        CrystalPhi[i] = 0;
      }
    }

    return 0;
  }
  
  std::vector<double> Configuration::GetPosition(int hole, int xtl, double x, double y, double z) const {
    /*
    const auto coords = CrystalMap.at({hole,xtl});
double ret_x = (coords[0][0] * x) + (coords[0][1] * y) + (coords[0][2] * z) + coords[0][3];
    double ret_y = (coords[1][0] * x) + (coords[1][1] * y) + (coords[1][2] * z) + coords[1][3];
    double ret_z = (coords[2][0] * x) + (coords[2][1] * y) + (coords[2][2] * z) + coords[2][3];
    */

    const Array* coords = &CrystalMap[4*(hole-1)+xtl];

    double ret_x = (coords->Get(0,0) * x) + (coords->Get(0,1) * y) + (coords->Get(0,2) * z) + coords->Get(0,3);
    double ret_y = (coords->Get(1,0) * x) + (coords->Get(1,1) * y) + (coords->Get(1,2) * z) + coords->Get(1,3);
    double ret_z = (coords->Get(2,0) * x) + (coords->Get(2,1) * y) + (coords->Get(2,2) * z) + coords->Get(2,3);

    if (coords->Get(1,3) > 0) { //beam left
      ret_y += BeamLeftOffset;
      ret_z += BeamLeftZOffset;
    }
    if (coords->Get(1,3) < 0) { //beam right
      ret_y -= BeamRightOffset;
      ret_z += BeamRightZOffset;
    }
    
    ret_x -= PosXOffset;
    ret_y -= PosYOffset;
    ret_z -= PosZOffset;

    return {ret_x, ret_y, ret_z};
  }

  void Configuration::PrintAngles() {
    std::cout << "ID       Theta        Phi" << std::endl;
    for (int hole = 1; hole <= 30; ++hole) {
      for (int xtl = 0; xtl < 4; ++xtl) {
        auto pos = GetPosition(hole, xtl, 0,0,0);    
        double PosX = pos[0];
        double PosY = pos[1];
        double PosZ = pos[2];

        double x2y2 = std::sqrt(PosX*PosX + PosY*PosY);
        double x2y2z2 = std::sqrt(PosX*PosX + PosY*PosY + PosZ*PosZ);
        double costhet = PosZ/x2y2z2;
        double Theta = std::acos(costhet)*180.0/3.141592635;
        double cosphi = PosX/x2y2;
        double Phi = std::acos(cosphi);
        if (PosY < 0) {
          Phi = 2.0*3.1415926536 - Phi;
        }
        Phi = Phi*180.0/3.1415926535;
      
        int id = (hole-1)*4 + xtl;
        std::cout << id << "    " << Theta << "     " << Phi << std::endl;
      }
    }
  }

  int Configuration::ReadCalibration(std::string fn) {
    std::ifstream file(fn.c_str());
    std::cout << fn.c_str() << std::endl;
    int id;
    double gain, offset, quad;
    while (file >> id >> offset >> gain >> quad) {
      std::cout << id << "   " << offset << "  " << gain << std::endl;
      CalGain[id] = gain;
      CalOffset[id] = offset;
    }
    return 0;
  }

  int Configuration::ReadEfficiencyCal(std::string fn) {
    if ( EfficiencyType == EffType::kLogLogPoly ) {
      std::ifstream file(fn.c_str());
      std::cout << fn.c_str() << std::endl;
      float eff;
      int i=0;
      while (file >> eff) {
        std::cout << eff << std::endl;
        EfficiencyCal[i] = eff;
        ++i;
      }
      return 0;
    }    
    else if ( EfficiencyType == EffType::kRadWare ) {
      std::ifstream file(fn.c_str());
      std::string line;
      getline(file, line);
      getline(file, line);
      std::stringstream ss(line);
      float abseff, a, aerr, b, berr, c, cerr, d, derr, e, eerr, f, ferr, g, gerr;
      ss >> abseff >> a >> aerr >> b >> berr >> c >> cerr >> d >> derr >> e >> eerr >> f >> ferr >> g >> gerr;
      AbsEff = abseff;
      EfficiencyCal[0] = a;
      EfficiencyCal[1] = b;
      EfficiencyCal[2] = c;
      EfficiencyCal[3] = d;
      EfficiencyCal[4] = e;
      EfficiencyCal[5] = f;
      EfficiencyCal[6] = g;
      return 0;
    }
    std::cout << "Efficiency Calibration not read successfully!" << std::endl;
    return 1;
  }  

  float Configuration::Efficiency(const double &e) const {
    if ( EfficiencyType == EffType::kLogLogPoly ) {
      double loge = std::log10(e);
      float eff = EfficiencyCal[0] +
        EfficiencyCal[1]*loge +
        EfficiencyCal[2]*std::pow(loge,2) +
        EfficiencyCal[3]*std::pow(loge,3) +
        EfficiencyCal[4]*std::pow(loge,4);
      return std::pow(10.0, eff);
    }
    else if (EfficiencyType == EffType::kRadWare) {
      double x = std::log(e/100.0);
      double y = std::log(e/1000.0);
      float eff = std::pow((EfficiencyCal[0] + EfficiencyCal[1]*x + EfficiencyCal[2]*x*x), -EfficiencyCal[6]) +
        std::pow((EfficiencyCal[3] + EfficiencyCal[4]*y + EfficiencyCal[5]*y*y), -EfficiencyCal[6]);
      eff = std::pow(eff, -1.0/EfficiencyCal[6]);
      eff = AbsEff*std::exp(eff);
      return eff;
    }

    return 0;
  }

  int Configuration::ReadCrystalSwap(std::string fn) {
    std::ifstream file(fn.c_str());
    int crst_i, crst_j;
    while (file >> crst_i >> crst_j) {
      CrystalSwap[crst_i] = crst_j;      
    }
    file.close();
    return 0;
  }

  int Configuration::ReadNoisyCrystals(std::string fn) {
    std::ifstream file(fn.c_str());
    int i, noisy;
    while (file >> i >> noisy) {
      NoisyCrystals[i] = noisy;
    }
    file.close();
    return 0;
  }

  int Configuration::ReadAddbackDets(std::string fn) {
    std::ifstream file(fn.c_str());
    int i, addback;
    while (file >> i >> addback) {
      AddbackDets[i] = addback;
    }
    file.close();
    return 0;
  }
}
