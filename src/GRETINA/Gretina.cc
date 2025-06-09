#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "GRETINA/Gretina.hh"
#include "GRETINA/GretinaHit.hh"
#include "GRETINA/GretinaConf.hh"
#include "GRETINA/Gamma.hh"

namespace Gret {
  Configuration Event::conf;

  void Event::Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length) {
    crys_intpts crystal;
    if (sizeof(crystal) != length) {
      printf("Severe Error: Inconsistent Gretina Payload length\n");
      return;
    }

    crys_intpts *crys = (crys_intpts*)data;
    AddHit(timestamp, crys);
  }

  void Event::PrintSummary(std::ostream &out) {
    out << "--------- Gretina Summary --------" << std::endl;
    out << "   " << ANSI_COLOR_YELLOW << nValidGretina << ANSI_COLOR_RESET << "/" << ANSI_COLOR_GREEN << nGretina << ANSI_COLOR_RESET
      << " valid Gretina (" << ANSI_COLOR_YELLOW << std::setprecision(4) << (float)nValidGretina/(float)nGretina*100.0 << "%" << ANSI_COLOR_RESET << ")"<< std::endl;
    if (nBadT0 > 0) {
      out << "      " << ANSI_COLOR_YELLOW << nBadT0 << ANSI_COLOR_RESET
        << " bad t0 (" << ANSI_COLOR_YELLOW << std::setprecision(4) << (float)nBadT0/(float)nGretina*100.0 << "%" << ANSI_COLOR_RESET << ")"<< std::endl;
    }
    if (nBadIntE > 0) {
      out << "      " << ANSI_COLOR_RED << nBadIntE << ANSI_COLOR_RESET
        << " bad interaction energy (" << ANSI_COLOR_RED << std::setprecision(4) << (float)nBadIntE/(float)nGretina*100.0 << "%" << ANSI_COLOR_RESET << ")"<< std::endl;
    }
    if (nBadPAD > 0) {
      out << "      " << ANSI_COLOR_RED << nBadPAD << ANSI_COLOR_RESET 
        << " bad PAD (" << ANSI_COLOR_RED << std::setprecision(4) << (float)nBadPAD/(float)nGretina*100.0 << "%" << ANSI_COLOR_RESET << ")"<< std::endl;
    }
    out << "   " << ANSI_COLOR_YELLOW << nTotalGammas << ANSI_COLOR_RESET << " total gammas, " << ANSI_COLOR_GREEN << nEvents << ANSI_COLOR_RESET << " events" << std::endl;
    out <<"  average gamma-ray multiplicity " << ANSI_COLOR_GREEN << (float)nTotalGammas/(float)nEvents << ANSI_COLOR_RESET << std::endl;
  }

  int Event::AddHit(const int64_t GEBtimestamp,
      const crys_intpts *data) {

    if (nhits == GRETINA_MAXHITS) { std::cerr << "Warning! More than " << nhits << " hits in a single event!" << std::endl; return -1; }
    hits.emplace_back();
    nGretina += 1;
    if (hits[nhits].Build(GEBtimestamp, data) >= 0) {
      if (hits[nhits].valid) {
        nValidGretina += 1;
      }
      if (hits[nhits].BadIntE) {
        nBadIntE += 1;
      }
      if (hits[nhits].BadT0) {
        nBadT0 += 1;
      }
      if (hits[nhits].PAD > 0) {
        nBadPAD += 1;
      }
      nhits += 1;
      return 0;
    }
    nhits += 1;
    return -1;
  }

  int Event::BuildGammas() {    
    int nAddbacks = 0;
    std::vector<int> addbacks(nhits);  
    if (conf.Addback == AddbackType::Cluster2Xtl || conf.Addback == AddbackType::Nearest2Xtl ) {
      for (int i=0; i<nhits; ++i) {
        auto hit_i = hits[i];
        if (hit_i.valid == false) { continue; }
        if (conf.AddbackDets[hit_i.CrystalID] == 0) { continue; }
        for (int j=i+1; j<nhits; ++j) {
          auto hit_j = hits[j];
          if (hit_j.valid == false) { continue; }
          if (conf.AddbackDets[hit_j.CrystalID] == 0) { continue; }
          if (addbacks[i] == 1 || addbacks[j] == 1) { continue; }
          if (std::abs(hit_i.Time - hit_j.Time) > conf.AddbackTDiff) { continue; }
          bool inside = false;

          if (conf.Addback == AddbackType::Cluster2Xtl) {
            double costhet = std::sin(hit_i.Theta)*std::sin(hit_j.Theta)*std::cos(hit_i.Phi - hit_j.Phi) + std::cos(hit_i.Theta)*std::cos(hit_j.Theta);
            if (costhet > 1.0) { continue; }
            if (costhet >= std::cos(conf.ClusterAngle*3.14159265/180.0)) { inside = true; }
          }
          else if (conf.Addback == AddbackType::Nearest2Xtl) {
            //crystal center positions
            double costhet = std::sin(conf.CrystalTheta[hit_i.CrystalID])*std::sin(conf.CrystalTheta[hit_j.CrystalID])*std::cos(conf.CrystalPhi[hit_i.CrystalID] - conf.CrystalPhi[hit_j.CrystalID]) + std::cos(conf.CrystalTheta[hit_i.CrystalID])*std::cos(conf.CrystalTheta[hit_j.CrystalID]);
            if (costhet > 1.0) { continue; }
            if (costhet >= std::cos(25.0*3.14159265/180.0)) { inside = true; }
          } 

          if (!inside) { continue; }

          ++nAddbacks;

          addbacks[i] = 1;
          addbacks[j] = 1;

          if (hit_i.TotalEnergy > hit_j.TotalEnergy) {
            gammas.emplace_back(hit_i, j);
            gammas[ngammas].AddHit(hit_j, j);
          }
          else {
            gammas.emplace_back(hit_j, j);
            gammas[ngammas].AddHit(hit_i, i);
          }
          ngammas +=1;
        } //loop j
      } //loop i
    } //2 Xtl logic
    if (conf.Addback == AddbackType::ClusterAnyXtl || conf.Addback == AddbackType::NearestAnyXtl) {
      std::sort(hits.begin(), hits.end());
      for (int i=nhits-1; i>=0; i--) { //go through hits in descending order of energy
        auto hit_i = hits[i];
        bool found = false;
        if (hit_i.valid == false) { continue; }
        for (int j=0; j<ngammas; ++j) { //search through list of gammas and check if this hit is within cone
          if (conf.AddbackDets[hit_i.CrystalID] == 0) { break; }
          auto &gam_j = gammas[j];
          if (conf.AddbackDets[gam_j.ID] == 0) { continue; }
          if (std::abs(hit_i.Time - gam_j.Time) > conf.AddbackTDiff) { continue; }
          bool inside = false;
          if (conf.Addback == AddbackType::ClusterAnyXtl) {
            double costhet = std::sin(hit_i.Theta)*std::sin(gam_j.Theta)*std::cos(hit_i.Phi - gam_j.Phi) + std::cos(hit_i.Theta)*std::cos(gam_j.Theta);
            if (costhet < std::cos(conf.ClusterAngle*3.14159265/180.0) || costhet > 1.0) { continue; }
            inside = true; 
          }
          else if (conf.Addback == AddbackType::NearestAnyXtl) {
            double costhet = std::sin(conf.CrystalTheta[hit_i.CrystalID])*std::sin(conf.CrystalTheta[gam_j.ID])*std::cos(conf.CrystalPhi[hit_i.CrystalID] - conf.CrystalPhi[gam_j.ID]) + std::cos(conf.CrystalTheta[hit_i.CrystalID])*std::cos(conf.CrystalTheta[gam_j.ID]);
            if (costhet < std::cos(25.0*3.14159265/180.0) || costhet > 1.0) { continue; }
            inside = true; 
          }
          if (!inside) { continue; }

          //it is in the cone, add to hit
          nAddbacks += 1;
          gam_j.AddHit(hit_i, i);
          found = true;
          break;
        } //loop over existing gammas
        if (found) { continue; }
        else {
          gammas.emplace_back(hit_i, i);
          ngammas += 1;
        }
      } //loop through hits i
    } //any Xtl logic

    if (conf.Addback == AddbackType::Cluster2Xtl || conf.Addback == AddbackType::Nearest2Xtl || conf.Addback == AddbackType::NoAddback) {
      for (int i=0; i<nhits; ++i) {
        if (addbacks[i] == 0) {
          auto hit_i = hits[i];
          if (hit_i.valid == false) { continue; }
          gammas.emplace_back(hit_i, i);
          ngammas += 1;
        }
      }
    }

    nTotalGammas += ngammas;
    return nAddbacks; 
  }  
} 
