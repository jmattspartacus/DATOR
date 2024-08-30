#include "S800Physics.hh"

namespace S800 {
  void PhysicsEvent::Process(unsigned long long int ts, unsigned short int *data, unsigned short int length) {
    s800_physicsdata phys_data;
    if (sizeof(phys_data) != length) {
      std::cerr << "Sever Error! S800 Physics data has inconsistent length" << std::endl;
      return;
    }

    timestamp = ts;

    s800_physicsdata *phys = (s800_physicsdata*)data;
    ata = phys->ata;
    bta = phys->bta;
    dta = phys->dta;
    yta = phys->yta;

    double x = std::sin(ata);
    double y = -std::sin(bta);
    double z = std::sqrt(1.0-x*x-y*y);
    pTheta = std::acos(z);
    pPhi = std::abs(y)/y * std::acos(x/std::sqrt(x*x + y*y));
    valid = true;
  }

  void PhysicsEvent::ProcessFinal() {}

  void PhysicsEvent::Reset() {
    valid = false;
  }
}
