#include "Resource.h"
//Resource is in DENSITY!

//member function definitions including Constructor
Resource::Resource(double density, double MaxDens, double Lit, std::string Rname) {
  Density = density;
  Kr = MaxDens;
  Liter = Lit;
  Name = Rname;
}

Resource::~Resource() {}

void Resource:: printInfo()
{
  std::cout << "Density of Resource " << Name <<": "<< Density << " Liters: " << Liter << std::endl << std::endl;
}

void Resource::Growth(double popIntake, double Rdep) {
  double Delta_R;
  Delta_R = rho * ((Kr - Rdep) - Density) - popIntake/Liter;
  //Delta_R = rho * ((Kr - Rdep) - Density); //For individualtest
  Density += Delta_R * delta_t;
  Density = std::max(0.0, Density); //make sure it does not drop below 0
  // std::cout << "Liters: " << Liter << std::endl;
  // std::cout << "rho: " << rho << std::endl;
  // std::cout << "Kr: "  << Kr << std::endl;
  // std::cout << "Popintake: " <<popIntake << std::endl;
  // std::cout << "Density: " <<Density << std::endl;
}
