#include "Pars.h"
#include <vector>
#include <iostream>
#include <algorithm>


class Resource {

  public:
    double Density;
    double Kr; //Carrying capacity
    double Liter;
    std::string Name;


  //Constructur for a Resource
  Resource(double density, double MaxDens, double Lit, std::string Rname = "Rx");

  //Deconstruc
  ~Resource();

  //member function declar
  void printInfo();
  void Growth(double popIntake, double Rdep);
};
