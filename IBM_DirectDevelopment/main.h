#include "Individual.h"
#include "Resource.h"
#include "Pars.h"


#include <vector>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <chrono>


#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>

using namespace std;
stringstream ss;

//Max exponent value
double MAX_EXP = 50;
//When is something zero?
double eps = 10e-9;


//Vectors for time-settings and parameters
std::vector<double> Setting(0);
std::vector<double> Parameter(0);


//Function to split a line. Needed to read the cvf file
  void split(const std::string &s, char delim, std::vector<std::string> &elems){
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);}
  }


//Mark dead individuals
bool IsMarkedToDelete(const Individual & o);

bool IsMarkedToDelete(const Individual & o)
{
  return o.Is_dead > 0;
}



////NOT NECESSARY BUT MAYBE GOOD TO STORE?
  //Function to sort list of pointers to objects
/*  template <typename T> bool PComp(const T * const & a, const T * const & b)
  {
     return *a < *b;
  }

  //Function to sort list of unique pointers to objects
  bool compare_by_uniqptr(const unique_ptr<Individual>& a,
                          const unique_ptr<Individual>& b) {
      return a->Age < b->Age;
  }

bool DEAD(const Individual* p);
bool DEAD(const Individual* p) {
  return p->Is_dead > 0;
 }

 bool UDEAD(const unique_ptr<Individual>& p);
 bool UDEAD(const unique_ptr<Individual>& p) {
   return p->Is_dead > 0;
 }*/
