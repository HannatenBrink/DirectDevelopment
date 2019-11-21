#include "Pars.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <cmath>


class Individual {

  public:
    double x; //Irreversible
    double y; //Reversible
    double w; //Weight
    double Age;
    double SizeBorn; //Size at which individual is born. This is not necessarily the same as the xb_trait. Since xb_trait determines offspring size of this individual, but it might itself have been born with a dfiferent size

    double Repro_buf; //How much energy for reproduction
    double Offspring; //How many offspring has it produced

    double Attack1; //Attack rate on R1
    double Attack2; //Attack rate on R2
    double Handling; //Handling time

    double Intake1;  //Intake of R1
    double Intake2;  //Intake of R2

    int Is_dead;  //Is it dead?
    int starve; //Is it starving?
    int Has_meta; //Has it metamorphosed?

    double delta_x;
    double delta_y;
    double delta_buf;
    double Total_buf;

    double meta_trait;
    double psi_trait;
    double xj_trait;
    double xb_trait;



  // Constructor for Individual
    Individual (double age, double Irr, double Rev, double has_meta, double meta_traits, double psi_traits, double xj_traits, double xb_traits);

  //Deconstructor for Individual
    ~Individual();

//

  // Member functions declarations
    void printInfo();
    bool operator <(Individual const & IndividualObj)const;
    void Grow();
    void R_Intake(double R_dens1, double R_dens2);
    void Die(double Rand);
    void Repro();
    bool   Is_death();
    void Metamorphose(double Rand);
};
