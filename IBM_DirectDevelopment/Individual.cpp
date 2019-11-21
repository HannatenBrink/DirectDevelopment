
#include "Individual.h"


//member function definitions including Constructor
Individual::Individual(double age, double Irr, double Rev, double has_meta, double meta_traits, double psi_traits, double xj_traits, double xb_traits){
  x = Irr;
  y = Rev;
  w = x + y;
  Age = age;
  //If individuals are already old in the initial run, they get trait size, otherwise the real size at birth. This is simplest way and shouldn't be a big problem
  if (Age > 0) {
    SizeBorn = xb_trait;
  } else {
    SizeBorn = Irr;
  }
  Repro_buf = 0;
  Offspring = 0;

  Attack1 = 0;
  Attack2 = 0;
  Handling = 0;

  Intake1 = 0;
  Intake2 = 0;

  Is_dead = 0;
  starve = 0;
  Has_meta = has_meta;


  meta_trait = meta_traits;
  psi_trait = psi_traits;
  xj_trait = xj_traits;
  xb_trait = xb_traits;

  delta_x = 0;
  delta_y = 0;
  delta_buf = 0;
  Total_buf = 0;

  if (Irr >= xj_trait) {
    Has_meta = 1;
  }


}

Individual::~Individual() {}

//For sorting of the vector
bool Individual::operator <(Individual const& IndividualObj)const
	{
		return Age < IndividualObj.Age;
	}

void Individual:: printInfo() {
  std::cout << "Age:    " << Age << std::endl;
  std::cout << "Irr: " << x << std::endl;
  std::cout << "Rev: " << y << std::endl;
  std::cout << "Intake 1: " << Intake1 << std::endl;
  std::cout << "Intake 2: " << Intake2 << std::endl;
  std::cout << "Attack 1: " << Attack1 << std::endl;
  std::cout << "Attack 2: " << Attack2 << std::endl;
  std::cout << "Standard weight: " << w << std::endl;
  std::cout << "Repro buffer: " << Repro_buf << std::endl;
  std::cout << "Offspring: " << Offspring << std::endl;
  std::cout << "Is it dead? " << Is_dead << std::endl;
  std::cout << "Has it metamorphosed? " << Has_meta << std::endl << std::endl;
  std::cout << "Trait values (meta, psi, xj, Xb): " << meta_trait <<", " << psi_trait << ", " <<xj_trait << ", " << xb_trait<< std::endl << std::endl;
}

void Individual::R_Intake(double R_dens1, double R_dens2) {

  double psi_A = std::min(1.0, psi_trait + meta_trait); //Specialization parameter after metamorphosis
  double psi;
  double W_min = x_min * (1 + qj);
  /*
      if (Has_meta)  //metamorphosed individuals
          psi = psi_A;
      else
          psi = psi_trait;*/

//Switch is a bit faster than if else.
switch(Has_meta) {
  case 0:
  psi = psi_trait;
  break;
  case 1:
  psi = psi_A;
  break;
  default:
  std::cerr<<"This individual has no information about metamorphosis\n";
  exit(1);
}

      //Maximum Attack rates
      double A1 = (1 - psi) * (A_max - A_min) + A_min;
      double A2 = psi * (A_max - A_min) + A_min;
      //std::cout << "In indiv" <<std::endl;
      //Handling time and Attack rates
      Handling = dzeta_1 + dzeta_2 * pow(w, -dzeta_3)*exp(w * dzeta_4);
      Attack1 = A1 * pow((w / W_0) * exp((1 - (w / W_0))), Alpha);
      //Attack1 = A1 * pow(((w / W_0) * exp(1 - (w / W_0))), Alpha);
      if (x <= x_min) {
          Attack2 = 0;
        }
      else {
          Attack2 = A2 * pow(((w - W_min) / W_0 ) * exp((1 - (w - W_min) / W_0)), Alpha);
        }
      //Attack2 = A1;
      //Preference for R1
      double pw = sigma * (Attack2 * R_dens2 - Attack1 * R_dens1);
      pw = std::max(pw, -MAX_EXP);
      pw = std::min(pw, MAX_EXP);
      pw = exp(pw);
      double theta = 1 / (1 + pw);

      //Food intake
      double denom = 1 + Handling * ((theta * Attack1 * R_dens1) + ((1 - theta) * Attack2 * R_dens2));
      Intake1 = std::max(0.0, theta * Attack1 * R_dens1 / denom);
      Intake2 = std::max(0.0, (1 - theta) * Attack2 * R_dens2 / denom);
      //Attack2 = theta;
      //Metabolic costs and biomass production
      double metabolic = p_1 * pow((x + y), p_2);
      double biomassproduction = k_e * (Intake1 + Intake2) - metabolic;

      //Allocation rules
      double kappa_j = y / ((1 + qj) * qj * x);
      double kappa_a = y / ((1 + qa) * qa * x);
      double kappa_r;



      if (y < x * qj)
          kappa_r = 1 - kappa_a;
      else
          kappa_r = (1 - kappa_j) * kappa_a / kappa_j;

      ////////growth or fat burning////////////

      if (biomassproduction > 0){
          starve = 0;
          if (x < x_A){
              delta_x = kappa_j * biomassproduction;
              delta_y = (1 - kappa_j) * biomassproduction;
            }
          else{
              delta_x = kappa_a * biomassproduction;
              delta_y = kappa_r * biomassproduction;
            }
      }
      else {
          starve = 1;
          delta_x = 0;
          delta_y = std::max(biomassproduction, - y);
        }

      ////////Reproduction//////////////
      if ((y >= x * qj) && (x >= x_A)) {
          if (xb_trait < xj_trait) {
          delta_buf = std::max(0.0, ((1 - kappa_a / kappa_j) * eta * biomassproduction / ((1 + qj) * xb_trait)));
        }
        else {
          delta_buf = std::max(0.0, ((1 - Meta_mu*meta_trait)*(1 - kappa_a / kappa_j) * eta * biomassproduction / ((1 + qj) * xb_trait)));
          //std::cout << "Xb < meta" << std::endl;
         } //Cost for meta if it is before birth.
       }
      else
          delta_buf = 0;
}

void Individual::Grow(){
  Age += delta_t;
  x += delta_x * delta_t;
  y += delta_y * delta_t;
  Repro_buf += delta_buf * delta_t;
  Total_buf += delta_buf * delta_t;
  w = x * (1 + qj);
}

void Individual:: Die(double Rand){
  double prob = Rand; //Random number
  double mu_starv = 0; //starvation mort
  if (y < eps) //If y is very small, die immediately.
    mu_starv = 1 / delta_t;
  else if  (y / x <= starv_qs) //If individuals have lost too much fat, they experience starvation mortality
      mu_starv = mu_s * (starv_qs * x / y - 1);
  double mu_total = (mu + mu_starv);
  if (mu_total * delta_t > prob){
     Is_dead = 1;}
}


void Individual:: Repro(){
     Offspring = std::trunc(Repro_buf/1);
     Repro_buf = std::fmod(Repro_buf, 1);
}

void Individual:: Metamorphose(double Rand) {
  y = y - meta_trait * std::max(0.0,(xj_trait - SizeBorn)) * (qj - meta_qs); ///Reduction in reversible mass. Mortality because of mass loss only later
  double mort_risk = meta_trait * Meta_mu; //Mortality risk
  double prob = Rand; //Random number
  //std::cout << "risk: " <<mort_risk << " prob" << prob << std::endl;
  if (mort_risk > prob) {
    Is_dead = 1;
  }
  Has_meta = 1;
}
