#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


extern double MAX_EXP;
extern double eps;

extern std::vector<double> pars;
extern std::vector<double> Setting;
extern std::vector<double> Parameter;

///If traits are in ISF file give TRAITS 1, else 0 and specifiy in cvf file
#define TRAITS 0


///All parameter and Setting definitions

#define rho Parameter[0]
#define Xmax1 Parameter[1]
#define Xmax2 Parameter[2]
#define xb Parameter[3]
#define x_min Parameter[4]
#define x_j Parameter[5]
#define x_A Parameter[6]
#define W_0 Parameter[7]
#define A_max Parameter[8]
#define A_min Parameter[9]
#define psi_L Parameter[10]
#define Meta Parameter[11]
#define sigma Parameter[12]
#define Alpha Parameter[13]
#define dzeta_1 Parameter[14]
#define dzeta_2 Parameter[15]
#define dzeta_3 Parameter[16]
#define dzeta_4 Parameter[17]
#define p_1 Parameter[18]
#define p_2 Parameter[19]
#define k_e Parameter[20]
#define qj Parameter[21]
#define qa Parameter[22]
#define eta Parameter[23]
#define mu Parameter[24]
#define meta_qs Parameter[25]
#define Meta_mu Parameter[26]
#define starv_qs Parameter[27]
#define mu_s Parameter[28]
#define Liters Parameter[29]
#define delta_t Parameter[30]
#define mut_rate_psi Parameter[31]
#define mut_rate_meta Parameter[32]
#define mut_rate_xj Parameter[33]
#define mut_rate_xb Parameter[34]
#define Mut_var Parameter[35]
#define Sr Parameter [36]

#define Total_time Setting[0]
#define Output_time Setting[1]
#define Pop_Output Setting[2]
