/*
 * Onto_scaling_TPB1998.h
 *
 * Header file specifying the life history of roach  (Ontogenetic scaling of foraging rates and the dynamics of a size-structured consumer-resource model. Persson et al 1998)
 *
 *
 *
 * Last modification: Htb - Nov 12 2019
 *  In the ms, we use parameter X_J to refer to the irreversible body mass at metamorphosis. However, for numerical reasons, in the model formulation of the pspm we use X_J=Xj+Xmin
     To contintue the ERP, it is sometimes necessary to make use of an alternative formulation (Loss_of_Meta_cwl.h), where the irreversible body mass at metamorphosis X_J = Xj. In the figures, total body mass at birth is calculated as (Xj+Xmin)*1.742 in case we used the standard model (Loss_of_meta_pspm.h), and as Xj*1.742 in case of the alternative model (Loss_of_Meta_cwl)
 */

/*
 *===========================================================================
 * 		DEFINITION OF PROBLEM DIMENSIONS AND NUMERICAL SETTINGS
 *===========================================================================
 */
// Dimension settings: Required
#define POPULATION_NR		    1
#define STAGES              	4
#define	I_STATE_DIM         	3
#if (PSPMIND == 1)
#define ENVIRON_DIM         3
#define Birthrate(n)        E[1+(n)]
#else
#define	ENVIRON_DIM         	2
#endif
#define INTERACT_DIM	       24
#define	PARAMETER_NR	       27
#define COHORT_NR               1000
#define FULLSTATEOUTPUT         100


// Numerical settings: Optional (default values adopted otherwise)
#define MIN_SURVIVAL	    	1.0E-9 		// Survival at which individual is considered dead
#define MAX_AGE                 100000		// Give some absolute maximum for individual age


#define DYTOL               	1.0E-7          // Variable tolerance
#define RHSTOL              	1.0E-7     		// Function tolerance
#define ESSTOL              	1.0E-12



/*
 *===========================================================================
 * 		DEFINITION OF ALIASES
 *===========================================================================
 */
// Define aliases for the istate variables
#define AGE         istate[0][0]
#define IR          istate[0][1]    	// Irreversible mass
#define REV         istate[0][2]    	// Reversible mass

// Define aliases for the environmental variables
#define R1          E[0]
#define R2          E[1]

// PARAMETERS
#define delta       parameter[ 0]	// Default:
#define R1MAX	    parameter[ 1]	// Default:
#define R2MAX	    parameter[ 2]	// Default:

#define AMIN        parameter[ 3]	// Default:
#define AMAX        parameter[ 4]
#define ALPHA       parameter[ 5]   	// Default:
#define W0          parameter[ 6]	// Default:

#define CHI1        parameter[ 7]  	// Default:
#define CHI2        parameter[ 8]	// Default:
#define	CHI3        parameter[ 9]	// Default:
#define	CHI4        parameter[10]	// Default:

#define M1          parameter[11]	// Default:
#define M2          parameter[12]	// Default:

#define K1          parameter[13]	// Default:
#define K2          parameter[14]	// Default:

#define X0          parameter[15]	// Default:
#define XMIN        parameter[16]
#define XJ          parameter[17]
#define XF          parameter[18]	// Default:

#define QJ          parameter[19]	// Default:
#define QA          parameter[20]	// Default:
#define QS          parameter[21]	// Default:

#define MU0         parameter[22]   	// Default:
#define SIGMA       parameter[23]
#define PSI         parameter[24]
#define META        parameter[25]   	// Meta extend
#define Rho         parameter[26]   	// Metamorphosing mortality


/*
 *===========================================================================
 * 		DEFINITION OF NAMES AND DEFAULT VALUES OF THE PARAMETERS
 *===========================================================================
 */
// At least two parameters should be specified in this array
char  *parameternames[PARAMETER_NR] =
{"Delta","R1max","R2max","Amin","Amax","Alpha","W0","ch1","chi2","chi3","chi4","m1","m2","k1","k2","x0","xmin","Xj","Xf","Qj","Qa","Qs","Mu0","sigma","PSI","Meta","Rho"};

// These are the default parameters values
double	parameter[PARAMETER_NR] =
{0.1,5,15,10000,100000,0.6,17.42,4E-6,8.19e-5,0.68,1.15e-3,0.033,0.77,6.71e-6,0.5,0.000804,1,1.1,5,0.742,1,0.2,0.01,10,1,1,1};

/*
 *===========================================================================
 * 		DEFINITION OF THE LIFE HISTORY MODELS FOLLOWS BELOW
 *===========================================================================
 * Specify the number of states at birth for the individuals in all structured
 * populations in the problem in the vector BirthStates[].
 *===========================================================================
 */

void SetBirthStates(int BirthStates[POPULATION_NR], double E[])
{
  BirthStates[0] = 1;
  return;
}


/*
 *===========================================================================
 * Specify all the possible states at birth for all individuals in all
 * structured populations in the problem. BirthStateNr represents the index of
 * the state of birth to be specified. Each state at birth should be a single,
 * constant value for each i-state variable.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 *===========================================================================
 */

void StateAtBirth(double *istate[POPULATION_NR], int BirthStateNr, double E[])
{
  AGE = 0.0;
  IR  = X0;
  REV = QJ*X0;

  return;
}


/*
 *===========================================================================
 * Specify the threshold determining the end point of each discrete life
 * stage in individual life history as function of the i-state variables and
 * the individual's state at birth for all populations in every life stage.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 *===========================================================================
 */



void IntervalLimit(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
        double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
        double limit[POPULATION_NR])
{
  switch (lifestage[0])
  {
    case 0:
      limit[0] = IR - XMIN; ///Second resource available
      break;
    case 1:
      limit[0] = IR - XMIN - XJ; //Metamorphosis
      break;
    case 2:
      limit[0] = IR - XF; //Maturation
      break;
    case 3:
      limit[0] = AGE - MAX_AGE; //Maximun age 
      break;
  }

  return;
}




/*
 *===========================================================================
 * Specify the individual development of individuals as function of i-state
 * and environment for all individuals of all populations in every life stage
 *
 * Notice that the first index of the variables 'istate[][]' and 'growth[][]'
 * refers to the number of the structured population, the second index refers
 * to the number of the individual state variable. The interpretation of the
 * latter is up to the user.
 *===========================================================================
 */

void Development(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
        double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
        double development[POPULATION_NR][I_STATE_DIM])
{
  const double	  MAX_EXP = 50.0;
  double	  pw;
  double          w, wr, h,a1,a2, gamma, Eg, Em,Ei,kappa_qj,kappa_qa;
  double          In1,In2,wmin, PHI,A1, A2, PSI_A;

  PSI_A = min(1, PSI+META);                       //Specialization R1 J&A
  w = (1+QJ)*IR;                                  // Effective body mass
  wmin = (1+QJ)*XMIN;                               // Min Weight for R2 (nu wanneer larval stage is beendig)
  h = (CHI1 + CHI2*pow(w, -CHI3)*exp(CHI4*w));    //handling time

   if ((lifestage[0] == 0) || (lifestage[0] == 1))
    {
      A1 = (1-PSI)*(AMAX-AMIN)+AMIN;
      A2 = PSI*(AMAX-AMIN)+AMIN;
    }
  else
    {
      A1 = (1-PSI_A)*(AMAX-AMIN)+AMIN;
      A2 = PSI_A*(AMAX-AMIN)+AMIN;
      
    }


  wr = w/W0;
  a1 = A1*pow(wr, ALPHA) * exp(ALPHA*(1 - wr));

  if (lifestage[0] == 0)
    {
      a2 = 0;
     
    }
  else
    {
      wr  = (w-wmin)/W0;
      pw  = log(wr + epsMach) + (1-wr);
      pw *= ALPHA;
      pw  = max(pw, -MAX_EXP);
      pw  = min(pw,  MAX_EXP);
      a2  = A2*exp(pw);
    }

  pw = -SIGMA * (a1  * R1 - a2  * R2);
  pw = max(pw, -MAX_EXP);
  pw = min(pw,  MAX_EXP);
  pw = exp(pw);

  PHI = 1 / (1 + pw);                                  //Time on R1
  gamma = (a1*PHI*R1+a2*(1-PHI)*R2) / (1 + h*(a1*PHI*R1+a2*(1-PHI)*R2));      //Intake rate
  In1 =(a1*PHI*R1) / (1 + h*(a1*PHI*R1+a2*(1-PHI)*R2));                           //Intake rate R1
  In2 = (a2*(1-PHI)*R2) / (1 + h*(a1*PHI*R1+a2*(1-PHI)*R2));                      //Intake rate R2
  Em = M1*pow(IR+REV, M2);                                                                    //maintenance cost
  Ei = (K1)*gamma;                                                                              //Conversion efficiency
  Eg = Ei - Em;                                                                               // Net production
  kappa_qj = 1 / ((1+QJ)*QJ) * (REV / IR);                                                      //fraction allocated to irreversible mass juve
  kappa_qa = 1 / ((1+QA)*QA) * (REV / IR);

  development[0][0] = 1.0;

  if (lifestage[0] < 3)
    {
      development[0][1] = kappa_qj*Eg;
      development[0][2] = (1-kappa_qj)*Eg;
    }
  else
    {
      double fecfrac = 1-(1+QJ)*QJ/((1+QA)*QA);    //1-(ka/kj)
      double s = (1000 * (REV / (QJ * IR)) - 999);
      s = (s < 0) ? 0 : s;
      s = (s > 1) ? 1 : s;
      fecfrac *= 3.0*s*s - 2.0*s*s*s;

      development[0][1] = kappa_qa*Eg;        			//ir
      development[0][2] = (1 - kappa_qa - fecfrac)*Eg;  	//rev
    }

  return;
}


/*
 *===========================================================================
 * Specify the possible discrete changes (jumps) in the individual state
 * variables when ENTERING the stage specified by 'lifestage[]'.
 *
 * Notice that the first index of the variables 'istate[][]' and 'growth[][]'
 * refers to the number of the structured population, the second index refers
 * to the number of the individual state variable. The interpretation of the
 * latter is up to the user.
 *===========================================================================
 */

void DiscreteChanges(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
        double *birthstate[POPULATION_NR], int BirthStateNr, double E[])
{
  if ((lifestage[0] == 2) &&     ((XMIN+XJ-X0) > 0)) 
    {
      REV = REV-META*(QJ-QS)*max(0,(XMIN+XJ - X0));
      SetSurvival(0,(1-META*Rho)*exp(istate[0][I_STATE_DIM]));
    }

  return;
}


/*
 *===========================================================================
 * Specify the fecundity of individuals as a function of the i-state
 * variables and the individual's state at birth for all populations in every
 * life stage.
 *
 * The number of offspring produced has to be specified for every possible
 * state at birth in the variable 'fecundity[][]'. The first index of this
 * variable refers to the number of the structured population, the second
 * index refers to the number of the birth state.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 *===========================================================================
 */

void Fecundity(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
               double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
               double *fecundity[POPULATION_NR])
{
    const double	  MAX_EXP = 50.0;
    double	  pw;
    double          w, wr, h,a1,a2, gamma, Eg, Em,Ei,kappa_qj,kappa_qa;
    double          In1,In2,wmin, PHI,A1,A2,PSI_A;
    
    PSI_A = min(1, PSI+META);                        //Specialization R1 J&A
    w = (1+QJ)*IR;                                  // Effective body mass
    wmin = (1+QJ)*XMIN;                               // Min Weight for R2 (nu wanneer larval stage is beendig)
    h = (CHI1 + CHI2*pow(w, -CHI3)*exp(CHI4*w));    //handling time
    
    if ((lifestage[0] == 0) || (lifestage[0] == 1))
    {
        A1 = (1-PSI)*(AMAX-AMIN)+AMIN;
        A2 = PSI*(AMAX-AMIN)+AMIN;
    }
    else
    {
        A1 = (1-PSI_A)*(AMAX-AMIN)+AMIN;
        A2 = PSI_A*(AMAX-AMIN)+AMIN;
    }
    
    wr = w/W0;
    a1 = A1*pow(wr, ALPHA) * exp(ALPHA*(1 - wr));
    
    if (lifestage[0] == 0)
    {
        a2 = 0;
    }
    else
    {
        wr  = (w-wmin)/W0;
        pw  = log(wr + epsMach) + (1-wr);
        pw *= ALPHA;
        pw  = max(pw, -MAX_EXP);
        pw  = min(pw,  MAX_EXP);
        a2  = A2*exp(pw);
    }
    
    pw = -SIGMA * (a1  * R1 - a2  * R2);
    pw = max(pw, -MAX_EXP);
    pw = min(pw,  MAX_EXP);
    pw = exp(pw);
    
    
    PHI = 1 / (1 + pw);                                    //Time on R1
    gamma = (a1*PHI*R1+a2*(1-PHI)*R2) / (1 + h*(a1*PHI*R1+a2*(1-PHI)*R2));      //Intake rate
    In1 =(a1*PHI*R1) / (1 + h*(a1*PHI*R1+a2*(1-PHI)*R2));                           //Intake rate R1
    In2 = (a2*(1-PHI)*R2) / (1 + h*(a1*PHI*R1+a2*(1-PHI)*R2));                      //Intake rate R2
    Em = M1*pow(IR+REV, M2);                                                                    //maintenance cost
    Ei = (K1)*gamma;                                                                              //Conversion efficiency
    Eg = Ei - Em;                                                                               // Net production
    kappa_qj = 1 / ((1+QJ)*QJ) * (REV / IR);                                                      //fraction allocated to irreversible mass juve
    kappa_qa = 1 / ((1+QA)*QA) * (REV / IR);
    
    
    
    if (lifestage[0] < 3)
    {
        fecundity[0][0] = 0.0;
    }
    else
    {
        double surv;
        double eggenergy;
        double fecfrac = 1-(1+QJ)*QJ/((1+QA)*QA); //
        double s = (1000 * (REV / (QJ * IR)) - 999);   ///stepfunctie betweeb 0 and 1 (0 if rev<qj*x, 1 if rev=>qj*x)
        s = (s < 0) ? 0 : s;
        s = (s > 1) ? 1 : s;
        fecfrac *= 3.0*s*s - 2.0*s*s*s;  ///sigmoid
        if ((XMIN+XJ-X0) <= 0) {  //Size at metamorphose is before size at birth
            surv = min(1,(1-META*Rho)); //Prob to survive Meta
            eggenergy = (1+QJ)*X0; //Cost of an egg
        }
        else {
            surv = 1;
            eggenergy =(1+QJ)*X0;
        }
        
        fecundity[0][0] = surv*K2*fecfrac*Eg/(eggenergy);//Fecundity of adults
    }
    
    return;
}



/*
 *===========================================================================
 * Specify the mortality of individuals as a function of the i-state
 * variables and the individual's state at birth for all populations in every
 * life stage.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 *===========================================================================
 */

void Mortality(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
        double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
        double mortality[POPULATION_NR])
{
 
    mortality[0] = MU0;
 
  return;
}


/*
 *===========================================================================
 * For all the integrals (measures) that occur in interactions of the
 * structured populations with their environments and for all the integrals
 * that should be computed for output purposes (e.g. total juvenile or adult
 * biomass), specify appropriate weighing function dependent on the i-state
 * variables, the environment variables and the current life stage of the
 * individuals. These weighing functions should be specified for all
 * structured populations in the problem. The number of weighing functions
 * is the same for all of them.
 *
 * Notice that the first index of the variables 'istate[][]' and 'impact[][]'
 * refers to the number of the structured population, the second index of the
 * variable 'istate[][]' refers to the number of the individual state variable,
 * while the second index of the variable 'impact[][]' refers to the number of
 * the interaction variable. The interpretation of these second indices is up
 * to the user.
 *===========================================================================
 */

void Impact(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
            double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
            double impact[POPULATION_NR][INTERACT_DIM])
{
    const double	  MAX_EXP = 50.0;
    double	  pw;
    double          w, wr, h,a1,a2, gamma, Eg, Em,Ei,kappa_qj,kappa_qa;
    double          In1,In2,wmin, PHI,A1,A2,PSI_A;
    
    
    PSI_A = min(1, PSI+META);                        //Specialization R1 J&A
    w = (1+QJ)*IR;                                  // Effective body mass
    wmin = (1+QJ)*XMIN;                               // Min Weight for R2 (nu wanneer larval stage is beendig)
    h = (CHI1 + CHI2*pow(w, -CHI3)*exp(CHI4*w));    //handling time
    
    if ((lifestage[0] == 0) || (lifestage[0] == 1))
    {
        A1 = (1-PSI)*(AMAX-AMIN)+AMIN;
        A2 = PSI*(AMAX-AMIN)+AMIN;
    }
    else
    {
        A1 = (1-PSI_A)*(AMAX-AMIN)+AMIN;
        A2 = PSI_A*(AMAX-AMIN)+AMIN;
    }
    
    wr = w/W0;
    a1 = A1*pow(wr, ALPHA) * exp(ALPHA*(1 - wr));
    
    
    if (lifestage[0] == 0)
    {
        a2 = 0;
        
    }
    else
    {
        wr  = (w-wmin)/W0;
        pw  = log(wr + epsMach) + (1-wr);
        pw *= ALPHA;
        pw  = max(pw, -MAX_EXP);
        pw  = min(pw,  MAX_EXP);
        a2  = A2*exp(pw);
    }
    
    pw = -SIGMA * (a1  * R1 - a2  * R2);
    pw = max(pw, -MAX_EXP);
    pw = min(pw,  MAX_EXP);
    pw = exp(pw);
    
    
    PHI = 1 / (1 + pw);                                     //Time on R1
    gamma = (a1*PHI*R1+a2*(1-PHI)*R2) / (1 + h*(a1*PHI*R1+a2*(1-PHI)*R2));      //Intake rate
    In1 =(a1*PHI*R1) / (1 + h*(a1*PHI*R1+a2*(1-PHI)*R2));                           //Intake rate R1
    In2 = (a2*(1-PHI)*R2) / (1 + h*(a1*PHI*R1+a2*(1-PHI)*R2));                      //Intake rate R2
    Em = M1*pow(IR+REV, M2);                                                                    //maintenance cost
    Ei = (K1)*gamma;                                                                              //Conversion efficiency
    Eg = Ei - Em;                                                                               // Net production
    kappa_qj = 1 / ((1+QJ)*QJ) * (REV / IR);                                                      //fraction allocated to irreversible mass juve
    kappa_qa = 1 / ((1+QA)*QA) * (REV / IR);
    
    
    
    
    double	maxage = 1, birthnul = 1;
    if (Birthrate(0) < 1E-9) {
        maxage = 0.0;
        birthnul = 0.0;
    } else {
        maxage = 1 / (Birthrate(0) * Survival(0));
        birthnul = 1 / Birthrate(0);
    }
    
    
    
    switch (lifestage[0])
    {
        case 0:  // Larvae
            impact[0][0] = In1;                                   // Ingestion R1 L
            impact[0][1] = In2;                                   // Ingestion R2 L
            impact[0][2] = IR;                                    // IR    L
            impact[0][3] = REV;                                   // REV   L
            impact[0][4] = 1;                                     // #L
            impact[0][5] = 0;
            impact[0][6] = 0;
            impact[0][7] = 0;
            impact[0][8] = 0;
            impact[0][9] = 0;
            impact[0][10] = 0;
            impact[0][11] = 0;
            impact[0][12] = 0;
            impact[0][13] = 0;
            impact[0][14] = 0;
            impact[0][15] = 0;
            impact[0][16] = 0;
            impact[0][17] = 0;
            impact[0][18] = 0;
            impact[0][19] = 0;
            impact[0][20] = maxage;                               //Age at metamorphosis
            impact[0][21] = maxage;                               //Age at maturation
            impact[0][22] = maxage;                               //Max Age
            impact[0][23] = (kappa_qj*Eg)*maxage;                 //Max size
            
            break;
        case 1:
            impact[0][0] = 0;                                   // Ingestion R1 L
            impact[0][1] = 0;                                   // Ingestion R2 L
            impact[0][2] = 0;                                    // IR    L
            impact[0][3] = 0;                                   // REV   L
            impact[0][4] = 0;                                     // #L
            impact[0][5] = In1;                                   // Ingestion R1 L
            impact[0][6] = In2;                                   // Ingestion R2 L
            impact[0][7] = IR;                                    // IR    L
            impact[0][8] = REV;                                   // REV   L
            impact[0][9] = 1;                                     // #L
            impact[0][10] = 0;
            impact[0][11] = 0;
            impact[0][12] = 0;
            impact[0][13] = 0;
            impact[0][14] = 0;
            impact[0][15] = 0;
            impact[0][16] = 0;
            impact[0][17] = 0;
            impact[0][18] = 0;                                   // Ingestion R1 L
            impact[0][19] = 0;                                   // Ingestion R2 L
            impact[0][20] = 0;                                    // IR    L
            impact[0][21] = 0;                                   // REV   L
            impact[0][22] = 0;
            impact[0][20] = maxage;                               //Age at metamorphosis
            impact[0][21] = maxage;                               //Age at maturation
            impact[0][22] = maxage;                               //Max Age
            impact[0][23] = (kappa_qj*Eg)*maxage;                 //Max size
            break;
            
        case 2:    // Juv
            impact[0][0] = 0;
            impact[0][1] = 0;
            impact[0][2] = 0;
            impact[0][3] = 0;
            impact[0][4] = 0;
            impact[0][5] = 0;
            impact[0][6] = 0;
            impact[0][7] = 0;
            impact[0][8] = 0;
            impact[0][9] = 0;
            impact[0][10] = In1;                                   // Ingestion R1 J
            impact[0][11] = In2;                                   // Ingestion R2 J
            impact[0][12] = IR;                                    // IR     J
            impact[0][13] = REV;                                   // REV    J
            impact[0][14] = 1;                                     // #Juv
            impact[0][15] = 0;
            impact[0][16] = 0;
            impact[0][17] = 0;
            impact[0][18] = 0;
            impact[0][19] = 0;
            impact[0][20] = 0;                                     //Age at metamorphosis
            impact[0][21] = maxage;                               //Age at maturation
            impact[0][22] = maxage;                               //Max Age
            impact[0][23] = (kappa_qj*Eg)*maxage;                 //Max size
            break;
            
        case 3:     // Ad
            impact[0][0] = 0;
            impact[0][1] = 0;
            impact[0][2] = 0;
            impact[0][3] = 0;
            impact[0][4] = 0;
            impact[0][5] = 0;
            impact[0][6] = 0;
            impact[0][7] = 0;
            impact[0][8] = 0;
            impact[0][9] = 0;
            impact[0][10] = 0;
            impact[0][11] = 0;
            impact[0][12] = 0;
            impact[0][13] = 0;
            impact[0][14] = 0;
            impact[0][15] = In1;                                   // Ingestion R1 A
            impact[0][16] = In2;                                   // Ingestion R2 A
            impact[0][17] = IR;                                    // IR    A
            impact[0][18] = REV;                                   // REV   A
            impact[0][19] = 1;
            impact[0][20] = 0;                                      //Age at metamorphosis
            impact[0][21] = 0;                                      //Age at maturation
            impact[0][22] = maxage;                                 //Max Age
            impact[0][23] = (kappa_qj*Eg)*maxage;                   //Max size
            
            break;
            
    }
    
    return;
}


/*
 *===========================================================================
 * Specify the type of each of the environment variables by setting
 * the entries in EnvironmentType[ENVIRON_DIM] to PERCAPITARATE, GENERALODE
 * or POPULATIONINTEGRAL based on the classification below:
 *
 * Set an entry to PERCAPITARATE if the dynamics of E[j] follow an ODE and 0
 * is a possible equilibrium state of E[j]. The ODE is then of the form
 * dE[j]/dt = P(E,I)*E[j], with P(E,I) the per capita growth rate of E[j].
 * Specify the equilibrium condition as condition[j] = P(E,I), do not include
 * the multiplication with E[j] to allow for detecting and continuing the
 * transcritical bifurcation between the trivial and non-trivial equilibrium.
 *
 * Set an entry to GENERALODE if the dynamics of E[j] follow an ODE and 0 is
 * NOT an equilibrium state of E. The ODE then has a form dE[j]/dt = G(E,I).
 * Specify the equilibrium condition as condition[j] = G(E,I).
 *
 * Set an entry to POPULATIONINTEGRAL if E[j] is a (weighted) integral of the
 * population distribution, representing for example the total population
 * biomass. E[j] then can be expressed as E[j] = I[p][i]. Specify the
 * equilibrium condition in this case as condition[j] = I[p][i].
 *
 * Notice that the first index of the variable 'I[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the interaction variable. The interpretation of the latter
 * is up to the user. Also notice that the variable 'condition[j]' should
 * specify the equilibrium condition of environment variable 'E[j]'.
 *===========================================================================
 */

const int EnvironmentType[ENVIRON_DIM] = {GENERALODE,GENERALODE};

void EnvEqui(double E[], double I[POPULATION_NR][INTERACT_DIM],
        double condition[ENVIRON_DIM])
{
    condition[0] = delta*(R1MAX - R1) - I[0][0] - I[0][5] - I[0][10] - I[0][15];
    condition[1] = delta*(R2MAX - R2) - I[0][1] - I[0][6] - I[0][11] - I[0][16];
    return;
}

/*==============================================================================*/

