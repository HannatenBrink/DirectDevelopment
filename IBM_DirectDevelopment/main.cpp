#include "main.h"



int main(int argc, char* argv[]) {
  //Check if correct input is given
  if(argc < 2 ){
    std::cerr << "Usage: " << argv[0] << " Name of csv/isf file" << std::endl;
    exit(1);
  }
  //Use second argument as filename
  string Filename = argv[1];
  cout << "Start of run: " << Filename << std::endl;
  //////////Read Cvf file/////////////////
  

  ss<<Filename<<".cvf";
  std::ifstream input(ss.str());

  if(!input){
    std::cerr<<"The parameterfile " <<ss.str()<< " does not exist.\n";
    exit(1);
  }

  std::string line;
  string search = "[ 0]"; // test variable to search in file
  unsigned int curLine = 0; //Var to keep track of line numbers
  unsigned int curLine2 = 0; //Var2 to keep track of line numbers
  int Start_read = 0;


  while(std::getline(input, line))
  {
    curLine++;
    vector<string> lines;
    if (Start_read == 0){
      split(line, '\t', lines);
      string b = lines.back();
      double AStand = strtof((b).c_str(),0);
      Setting.push_back(AStand); //Read in time settings data in vector 'Setting'
      if (line.find(search) != string::npos) {
        Start_read = 1; //Start with reading parameter values
      }}
      if ((Start_read == 1)){
        curLine2++;
        vector<string> row_values;
        split(line, '\t', row_values);
        string a =  row_values.back();
        double aNumero=strtof((a).c_str(),0);
        Parameter.push_back(aNumero); //Add parameter value to the vector 'Parameter'
        //cout << "Par is: " << aNumero << std::endl; //Print all parameter values
      }
    }
    ss.str("");

    //////Finished reading cvf file////////////
    cout << "psi_L is: " << psi_L << std::endl;
    cout << "Meta is: " << Meta << std::endl;
    cout << "x_j is: " << x_j  << std::endl;
    cout << "xb is: " << xb  << std::endl;
    //////////INITIALIZE THE POPULATION/////////////////////////////////////////
    //Init the population with an isf file
    ss<<Filename<<".isf";
    std::ifstream initPop(ss.str());
    //Check if file exist

   if(!initPop){
      std::cerr<<"The init popfile " <<ss.str()<< " does not exist.\n";
      exit(1);
    }
    ss.str("");


    //Create a vector of individuals, called Population
    std::vector<Individual> Population;
    //Init all the values that will be in isf file///
    double R1_init;
    double R2_init;
    double Liters_init;
    double Starttime_init;
    double density_init;
    double age_init;
    double x_init;
    double y_init;
    double has_meta_init;

    double empty; //dummy variable to ignore traits from isf file
    double Vol = 1; //variable to increase the density when volume in isf file is < cvf file. When the volume in isf file > cvf file, there are too many individuals. They will die automatically...

    int i = 0;
    while (std::getline(initPop, line)){
      std::stringstream linestream(line);
      if (i == 0) {
        linestream >> Starttime_init >> R1_init >> R2_init >> Liters_init;
        if (Liters_init < Liters) {
          Vol = (Liters/Liters_init);
        }
        else
        Vol = 1;
      }
      if (i > 1) {
        if (TRAITS){ //Traits of each individual in the ISF file, parameters of cvf file replaced with these values
        linestream >> density_init >> age_init >> x_init >> y_init >> has_meta_init >> Meta >> psi_L >> x_j >>xb;
      } else {//Start with the values of the cvf file
        linestream >> density_init >> age_init >> x_init >> y_init >> has_meta_init >> empty >> empty >> empty >>empty;
      }
        for (int p = 0; p < round(density_init*Vol); p++) {
          Population.push_back(Individual(age_init, x_init, y_init, has_meta_init, Meta, psi_L, x_j, xb));
        }
      }
      i += 1;
    }
    //close file
    initPop.close();
    //Check if the population is present//

    if(Population.size()<=0){
      std::cerr<<"The consumer population is absent in the isf file\n";
      exit(1);
    }
    cout << "Initial size of the population: " << Population.size() << std::endl;

    std::sort(Population.begin(),Population.end());


    /////////////////////////////////////////////////////////////////////////////////


    //Create names for writing output
    ofstream Timefile;
    ss<<Filename<<"_time"<<".txt";
    Timefile.open(ss.str());
    ss.str("");
    ofstream EndPopfile;
    ss<<Filename<<".esf";
    EndPopfile.open(ss.str());
    ss.str("");
    ofstream BinFile("foobar.bin", std::ios::binary);
    ofstream writer;



    //Set seed for random and create random distributions
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    mt19937 mt_rand(seed);

    //mt19937 mt_rand(1); This one instead of the two above in case you want to repeat numbers
    std::uniform_real_distribution<double> unif(0, 1);
    //Normal distribution
    std::normal_distribution<double> Norm(0, Mut_var);
    //heat up the random number generator (no idea why/if this is necessary)
    for (int n = 0; n < 10000; n++){
      unif(mt_rand);
      Norm(mt_rand);
    }

    ///Getting random numbers is relatively slow. Therefore, produce vectors with random values and iterate over these vectors to get random values
        int RandSize = 100000;
        std::vector<double> MutVal;
        //Make a vector with Mutation values
        for (int i=0; i < RandSize; i++){
          MutVal.push_back(Norm(mt_rand));
        }

        //Make a vector with values between 0 and 1
        std::vector<double> UniVal;
        for (int i=0; i < RandSize; i++){
          UniVal.push_back(unif(mt_rand));
        }

        //Integer to iterate over random vectors
        int NormInt = 0;
        int UniInt = 0;

    //Init the environment
    Resource R1(R1_init, Xmax1, Liters, "R1");
    Resource R2(R2_init, Xmax2, Liters, "R2");
    R1.printInfo();
    R2.printInfo();

    //Create an Iterator
    std::vector<Individual>::iterator it;  //Iterator

    //Create an empty vector of individuals, to store newborns
    std::vector<Individual> New_Cohort;

    //Define the newtime
    double newtime = Total_time / delta_t;

    //start the loop
    int time = Starttime_init;

    double Pop_Intake1 = 0;
    double Pop_Intake2 = 0;
    double A = 0;
    double Am = 0;
    double L = 0;
    double Lm = 0;
    double J = 0;
    double Jm = 0;
    double Size = 0;
    double Size2 = 0;
    double meanPsi = 0;
    double meanMeta = 0;
    double meanXj = 0;
    double meanXb = 0;
    //double meanIrr = 0;
    //double meanRev = 0;

    double psi_mut = 0;
    double meta_mut = 0;
    double xj_mut = 0;
    double xb_mut = 0;

    double newPsi;
    double newMeta;
    double newXj;
    double newXb;
    double T_ESF = 0; //For output
    double T_Time = 0; //For output
    double Rmax1New = Xmax1; //For output
    while (time < (newtime)) {
      //Set everything back to zero

      Pop_Intake1 = 0;
      Pop_Intake2 = 0;
      A = 0;
      Am = 0;
      L = 0;
      Lm = 0;
      J = 0;
      Jm = 0;
      Size = 0;
      Size2 = 0;
      meanPsi = 0;
      meanMeta = 0;
      meanXj = 0;
      meanXb = 0;

      for (it = Population.begin(); it != Population.end(); ++it){
        it->R_Intake(R1.Density, R2.Density); //Eat resources
        Pop_Intake1 += it->Intake1;
        Pop_Intake2 += it->Intake2;
        it->Grow(); //Let individuals grow
        Size = it->x;
        if ((Size >= it->xj_trait) && (it->Has_meta == 0)) { //Let individuals metamorphose
          it->Metamorphose(UniVal[UniInt]);
          UniInt +=1;
          if (UniInt >= RandSize){
            UniInt = 0;
            std::random_shuffle(UniVal.begin(), UniVal.end());
          }
        }
        if (it->Is_dead == 0) {  //Kill individuals that are not dead yet
          it->Die(UniVal[UniInt]);
          UniInt +=1;
          if (UniInt >= RandSize){
            UniInt = 0;
            std::random_shuffle(UniVal.begin(), UniVal.end());
          }
        }
        if (it->Is_dead == 0) {
        if (Size < x_j){
          L += 1;
          Lm += it->w;
        }
        else if (Size < x_A){
          J += 1;
          Jm += it->w;
        }
        else {
          A += 1;
          Am += it->w;
        }
        meanPsi += it->psi_trait;
        meanXj += it->xj_trait;
        meanMeta += it->meta_trait;
        meanXb += it->xb_trait;
        //meanIrr += it->x;
        //meanRev += it->y;
        //it->printInfo();
      }}

      //Remove dead consumers
      Population.erase(std::remove_if(Population.begin(),Population.end(), IsMarkedToDelete), Population.end());


      meanPsi = meanPsi/Population.size();
      meanMeta = meanMeta/Population.size();
      meanXj = meanXj/Population.size();
      meanXb = meanXb/Population.size();
      //meanIrr = meanIrr/Population.size();
      //meanRev = meanRev/Population.size();



      //Resource growth
      R1.Growth(Pop_Intake1, Sr * time * delta_t); //Let Resources grow
      R2.Growth(Pop_Intake2, 0);  //Let resources grow
      Rmax1New = Xmax1 - Sr * time * delta_t; //New Rmax1value, needed for output

      //Reproduction and mutation
      for (it = Population.begin(); it != Population.end(); ++it){
        Size2 = it->x;
        if (Size2 > x_A){  //Check if individual has matured
          it->Repro(); //Let it reproduce
          double Off_num = it->Offspring;  //Get number of offspring produced
          if (Off_num > 0 ) {
            for (int i = 0; i < Off_num; i++){  //For each newborn determine if it will get a mutation & determine new trait
              psi_mut = UniVal[UniInt];
              UniInt+=1;
              if (UniInt >= RandSize){
                UniInt = 0;
                std::random_shuffle(UniVal.begin(), UniVal.end());
              }
              meta_mut = UniVal[UniInt];
              UniInt+=1;
              if (UniInt >= RandSize){
                UniInt = 0;
                std::random_shuffle(UniVal.begin(), UniVal.end());
              }
               xj_mut = UniVal[UniInt];
              UniInt+=1;
              if (UniInt >= RandSize){
                UniInt = 0;
                std::random_shuffle(UniVal.begin(), UniVal.end());
              }
               xb_mut = UniVal[UniInt];
              UniInt+=1;
              if (UniInt >= RandSize){
                UniInt = 0;
                std::random_shuffle(UniVal.begin(), UniVal.end());
              }

              if (mut_rate_psi > psi_mut){
                newPsi = min(1.0, max(0.0, it->psi_trait + MutVal[NormInt]));
                NormInt +=1;
                if (NormInt >= RandSize) {
                  NormInt = 0;
                  std::random_shuffle (MutVal.begin(), MutVal.end() );
                }

              }
              else {
                newPsi = it->psi_trait;
              }
              if (mut_rate_meta > meta_mut) {
                newMeta = min(1.0, max(0.0, it->meta_trait +  MutVal[NormInt]));
                NormInt +=1;
                if (NormInt >= RandSize) {
                  NormInt = 0;
                  std::random_shuffle (MutVal.begin(), MutVal.end() );
                }

              }
              else {
                newMeta = it->meta_trait;
              }
              if (mut_rate_xj > xj_mut){
                newXj = max(xb, it->xj_trait + MutVal[NormInt]);
                NormInt +=1;
                if (NormInt >= RandSize) {
                  NormInt = 0;
                  std::random_shuffle (MutVal.begin(), MutVal.end() );
                }
              }
              else {
                newXj = it->xj_trait;
              }
              if (mut_rate_xb > xb_mut){
                newXb = max(0.00001, it->xb_trait + MutVal[NormInt]);
                NormInt +=1;
                if (NormInt >= RandSize) {
                  NormInt = 0;
                  std::random_shuffle (MutVal.begin(), MutVal.end() );
                }
              }
              else {
                newXb = it->xb_trait;
              }
              //New individuals get the size of the xb_trait of the parents! But will produce offspring with size newXB
              New_Cohort.push_back(Individual(0, it->xb_trait, it->xb_trait * 0.742, 0, newMeta, newPsi, newXj, newXb)); //Add to newborn cohort
            }
          }
        }}

        //Add newborns to the population
        Population.insert(std::end(Population), std::begin(New_Cohort), std::end(New_Cohort));


        //Clear the newborn vector
        New_Cohort.clear();
        //Write time data

        if (round(fmod(T_Time, (Output_time / delta_t))) <= 0){
          T_Time = 0;
          /*for (it = Population.end() - 1; it != Population.begin() - 1; it--){
            Timefile << 1 << '\t' << it->Age << '\t' << it->x << '\t' << it->y << '\t' << it->Has_meta << '\t' << it->meta_trait << "\t" << it->psi_trait << "\t" << it->xj_trait << "\t" << it->xb_trait<< "\t" << it->Intake1 << "\t" << it->Intake2 << "\t" << it->Attack1 << "\t" << it->Attack2 <<'\n';
          }*/ //Outputdata in case I want to check growth of one individual
        Timefile << time * delta_t << "\t" << R1.Density << "\t" << R2.Density << "\t" << L << "\t" << Lm << "\t" << J << "\t" << Jm << "\t" << A << "\t" << Am << "\t" << meanPsi << "\t" <<meanMeta<< "\t" <<meanXj<< "\t" <<meanXb<< "\t" << Rmax1New << "\n";
         }

        //Exit when population is extinct
      if(Population.size()<=0){
          std::cerr<<"The consumer population is extinct at time " << time * delta_t << "\n";
          exit(1);
        }

        //Write POPstate output//
        if ((Pop_Output > 0) && (round(fmod(T_ESF, (Pop_Output / delta_t))) == 0)){
          T_ESF = 0;
          int timename = int(time * delta_t); //This one if you want to create new timefiles
          //string timename = "std"; //this one just as a sort of backup
          ss<<Filename<<"_Popfile_"<<timename<<".esf";
          writer.open(ss.str());
          writer << time * delta_t << "\t" << R1.Density << "\t"<< R2.Density << "\t" << Liters <<"\n\n";
          for (it = Population.end() - 1; it != Population.begin() - 1; it--){
            writer << 1 << '\t' << it->Age << '\t' << it->x << '\t' << it->y << '\t' << it->Has_meta << '\t' << it->meta_trait << "\t" << it->psi_trait << "\t" << it->xj_trait << "\t" << it->xb_trait<< '\n';
          }
          writer.close();
          //Clear the string stream
          ss.str("");

        }

        if(Rmax1New <= 1*pow(10,-16)){
           std::cerr<<"Rmax1 equals zero at time " << time * delta_t << "\n";
           EndPopfile << time * delta_t << "\t" << R1.Density << "\t"<< R2.Density << "\t" << Liters <<"\n\n";
           for (it = Population.end() - 1; it != Population.begin() - 1; it--){
             EndPopfile << 1 << '\t' << it->Age << '\t' << it->x << '\t' << it->y << '\t' << it->Has_meta << "\t" << it->meta_trait << "\t" << it->psi_trait << "\t" << it->xj_trait << "\t" << it->xb_trait<< '\n';
           }

           Timefile.close();
           EndPopfile.close();
           exit(1);
         }
        //Next step in loop
        time += 1;
        T_Time +=1;
        T_ESF +=1;
      }
      cout << "Run is finished at time: " << time <<std::endl;
      ///Create ESF file
      EndPopfile << time * delta_t << "\t" << R1.Density << "\t"<< R2.Density << "\t" << Liters <<"\n\n";
      for (it = Population.end() - 1; it != Population.begin() - 1; it--){
        EndPopfile << 1 << '\t' << it->Age << '\t' << it->x << '\t' << it->y << '\t' << it->Has_meta << "\t" << it->meta_trait << "\t" << it->psi_trait << "\t" << it->xj_trait << "\t" << it->xb_trait<< '\n';
      }

      Timefile.close();
      EndPopfile.close();

      return 0;
    }
