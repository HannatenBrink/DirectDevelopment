The command 'make' will compile the program "RunIBM".
The program needs two files, a .csv file with all the parameters of the model, and an .isf file with the initial state of the population. 
The command './RunIBM NAME' will run the IBM, where NAME is the name of the cvs & isf file. 
For example, './RunIBM Default' will run the ecological dynamics for 10000 days, for the parameters specified in Default.csv. Since the mutation probabilities are zero in this file, there is no evolution.

The first row of the ISF file has 4 values, the starttime, density of R1, density of R2, and the size of the system (in liters), separated by a tab.
The rest of the file contains the data for the initial population. Each line represents one (or more) individual(s) with its individual states (number of individuals (1 or higher), its age, irreversible and reversible mass, if it has metamorphosed, and its four traits: extend of metamorphosis, psi_L, irreversible mass at metamorphosis, Irreversible mass at birth), separated by a tab. 
In the Default.isf file, for example, the initial population consist of 100 newborn individuals. 

The IBM creates a minimum of 2 output files, one file with the population & trait dynamics over time (NAME_time.txt), and one file with the state of the population (NAME.ESF) at the end of the run. If the third parameter of the .csv file ("Complete state output interval, 0 for none") has a value larger than zero, the program will create every X time steps a file (NAME_Popfile_CURRENTTIME.esf) with the current state of the population, where X can be specified in the .csv file. 

The file with the dynamics over time contains the following information: "time", "Density R1", "Density R2", "Number of larvae", "Biomass of Larvae", "Number of J", "Biomass of J", "Number of A", "Biomass of A" ,"mean trait Psi_L", "mean trait Meta", "mean trait X_J", "mean trait X_b", "Xmax1".

The file(s) with the state of the population (.ESF) has the same structure as the .isf file, with the first line of the file containing the environmental values, and the rest of the file containing the information of the population. These .esf files can be used for a new run.
The program will normally use the trait values specified in the .csv file to initialise the population. However, this can be changed in the 'Pars.h' file, by changing the value of TRAITS to 1, and recompiling the code. In this case, the program will use the trait values specified in the .isf file to initialise the population. This is, for example, useful to continue an evolutionary run.  



