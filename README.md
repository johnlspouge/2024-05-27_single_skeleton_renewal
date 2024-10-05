# 2024-05-27_single_skeleton_renewal

**To download a single directory from a GitHub project, 
<br>clone the project and then delete the other directories.** 

The project computes values relevant to a branching process that corresponds to an SEIR process with gamma-distributed incubation (exposed) and infectious periods.

Contains 3 directories:
1. **Data/**
2. **Executable/**
3. **Output/**

**Executable/** contains 4 subdirectories:

0. modules/ : reusable subroutines for the other executables below.
modules/ contains a test file test_make.py.
> python test_make.py

tests the installation and functionalities of the Python files in modules/. If all is well, the program exits without complaint.

The remaining subdirectories contain executables with specific functions. 

1. Durations/
2. Generations/
3. Fecundity/

Each of the 3 subdirectories contains 4 files.

0. requirements.txt : contains the version numbers of the Python packages used by the program.
1. run_[executable].py : A program built with the Python argparse package, so '-h' displays program arguments.
<br>'-h' displays the purpose of the program (also summarized below) and its arguments. 
2. run_[executable]_make.py : A Makefile that runs the program.
<br>Each Makefile therefore contains an example of the Python command to run the program, along with the input and output filenames relevant to the program. My programming conventions use corresponding **Data/** and **Output/** directories, but the filenames may be arbitrary. Each Makefile provides an example of program usage, so examining it shows how to run the program, the program arguments, and the input and output filenames for the example.
3. run_[executable].log : A plain text file of system output, sometimes blank. The file received the print statements from the Python program.
<br>In a long running program, the text file can indicate where a failure occurred or how long execution took.

The executables 1, 2, and 3 above each have corresponding subdirectories in **Data/** and **Output/** for their input and output, reflecting my programming conventions.  

The 3 subdirectories of **Executable/** contain executables with the following purposes, also discoverable by running the executable with the '-h' option.

1. Generations/run_ui_generations.py :
<br> Simulates realizations of the generation counts in the single-skeleton renewal for an SEIR model where E and I are gamma-distributed.
2. Durations/run_ui_durations.py :
<br> The single-skeleton renewal for an SEIR model where E and I are gamma-distributed corresponds to a Galton-Watson (GW) process with a Negative_Binomial( k, p ) offspring distribution. The executable simulates realizations of the generation counts for the GW process.
3. Fecundity/run_fecundity_negative_binomial.py :
<br> In a GW process with Negative_Binomial( k, p ) offspring distribution and different ( k, p ), the executable computes the following for an extinct lineage: the mean number of offspring, the extinction probability, and the total derivative of the extinction probability q with respect to the negative binomial parameter p. 

**Data/**

1. **Input Files for Executable/Generations/run_ui_generations.py**

The Makefile **run_ui_generations_make.py** specifies the following input for the executable **run_ui_generations.py**.  

**Data/Generations/generations.csv** specifies the parameters of the waiting times in the compartments E and I. 
<br> The input is a DataFrame with 5 columns corresponding to E mean, E dispersion, I mean, I dispersion, and the basic reproduction number R_0.  

2. **Input Files for Executable/Durations/run_ui_durations.py**

The Makefile **run_ui_durations_make.py** specifies the following input for the executable **run_ui_durations.py**. 

**Data/Durations/generations.csv** is identical to Data/Generations/generations.csv but could be any CSV with a comparable format). 

3. **Input Files for Executable/Fecundity/run_fecundity_negative_binomial.py**

**run_fecundity_negative_binomial.py** does not require any input files.

**Output/**

1. **Output/Generations/generations.csv**

The output is a comma-separated file containing a DataFrame. The first 5 columns are input; the remainder are output.

The headings follow:
A. e->i_mean : mean of the gamma-distributed incubation period
B. e->i_dispersion : dispersion of the gamma-distributed incubation period
C. i->r_mean : mean of the gamma-distributed infectious period
D. i->r_dispersion : dispersion of the gamma-distributed infectious period
E. s->e:i_R_0 : R_0 for an exposed individual
-. The remaining columns count the realizations with the given number of generations between the primary infection and the long-time most recent common ancestor (MRCA).

2. **Output/Durations/durations.csv**

The output is a comma-separated file containing a DataFrame. The first 5 columns are input, identical to the first 5 columns in **Output/Generations/generations.csv** above.

The remainder are output, as follows:
-. The remaining columns list the durations (days) between the primary infection and the long-time most recent common ancestor (MRCA).

4. **Output/Fecundity/**

All output files pertain Negative_Binomial( k, p ) offspring distributions. The headings in the 3 output files are the same, as follows:
A. the parameter k
-. The remaining headings are values of p.

1. **Output/Fecundity/negative_binomial_q.csv** :
The values in the cross-table for ( k, p ) are 0 if the GW process is not supercritical. Otherwise, they are the total derivative of the extinction probability q.

2. **Output/Fecundity/negative_binomial_derivative.csv** :
The values in the cross-table for ( k, p ) are 1 if the GW process is not supercritical. Otherwise, they are the total derivative of the extinction probability q with respect to the parameter p.

3. **Output/Fecundity/negative_binomial_gamma.csv** :
The values in the cross-table for ( k, p ) are 1 if the GW process is not supercritical. Otherwise, they are the total derivative of the extinction probability q with respect to the negative binomial parameter p.

### Code Verification ###

1. The modules in Executable/modules/ contain unit tests that run by entering
> python [filename]

> python test_make.py
runs all the unit tests in Executable/modules/.
