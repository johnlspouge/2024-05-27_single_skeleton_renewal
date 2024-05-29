# 2024-05-27_single_skeleton_renewal

**To download a single directory from a GitHub project, 
<br>clone the project and then delete the other directories.** 

Contains 3 directories:
1. **Data/**
2. **Executable/**
3. **Output/**

**Executable/** contains 3 subdirectories:

0. modules/ : reusable subroutines that the executables below require.
modules/ contains a test file test_make.py.
> python test_make.py

tests the installation and functionalities of the Python files in modules/. If all is well, the program exits without complaint.

The remaining subdirectories contain executables with specific functions, each containing analogous files. 

1. Durations/
2. Generations/

The **Data/** and **Output/** directories contain these 2 subdirectories, too, so they can contain the input and output from the corresponding executables if desired. Each subdirectory for executables generally contains 4 files, as follows. 

1. run_ui_[executable].py : A program built with the Python argparse package, so '-h' displays program arguments.
<br>'-h' displays the purpose of the program (also summarized below) and its arguments. 
2. run_ui_[executable]_make.py : A Makefile that runs the program.
<br>The Makefile therefore contains an example of the Python command to run the file, along with the input and output filenames relevant to the program, which need not be in the corresponding **Data/** and **Output/** directories. Each Makefile provides an example of program usage, so examining it will yield an example of how to run the program, the program arguments, and the input and output filenames for the example.
3. run_ui_[executable].log : A plain text file of system output, sometimes blank. The file received the print statements from the Python program.
<br>In a long running program, the text file often indicates where a failure occurred or how long execution took.
4. requirements.txt : contains the version numbers of the Python packages used by the program.

The subdirectories of **Executable/** contain the following executables with the following purposes, also discoverable by running the Python program with the '-h' option.

1. Generations/run_ui_generations.py :
<br> Simulates realizations of the generations in the single-skeleton renewal for an SEIR model where E and I are gamma-distributed.
2. Durations/run_ui_durations.py :
<br> Simulates realizations of the durations of the single-skeleton renewal for an SEIR model where E and I are gamma-distributed.

**Input Files for Executable/Generations/run_ui_generations.py**

See **Executable/Generations/**. The Makefile **run_ui_generations_make.py** specifies the following input for **run_ui_generations.py**. The input is in **Data/Generations/**.  

1. **generations.csv** specifies the parameters of the waiting times in the compartments E and I. 
<br> The input is a DataFrame with 5 columns corresponding to E mean, E dispersion, I mean, I dispersion, and the basic reproduction number R_0.  

**Input Files for Executable/Durations/run_ui_durations.py**

See **Executable/Durations/**. The Makefile **run_ui_durations_make.py** specifies the following input for **run_ui_durations.py**. The input is in **Data/Durations/**. 

1. **generations.csv** is identical to Data/Generations/generations.csv. 

### Code Verification ###

1. The modules in Executable/modules/ contain unit tests that run by entering
> python [filename]

> python test_make.py
runs all the unit tests in Executable/modules/.
