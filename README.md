# steering_magnetic_robots
Code supporting IROS 2020 paper "Steering Magnetic Robots in Two Axes with One Pair of Maxwell Coils" 

## MAIN Scripts

These files are the main functions for performing different calculations. Please start with MAIN.m which establishes the codebase workflow. Further details for individual scripts are below.    

**(1) MAIN_01_BFieldDataProcessing :-** Computes the mean and standard deviation for a set of data characterizing the magnetic field. Saves the values to a csv file. Plots the data and compares to theoretical values.      

**(2) MAIN_02_ExpDataProcessing :-** Process the experimental data and save the average values and standard deviations to a csv file for future use.     

**(3) MAIN_03_GenerateFieldData :-** Generate an array of values for the magnetic field and magnetic field gradient.     

**(4) MAIN_04_DragGradDescent :-** Uses gradient descent to find value of drag coefficient that best fits experimental data.   

**(5) MAIN_05_DynamicsModel :-** Uses ode15s to simulate the dynamics of the swimmer given the drag coefficient and magnetic field model computed earlier.    

**(6) MAIN_06_PlotDynamics :-** Script that contains code to generate multiple plots as desired for paper.    

**MAIN_SwitchingTimeOpt_SingleLinkModelCoil_WithCost :-** Implementation of switching time optimization.   

**MAIN_SwitchingTimeOpt_SingleLinkModelCoil_WithCost_Case3 :-** Implementation of switching time optimization for test case 3.    

**MAIN_SwitchingTimeOpt_SingleLinkModelCoil_WithCost_Manual_dJdt :-** Implementation of switching time optimization with manual check of gradient of cost function.

## Additional Scripts

**ReynoldsNumber :-** Used to compute the Reynolds number for this work and comparative examples in literature.    

**FitCoilModeltoData :-** Used to apply an optional multivariate polynomial to experimental magnetic field data.   

## Folders

**(1) active_functions :-** contains the functions currently being used in the various MAIN scripts.       

**(2) script_archive :-** contains old scripts no longer in use.    

**(3) experimental_data :-** contains all experimental data used in this work.    

## VAR Files

**VARIABLES :-** These parameters describe the Maxwell coil parameters and other relevant values of commonly used variables in the models.    

## For Questions/Comments/Concerns

Please contact Emma Benjaminson at ebenjami@andrew.cmu.edu     

Please cite our paper in IROS if you draw upon this codebase. Thank you.   
