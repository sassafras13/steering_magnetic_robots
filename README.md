# steering_magnetic_robots
Code supporting IROS 2020 paper "Steering Magnetic Robots in Two Axes with One Pair of Maxwell Coils" 

## MAIN Scripts

These files are the main functions for performing different calculations.    

**(1) MAIN_BFieldDataProcessing :-** Computes the mean and standard deviation for a set of data characterizing the magnetic field. Saves the values to a csv file. Plots the data and compares to theoretical values. 

**(2) MAIN_ExpDataProcessing :-** Process the experimental data and save the average values and standard deviations to a csv file for future use.

**(3) MAIN_CoilModel :-** Characterizes the magnetic field and magnetic field gradient produced by a pair of Helmholtz or Maxwell coils. For comparison with other models to prove it is a robust model. Also compares model to some experimental data.    

**(4) MAIN_FitCoilModeltoData :-** Fit the model of the Maxwell coils and their resultant magnetic field to the experimental data.


**(5) MAIN_GenerateFieldData :-** Generate an array of values for the magnetic field and magnetic field gradient after the optimal fit parameters have been found. 

**(6) MAIN_DragGradDescent :-** Uses gradient descent to find value of drag coefficient that best fits experimental data.   

**(7) MAIN_DynamicsModel :-** Uses ode15s to simulate the dynamics of the swimmer given the drag coefficient and magnetic field model computed earlier. 

**(8) MAIN_PlotDynamics :-** Script that contains code to generate multiple plots as desired for paper.    

**MAIN_SwitchingTimeOpt_SingleLinkModelCoil_WithCost :-** Implementation of switching time optimization.   

**MAIN_SwitchingTimeOpt_SingleLinkModelCoil_WithCost_Case3 :-** Implementation of switching time optimization for test case 3.    

**MAIN_SwitchingTimeOpt_SingleLinkModelCoil_WithCost_Manual_dJdt :-** Implementation of switching time optimization with manual check of gradient of cost function.

**MAIN_ReynoldsNumber :-** Used to compute the Reynolds number for this work and comparative examples in literature. 

## Folders

**(1) active_functions :-** contains the functions currently being used in the various MAIN scripts.       

**(2) script_archive :-** contains old scripts no longer in use.

## VAR Files

Scripts beginning with VAR contain parameters used by multiple scripts. These parameters describe the Maxwell coil parameters and other relevant values of commonly used variables in the models.

## Important Parameters

**(1) Fit Parameters for Magnetic Field Model :-** c1 = 0.8014 ; c2 = 0.2144 

**(2) Drag Coefficient :-** 0.0250

## Workflow for Analyzing Motion Primitives

This section describes how we can take raw video footage of motion primitives and obtain a dataset ready for analysis or use in gradient descent algorithms

**(1)** Process the videos using ImageJ and TrackMate (there is a different .md file explaining this process in more detail). The final output of this step should be a folder containing files ending in "n_spotstats.csv" where n is the iteration number. 

**(2)** Create a file with a name ending in "-indor.csv", containing the measured locations of the origins of all the video files. It should be located in the same folder that contains all of the spotstats files. 

**(3)** Open MAIN_ExpDataProcessing.m and run the code through the functions extractIndOr() and processExpData(). 

**(4)** Open MAIN_DragGradDescent.m and edit the directory name to import the experimental data correctly. Then run the script to iterate on the drag coefficient until the model is very close to the simulation.
