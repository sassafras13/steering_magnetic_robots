# steering_magnetic_robots
Code supporting IROS 2020 paper "Steering Magnetic Robots in Two Axes with One Pair of Maxwell Coils" 

## MAIN Scripts

These files are the main functions for performing different calculations.    

**(1) MAIN_CoilModel :-** Characterizes the magnetic field and magnetic field gradient produced by a pair of Helmholtz or Maxwell coils. For comparison with other models to prove it is a robust model.    

**(2) MAIN_SingleLinkModelCoil :-** Model the motion of a single link microswimmer in a magnetic field produced by a pair of Helmholtz or Maxwell coils. For comparison with experimental results to show model is accurate.   

**(3) MAIN_DragGradDescent :-** Uses gradient descent to find value of drag coefficient that best fits experimental data.   

**(4) MAIN_PlotDynamics :-** Script that contains code to generate multiple plots as desired for paper.    

**(5) MAIN_SwitchingTimeOpt_SingleLinkModelCoil_WithCost :-** Implementation of switching time optimization.   

**(6) MAIN_SwitchingTimeOpt_SingleLinkModelCoil_WithCost_Case3 :-** Implementation of switching time optimization for test case 3.    

**(7) MAIN_SwitchingTimeOpt_SingleLinkModelCoil_WithCost_Manual_dJdt :-** Implementation of switching time optimization with manual check of gradient of cost function.

## Folders

**(1) active_functions :-** contains the functions currently being used in the various MAIN scripts.    

**(2) Experiment Data Processing :-** contains scripts for processing experimental data obtained in experiments.    

## VAR Files

Scripts beginning with VAR contain parameters used by multiple scripts. These parameters describe the Maxwell coil parameters and other relevant values of commonly used variables in the models.

## Workflow for Analyzing Motion Primitives

This section describes how we can take raw video footage of motion primitives and obtain a dataset ready for analysis or use in gradient descent algorithms

**(1)** Process the videos using ImageJ and TrackMate (there is a different .md file explaining this process in more detail). The final output of this step should be a folder containing files ending in "n_spotstats.csv" where n is the iteration number. 

**(2)** Create a file with a name ending in "-indor.csv", containing the measured locations of the origins of all the video files. It should be located in the same folder that contains all of the spotstats files. 

**(3)** Open MAIN_ExpDataProcessing.m and run the code through the functions extractIndOr() and processExpData(). 

**(4)** Open MAIN_DragGradDescent.m and edit the directory name to import the experimental data correctly. Then run the script to iterate on the drag coefficient until the model is very close to the simulation.
