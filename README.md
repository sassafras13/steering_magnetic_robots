# steering_magnetic_robots
Code supporting IROS 2020 paper "Steering Magnetic Robots in Two Axes with One Pair of Maxwell Coils" 

This folder contains scripts for modeling gradient steering using Maxwell coils on microswimmers. 

## MAIN Scripts

These files are the main functions for performing different calculations.    

**(1) MAIN_CoilModel :-** Characterizes the magnetic field and magnetic field gradient produced by a pair of Helmholtz or Maxwell coils. For comparison with other models to prove it is a robust model.    

**(2) MAIN_SingleLinkModelCoil :-** Model the motion of a single link microswimmer in a magnetic field produced by a pair of Helmholtz or Maxwell coils. For comparison with experimental results to show model is accurate.   

**(3) MAIN_SingleLinkModelCoil_ExpData :-** Model the motion of a single link microswimmer in a magnetic field produced by a pair of Helmholtz or Maxwell coils. Also process experimental results for comparison with experimental results to show model is accurate.   

## Folders

**(1) active_functions :-** contains the functions currently being used in the various MAIN scripts.    

**(2) Archived Functions :-** contains old functions no longer used in any of the MAIN scripts above.    

**(3) Images :-** images produced using MAIN scripts.    

**(4) DeTroye Source Code :-** code written by David J. DeTroye and Ronald J. Chase, "The Calculation and Measurement of Helmholtz Coil Fields," Army Research Laboratory, 1994, <https://pdfs.semanticscholar.org/e294/9553610fa6051a3af3420271eda0135427a6.pdf>.    

**(5) Emmas_linearization :-** code written by Steve Crews to demonstrate how to perform partial derivatives with forward and central difference methods.  

## Workflow for Analyzing Motion Primitives

This section describes how we can take raw video footage of motion primitives and obtain a dataset ready for analysis or use in gradient descent algorithms

**(1)** Process the videos using ImageJ and TrackMate (there is a different .md file explaining this process in more detail). The final output of this step should be a folder containing files ending in "n_spotstats.csv" where n is the iteration number. 

**(2)** Create a file with a name ending in "-indor.csv", containing the measured locations of the origins of all the video files. It should be located in the same folder that contains all of the spotstats files. 

**(3)** Open MAIN_ExpDataProcessing.m and run the code through the functions extractIndOr() and processExpData(). 

**(4)** Open MAIN_DragGradDescent.m and edit the directory name to import the experimental data correctly. Then run the script to iterate on the drag coefficient until the model is very close to the simulation.

'/home/emma/repos/microswimmers/Gradient Steering/MAIN_CoilModel.m' 
