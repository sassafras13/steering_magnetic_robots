% MAIN.m
% Purpose: Runs through all parts of the gradient_swimmers codebase in
% order. Demonstrates how to use code on data obtained as part of [1].
% Author: Emma Benjaminson
% References:
% [1] Benjaminson, E., Travers, M., Taylor, R. E. "Steering Magnetic Robots
% in Two Axes with One Pair of Maxwell Coils." IROS 2020. 

%% Clear Workspace

clc ; clear all ; close all ; 

%% Directories

localDir = fileparts(mfilename('fullpath')) ;
restoredefaultpath ;
addpath(fullfile(localDir, 'active_functions')) ;
addpath(fullfile(localDir, 'experimental_data')) ; 

%% Load Variables

run(fullfile(localDir,'VARIABLES.m')) ; 

%% 1: Compare Model of Magnetic Field to Experimental Data
% This section plots experimentally collected magnetic field data and
% compares it to the model.

run(fullfile(localDir,'MAIN_01_BFieldDataProcessing.m')) ; 

%% 2: Process Experimental Data
% This section computes average trajectories over multiple runs. The data
% is obtained using the particle tracking tools in ImageJ. 

% run(fullfile(localDir,'MAIN_02_ExpDataProcessing.m')) ; 

%% 3: Generate Magnetic Field Data
% This section generates the magnetic field data for each coil
% configuration used in the control sequence. 

run(fullfile(localDir, 'MAIN_03_GenerateFieldData.m')) ; 

%% 4: Find Drag Coefficient
% This section runs gradient descent to fit the drag coefficient to the
% experimental data. 

% run(fullfile(localDir, 'MAIN_04_DragGradDescent.m')) ; 

%% 5: Compute Swimmer Dynamics
% Plot the swimmers' motion for the chosen control sequence. 

run(fullfile(localDir, 'MAIN_05_DynamicsModel.m')) ; 

%% 6: Plot Results
% Plot the results from the dynamics simulation. 

run(fullfile(localDir, 'MAIN_06_PlotDynamics.m')) ; 
