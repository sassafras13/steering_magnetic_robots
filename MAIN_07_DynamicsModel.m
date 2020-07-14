% MAIN_DynamicsModel.m
% Purpose: Model the motion of a single link microswimmer in a magnetic 
% field produced by a pair of Helmholtz or Maxwell coils. For comparison 
% with experimental results to show model is accurate.
% Author: Emma Benjaminson
% References: 

%% 

clc ; clear all ; close all ; 

%% Directories

localDir = fileparts(mfilename('fullpath')) ;
restoredefaultpath ;
addpath(fullfile(localDir, 'active_functions')) ;

%% Load Variables

% load macroswimmer data
run VAR_SingleLinkModelCoil.m ; 
load('MAT_gBarray.mat') ; 
load('MAT_Barray.mat') ; 

% load experimental data
% directory = sprintf('%s/TP12','/home/emma/repos/microswimmers/Gradient Steering') ; 
directory = '/home/emma/Documents/Research/TP12/Composites/Case 1/Mod' ; 
filename = 'c1' ; 
[indices, origins] = extractIndOr(directory,filename) ; 
filename = 'Case1-CompositeData' ; 
loadExpData(indices, origins, directory, filename, 0) ; 

% % load experimental data indices of origins
% directory = '/home/emma/Documents/Research/TP12/Motion Primitives/For Analysis' ; 
% % directory = sprintf('%s/TP12',localDir) ; 
% % directory = '/home/biorobotics/microswimmers/Gradient Steering/TP12' ; 
% filename = 'mpA' ; 
% [indices, origins] = extractIndOr(directory,filename) ; 
% 
% % load all experimental data for a particular motion primitive
% % filename_exp_data = 'tp12_mpA_iter' ; 
% filename_exp_data = 'mpA-CompositeData' ; 
% loadExpData(indices, origins, directory, filename_exp_data, 0)

%%%%%%%%%% USER DEFINED VALUE %%%%%%%%%%%%%%%%%%%%%%%%%%%
% select experimental data
index = 1 ; % index of the trial to plot

% waypoint times - when the swimmer changes direction
% waypointT = [0, 170] ; % [s] for mpA
% waypointT = [0,225,470,670] ; % [s] 
% waypointT = [0,130,255,365] ; % [s]
waypointT = [0,65,110,200] ; % [s] Case 1
% waypointT = [0,140,260,400,500,600,750,900,1050,1200] ; 

% number of primitives in trajectory
niter = 3 ; 

% initial conditions (USER DEFINES THETA!)
% y0 = [x ; y ; theta ; xdot ; ydot ; thetadot] ; 
y0 = [meanposx(1,1) ; meanposy(1,1) ; pi/2 ; 0 ; 0 ; 0 ]  ;
% y0 = [posx{index}(1,1) ; posy{index}(1,1) ; 0.7854 ; 0 ; 0 ; 0] ; 

% control inputs
I1 = [1.9  ;
      0 ; 
      1.9 ]; 
%       0 ; 
%       1.9 ; 
%       0 ; 
%       1.9 
%       0 ; 
%       1.9 ];  % [A] current
I2 = [0 ; 
      1.9 ; 
      0 ]; 
%       1.9 ; 
%       0 ; 
%       1.9 ; 
%       0 ; 
%       1.9 ; 
%       0 ] ; % [A] current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Dynamics Simulation
  
yarray = {} ; 
timearray = {} ; 

for i = 1:niter
    
    % time 
    T0 = waypointT(i) ; 
    Tf = waypointT(i+1) ; 
    Tspan = [T0 Tf] ; 

    options = odeset('RelTol',1E-3,'AbsTol',1E-6,'Stats','on','OutputFcn',...
        @(t,y,flag) myOutputFcnCoil(t,y0,flag,mag,mass,w,L,ct,cn,Barray{i},gBarray{i})) ;

    tic 
    [t,y] = ode15s(@(Tspan, y0) singleDynamicsCoil(Tspan,y0,mag,...
        mass,w,L,ct,cn,Barray{i},gBarray{i}), Tspan, y0, options) ;  
    toc
    
    yarray{i} = y ; 
    timearray{i} = t ; 
    
    y0 = [y(end,1) ; y(end,2) ; y(end,3) ; 0 ; 0 ; 0] ;
    
end

% save all data to a .mat file
save('MAT_Dynamics_case1') ; 
