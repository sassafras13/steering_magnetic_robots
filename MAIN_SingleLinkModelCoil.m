% MAIN_SingleLinkModelCoil.m
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

% load experimental data
directory = sprintf('%s/TP12','/home/emma/repos/microswimmers/Gradient Steering') ; 
filename = 'c1' ; 
[indices, origins] = extractIndOr(directory,filename) ; 
filename = 'Case1-CompositeData' ; 
loadExpData(indices, origins, directory, filename, 0) ; 

%%%%%%%%%% USER DEFINED VALUE %%%%%%%%%%%%%%%%%%%%%%%%%%%
% select experimental data
index = 1 ; % index of the trial to plot

% waypoint times - when the swimmer changes direction
waypointT = [0, 200] ; % [s]
% waypointT = [0,225,470,670] ; % [s] 
% waypointT = [0,130,255,365] ; % [s]
% waypointT = [0,65,110,200] ; % [s]
% waypointT = [0,140,260,400,500,600,750,900,1050,1200] ; 

% number of primitives in trajectory
niter = 1 ; 

% initial conditions (USER DEFINES THETA!)
% y0 = [x ; y ; theta ; xdot ; ydot ; thetadot] ; 
y0 = [meanposx(21,1) ; meanposy(21,1) ; pi/2 ; 0 ; 0 ; 0 ]  ;
% y0 = [posx{index}(1,1) ; posy{index}(1,1) ; 0.7854 ; 0 ; 0 ; 0] ; 

% control inputs
I1 = [1.9 ] ;
%       0 ; 
%       1.9 ]; 
%       0 ; 
%       1.9 ; 
%       0 ; 
%       1.9 
%       0 ; 
%       1.9 ];  % [A] current
I2 = [0 ]; 
%       1.9 ; 
%       0 ]; 
%       1.9 ; 
%       0 ; 
%       1.9 ; 
%       0 ; 
%       1.9 ; 
%       0 ] ; % [A] current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate Magnetic and Gradient Fields

Barray = {} ; 
gBarray = {} ; 

for i = 1:niter
    
    % coils
    coils = [I1(i) ; I2(i) ; a ; d ; nturns] ; 

    %%%%%%%%%%%%%%%%
    % magnetic field
    %%%%%%%%%%%%%%%%
%     [x_B, y_B, Bx, By, B, Bxcheck] = magFieldCoil(coils,mu0,positionArrayField,Bmax) ; % [T]
    directory = '/home/emma/Documents/Research/TP11/Field Data' ; 
    filenameGeneric = 'u1_1_9_u2_0' ; 
    nfiles = 5 ;
    rawOn = 0 ; 
    nfig = 1 ; 
    [X,Y,avgdata] = plotBFieldData(directory,filenameGeneric,nfiles,rawOn,nfig) ; 
    [x_B, y_B, Bx, By, B] = magFieldCoil_exp(X, Y, avgdata) ; % [T]
    
%     if rem(i, 2) == 1
%         Bx = flip(Bx, 2) ; 
%         By = flip(By, 2) ; 
%         B = flip(B, 2) ; 
%     end        
    
    temp_B = zeros(size(x_B,1),size(x_B,2),5) ; % create a 3D array 
    
    temp_B(:,:,1) = x_B ; 
    temp_B(:,:,2) = y_B ; 
    temp_B(:,:,3) = Bx ; 
    temp_B(:,:,4) = By ; 
    temp_B(:,:,5) = B ; 
    Barray{i} = temp_B ; 
    
%     figure(1)
%     surf(x_B,y_B,B) ; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % magnetic field gradient 
    %%%%%%%%%%%%%%%%%%%%%%%%%

    n = 26 ; % size of position vector
%     [x_gB, y_gB, gBx, gBy, gB] = magGradientCoil(coils,mu0,positionArrayGrad,n,gBmax) ; 
    [x_gB, y_gB, gBx, gBy, gB] = magGradientCoil_exp(X, Y, avgdata) ; 
    
%     if rem(i, 2) == 1
%         gBx = flip(gBx, 2) ; 
%         gBy = flip(gBy, 2) ; 
%         gB = flip(gB, 2) ; 
%     end        
    
    temp_gB = zeros(size(x_gB,1),size(x_gB,2),5) ; % create a 3D array 
    
    temp_gB(:,:,1) = x_gB ; 
    temp_gB(:,:,2) = y_gB ; 
    temp_gB(:,:,3) = gBx ; 
    temp_gB(:,:,4) = gBy ; 
    temp_gB(:,:,5) = gB ; 
    gBarray{i} = temp_gB ; 
    
%     figure(2)
%     quiver(x_gB,y_gB,gBx,gBy,3) ;
%     hold on
%     axis([-(d+0.01),(d+0.01),-(d+0.01),(d+0.01)])
    
end

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
save('MAT_SingleLinkModelCoil_case1') ; 
