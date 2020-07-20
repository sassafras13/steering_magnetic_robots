% MAIN_05_GenerateFieldData.m
% Purpose: Generate an array of values for the magnetic field and magnetic
% field gradient after the optimal fit paramters have been found. 
% Author: Emma Benjaminson

%% 

clc ; clear all ; close all ; 

%% Directories

localDir = fileparts(mfilename('fullpath')) ;
restoredefaultpath ;
addpath(fullfile(localDir, 'active_functions')) ;

%% Load Variables

run(fullfile(localDir,'VAR_SingleLinkModelCoil.m')) ; 

Bmax = 2E-2 ; % [T] highest value of magnetic field allowed

%% Load Variables

%%%%%%%%%% USER DEFINED VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of primitives in trajectory
niter = 3 ; 

% control inputs
I1 = [1.9 ;
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

%% Calculate Magnetic and Gradient Fields

[x_fit, y_fit] = meshgrid(-0.0883:0.005:0.0883, -0.1:0.005:0.1) ; 
x_fit = reshape(x_fit, size(x_fit, 1)*size(x_fit, 2), 1) ; 
y_fit = reshape(y_fit, size(y_fit, 1)*size(y_fit, 2), 1) ; 

Barray = {} ; 
gBarray = {} ; 

for i = 1:niter
    
    % coils
    coils = [I1(i) ; I2(i) ; a ; d ; nturns] ; 

    %%%%%%%%%%%%%%%%
    % magnetic field
    %%%%%%%%%%%%%%%%
    
    % best fit found by fit_3()

    % for u1_1_9_u2_0
    c0 = 0.8026 ; 
    c1 = 2.347 ; 
    c2 = 0.7705 ; 
    c3 = 0.2539 ; 
    c4 = 51.86 ; 
    c5 = -153.3 ; 
    c6 = 1.412 ; 
    c7 = -4.613 ; 
    c8 = 1.052 ; 
    c9 = 30.29 ; 
    c10 = -83.28 ; 
    c11 = -179.9 ; 

    [x_B, y_B, Bx, By, B, gBx, gBy, gB] = ...
        magFieldCoil_fit_v3_forMATfile(c0, c1, c2, c3, c4, c5, c6, c7, ...
        c8, c9, c10, c11, I1(i), I2(i), a, nturns, mu0, Bmax, x_fit, y_fit) ; 

%     % for u1_1_9_u2_0
%     c1 = 0.8014 ; 
%     c2 = 0.2144 ; 
% 
%     [x_B, y_B, Bx, By, B, gBx, gBy, gB] = ...
%         magFieldCoil_fit_v2_forMATfile( c1, c2, I1(i), I2(i), a, nturns,...
%         mu0, Bmax, x_fit, y_fit ) ; 
    
    temp_B = zeros(size(x_B,1),size(x_B,2),5) ; % create a 3D array 
    
    temp_B(:,:,1) = x_B ; 
    temp_B(:,:,2) = y_B ; 
    temp_B(:,:,3) = Bx ; 
    temp_B(:,:,4) = By ; 
    temp_B(:,:,5) = B ; 
    Barray{i} = temp_B ; 
    
    temp_gB = zeros(size(x_B,1),size(x_B,2),5) ; % create a 3D array 
    temp_gB(:,:,1) = x_B ; 
    temp_gB(:,:,2) = y_B ; 
    temp_gB(:,:,3) = gBx ; 
    temp_gB(:,:,4) = gBy ; 
    temp_gB(:,:,5) = gB ; 
    gBarray{i} = temp_gB ; 
%     
%     figure(1)
%     surf(x_B,y_B,B) ; 
%     
    % check the gradient is similar to my model with no fitting
    n = 26 ; % size of position vector
    [x_gB, y_gB, gBxtest, gBytest, gBtest] = magGradientCoil(coils,mu0,positionArrayGrad,n,gBmax) ; 

%     figure(2)
%     quiver(x_gB,y_gB,gBxtest,gBytest,3) ;
%     hold on
%     axis([-(d+0.01),(d+0.01),-(d+0.01),(d+0.01)])
%     hold on
%     quiver(x_B, y_B, gBx, gBy, 3) ; 
%     hold on 
%     legend('Model with No Fitting', 'Model with Fitting')
%     title('Magnetic Field Gradient Check')
%     xlabel('X (m)') ; ylabel('Y (m)') 
%     hold on
%     
%     pause
end

%% Save Data

save('MAT_gBarray.mat', 'gBarray') ; 
save('MAT_Barray.mat', 'Barray') ; 

