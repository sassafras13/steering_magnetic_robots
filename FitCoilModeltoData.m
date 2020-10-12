% FitCoilModeltoData.m
% Purpose: Fit the model of the Maxwell coils and their resultant magnetic
% field to the experimental data using a multivariate polynomial model.
% Author: Emma Benjaminson

%% Clear Workspace

clc ; clear all ; close all ; 

%% Directories

localDir = fileparts(mfilename('fullpath')) ;
restoredefaultpath ;
addpath(fullfile(localDir, 'active_functions')) ;
addpath(fullfile(localDir, 'experimental_data')) ; 

%% Load Variables

run(fullfile(localDir,'VARIABLES.m')) ; 

%% Setup

% coil dimensions 
I1 = 1.9 ; % [A] current
I2 = 0 ; % [A] current
coils = [I1 ; I2 ; a ; d ; nturns] ; 
Bmax = 2E-2 ; % [T] highest value of magnetic field allowed

%% Get Experimental Data

% import all csvs for a particular configuration
directory = fullfile(localDir, '/experimental_data') ; 
filenameGeneric = 'u1_1_9_u2_0' ; 
% filenameGeneric = 'u1_1_9_u2_1_9' ; 
% filenameGeneric = 'u1_1_9_u2_-1_9' ; 

nfiles = 5 ; % for u1_1_9_u2_0
% nfiles = 4 ; % for u1_1_9_u2_1_9

% calculate standard deviation and average values for each point
processBFieldData(directory,filenameGeneric,nfiles) ; 

rawOn = 1 ; % show all data 
nfig = 1 ; 

% plot the average data 
[x,y,avgdata] = plotBFieldData(directory,filenameGeneric,nfiles,rawOn,nfig) ; 
B = avgdata(:, 1) ; 

%% Use Optimizer to Match Model to Data 

ft = fittype( 'magFieldCoil_fit(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, I1, I2, a, nturns, mu0, Bmax, x, y)',...
    'independent', {'x', 'y'}, 'dependent', 'B_final', ...
    'problem', {'I1', 'I2', 'a', 'nturns', 'mu0', 'Bmax'}) ;

f = fit( [x, y], B, ft, 'problem', {I1, I2, a, nturns, mu0, Bmax}) ;

%% Plot Results

% best fit found by fit()
[x_fit, y_fit] = meshgrid(-0.0883:0.005:0.0883, -0.1:0.005:0.1) ; 
x_fit = reshape(x_fit, size(x_fit, 1)*size(x_fit, 2), 1) ; 
y_fit = reshape(y_fit, size(y_fit, 1)*size(y_fit, 2), 1) ; 

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

B_fit = magFieldCoil_fit(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, I1, I2, a, nturns, mu0, Bmax, x_fit, y_fit) ; 

% model with no fit
[gx, gy, Bx, By, B_model, Bxcheck] = magFieldCoil(coils,mu0,...
    positionArrayField,Bmax) ; % [T]

figure(2)
plot3(x, y, B, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r') % mean experimental data
hold on
surf(gx, gy, B_model, 'FaceAlpha',0.5,'EdgeColor','k',...
    'FaceColor','#52639e') % model with no fit
hold on
surf(reshape(x_fit, 41, 36), reshape(y_fit, 41, 36), reshape(B_fit, 41, 36),...
    'FaceAlpha',0.5,'EdgeColor','k','FaceColor','#EA7AE3') % best fit
hold on
legend('Mean Exp Data', 'Model with No Fit', 'Model with Scaled Fit')
xlabel('Z (m)','Interpreter','latex','FontSize',32) ; 
ylabel('X (m)','Interpreter','latex','FontSize',32) ; 
zlabel('B (T)','Interpreter','latex','FontSize',32) ; 
