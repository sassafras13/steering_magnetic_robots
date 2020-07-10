% MAIN_04_FitCoilModeltoData.m
% Purpose: Fit the model of the Maxwell coils and their resultant magnetic
% field to the experimental data.
% Author: Emma Benjaminson

%% 

clc ; clear all ; close all ; 

%% Directories

localDir = fileparts(mfilename('fullpath')) ;
restoredefaultpath ;
addpath(fullfile(localDir, 'active_functions')) ;

%% Load Variables

run(fullfile(localDir,'VAR_SingleLinkModelCoil.m')) ; 

%% Setup

% coil dimensions 
I1 = 1.9 ; % [A] current
I2 = 0 ; % [A] current
coils = [I1 ; I2 ; a ; d ; nturns] ; 
Bmax = 2E-2 ; % [T] highest value of magnetic field allowed

%% Get Experimental Data

% import all csvs for a particular configuration
directory = '/home/emma/Documents/Research/TP11/Field Data' ; 
filenameGeneric = 'u1_1_9_u2_0' ; 
% filenameGeneric = 'u1_1_9_u2_1_9' ; 
% filenameGeneric = 'u1_1_9_u2_-1_9' ; 

nfiles = 5 ;

% calculate standard deviation and average values for each point
processBFieldData(directory,filenameGeneric,nfiles) ; 

rawOn = 1 ; % show all data 
nfig = 1 ; 

% plot the average data 
[x,y,avgdata] = plotBFieldData(directory,filenameGeneric,nfiles,rawOn,nfig) ; 
B = avgdata(:, 1) ; 

%% Use Optimizer to Match Model to Data 

ft = fittype( 'magFieldCoil_fit( I1, I2, a, nturns, mu0, Bmax, x, y )',...
    'independent', {'x', 'y'}, 'dependent', 'B_final',...
    'problem', {'mu0', 'Bmax'})

f = fit( [x, y], B, ft, 'problem', {mu0, Bmax},...
    'StartPoint', [1.9, 0, 0.0680, 320])

ft2 = fittype( 'magFieldCoil_fit_v2( c1, c2, I1, I2, a, nturns, mu0, Bmax, x, y )', ...
    'independent', {'x', 'y'}, 'dependent', 'B_final', ...
    'problem', {'I1', 'I2', 'a', 'nturns', 'mu0', 'Bmax'})

f2 = fit( [x, y], B, ft, 'probem', {I1, I2, a, nturns, mu0, Bmax})
%% Plot Results

% best fit found by fit()
[x_fit, y_fit] = meshgrid(-0.0883:0.005:0.0883, -0.1:0.005:0.1) ; 
x_fit = reshape(x_fit, size(x_fit, 1)*size(x_fit, 2), 1) ; 
y_fit = reshape(y_fit, size(y_fit, 1)*size(y_fit, 2), 1) ; 
B_fit = magFieldCoil_fit( 1.901, -0.05156, 0.08196, 320, mu0, Bmax,...
    x_fit, y_fit ) ; 

% best fit found by fit_2()
B_fit2 = magFieldCoil_fit_v2(c1, c2, I1, I2, a, nturns, mu0, Bmax, x, y)

% model with no fit
[gx, gy, Bx, By, B_model, Bxcheck] = magFieldCoil(coils,mu0,...
    positionArrayField,Bmax) ; % [T]

figure(2)
plot3(x, y, B, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r') % mean experimental data
hold on
surf(reshape(x_fit, 41, 36), reshape(y_fit, 41, 36), reshape(B_fit, 41, 36),...
    'FaceAlpha',0.5,'EdgeColor','k','FaceColor','#7AEA86') % best fit
hold on
surf(gx, gy, B_model, 'FaceAlpha',0.5,'EdgeColor','k',...
    'FaceColor','#52639e') % model with no fit
hold on
legend('Mean Exp Data', 'Model with Best Fit', 'Model with No Fit')
xlabel('Z (m)','Interpreter','latex','FontSize',32) ; 
ylabel('X (m)','Interpreter','latex','FontSize',32) ; 
zlabel('B (T)','Interpreter','latex','FontSize',32) ; 
