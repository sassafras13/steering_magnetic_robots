% MAIN_01_BFieldDataProcessing.m
% Purpose: Computes the mean and standard deviation for a set of data 
% characterizing the magnetic field. Saves the values to a csv file. 
% Plots the data and compares to theoretical values.
% Author: Emma Benjaminson


%% Process Experimental Data

% import all csvs for a particular configuration
directory = fullfile(localDir, '/experimental_data') ; 
filenameGeneric = 'u1_1_9_u2_0' ; 
% filenameGeneric = 'u1_1_9_u2_1_9' ; 
% filenameGeneric = 'u1_1_9_u2_-1_9' ; 

nfiles = 5 ;

% calculate standard deviation and average values for each point
processBFieldData(directory,filenameGeneric,nfiles) ; 

%% Figure 1: Plot all raw data

rawOn = 1 ; % show all data 
nfig = 1 ; 

% plot all the data (raw and avg +/- std dev) in a 3D plot
plotBFieldData(directory,filenameGeneric,nfiles,rawOn,nfig) ; 

% format plot for raw data display
xlabel('X (cm)','Interpreter','latex','FontSize',12) ; 
ylabel('Y (cm)','Interpreter','latex','FontSize',12) ; 
zlabel('B (G)','Interpreter','latex','FontSize',12) ; 
title('Coil Configuration: U1 = 1.9A, U2 = 0A','Interpreter','latex','FontSize',14) ; 


%% Figure 2: Compare to B field model

rawOn = 0 ; % only show average data
nfig = 2 ; 

% plot the average data 
[X,Y,avgdata] = plotBFieldData(directory,filenameGeneric,nfiles,rawOn,nfig) ; 

% compute the model of the magnetic field
Bmax = 2E-2 ; % [T] highest value of magnetic field allowed
I1 = 1.9 ; 
I2 = 0 ; 
coils = [I1 ; I2 ; a ; d ; nturns] ; 

[gx, gy, Bx, By, B, Bxcheck] = magFieldCoil(coils,mu0,positionArrayField,Bmax) ; % [T]
% B = B * 1E3; % convert to mT
% B = B * 10 ; % convert to G

figure(nfig)
surf(gx,gy,B,'FaceAlpha',0.5,'EdgeColor','k','FaceColor','#52639e')
hold on
ax = gca ; 
ax.FontSize = 14 ; 
xlabel('Z (m)','Interpreter','latex','FontSize',32) ; 
ylabel('X (m)','Interpreter','latex','FontSize',32) ; 
zlabel('B (T)','Interpreter','latex','FontSize',32) ; 
title('Magnetic Field Produced by Maxwell Coils','Interpreter','latex','FontSize',32) ; 
% axis([0.035,0.045,-0.15,0.15,0,0.01]) ; 
legend('Experimental Data','Model','interpreter','latex','FontSize',20) ; 
grid on

%% Figure 3: Compare to B field gradient model 

% magnetic field gradient
n = 14 ; % number of points in position vector for gradient calculation
gBmax = 10 ; % [mT] highest value of magnetic field gradient allowed

% compute the gradient from the model
[gx2, gy2, gBx, gBy, ~] = magGradientCoil(coils,mu0,positionArrayGrad,n,gBmax) ; 

% compute the gradient from the experimental data
[x_gB, y_gB, gBx1, gBy1, gB] = magGradientCoil_exp(X, Y, avgdata) ; 

% check the gradient in x and y
X = reshape(X,[14,11]) ; 
Y = reshape(Y,[14,11]) ; 
avgdata = reshape(avgdata(:,1),[14,11]) ; 
[dx,dy] = gradient( avgdata ) ; 

figure(3)
quiver(gx2,gy2,gBx,gBy,3,'-b') % model
hold on 
quiver(x_gB, y_gB, gBx1, gBy1,3,'-r') % data with new function
hold on
xlabel('X [m]','interpreter','latex','FontSize',12) ; 
ylabel('Y [m]','interpreter','latex','FontSize',12) ;  
title('Magnetic Field Gradient Produced by Maxwell Coils','Interpreter','latex') ; 
legend('Model','Exp. Field Data','interpreter','latex') ; 