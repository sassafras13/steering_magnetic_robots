% MAIN_BFieldDataProcessing.m
% Purpose: Process the magnetic field data and save the average values and
% standard deviations to a csv file for future use.
% Author: Emma Benjaminson
% References: 

%% 

clc ; clear all ; close all ; 

%% Directories

localDir = fileparts(mfilename('fullpath')) ;
restoredefaultpath ;
addpath(fullfile(localDir, 'active_functions')) ;

%% load variables

run VAR_BFieldDataProcessing.m ; 

%% process experimental data

% import all csvs for a particular configuration
directory = '/home/emma/Documents/Research/TP11/Field Data' ; 
filenameGeneric = 'u1_1_9_u2_0' ; 
% filenameGeneric = 'u1_1_9_u2_1_9' ; 
% filenameGeneric = 'u1_1_9_u2_-1_9' ; 

nfiles = 5 ;

% calculate standard deviation and average values for each point
processBFieldData(directory,filenameGeneric,nfiles) ; 

% %% plot all raw data
% 
% rawOn = 1 ; % show all data 
% nfig = 1 ; 
% 
% % plot all the data (raw and avg +/- std dev) in a 3D plot
% plotBFieldData(directory,filenameGeneric,nfiles,rawOn,nfig) ; 
% 
% % format plot for raw data display
% xlabel('X (cm)','Interpreter','latex','FontSize',12) ; 
% ylabel('Y (cm)','Interpreter','latex','FontSize',12) ; 
% zlabel('B (G)','Interpreter','latex','FontSize',12) ; 
% title('Coil Configuration: U1 = 1.9A, U2 = 1.9A','Interpreter','latex','FontSize',14) ; 


%% compare to B field model

% only show average data
rawOn = 0 ; 
nfig = 2 ; 
[X,Y,avgdata] = plotBFieldData(directory,filenameGeneric,nfiles,rawOn,nfig) ; 

% magnetic field
Bmax = 2E-2 ; % [T] highest value of magnetic field allowed
[gx, gy, Bx, By, B, Bxcheck] = magFieldCoil(coils,mu0,positionArrayField,Bmax) ; % [T]
% B = B * 1E3; % convert to mT
% B = B * 10 ; % convert to G

figure(2)
surf(gx,gy,B,'FaceAlpha',0.5,'EdgeColor','k','FaceColor','#52639e')
hold on
ax = gca ; 
ax.FontSize = 14 ; 
xlabel('Z (m)','Interpreter','latex','FontSize',32) ; 
ylabel('X (m)','Interpreter','latex','FontSize',32) ; 
zlabel('B (T)','Interpreter','latex','FontSize',32) ; 
% title('Magnetic Field Produced by Maxwell Coils','Interpreter','latex','FontSize',32) ; 
axis([0.035,0.045,-0.15,0.15,0,0.01]) ; 
% legend('Experimental Data','Model','interpreter','latex','FontSize',20) ; 
grid on

% %% compare to B field gradient model 
% 
% % magnetic field gradient
% n = 14 ; % number of points in position vector for gradient calculation
% gBmax = 10 ; % [?] highest value of magnetic field gradient allowed
% [gx2, gy2, gBx, gBy, gB] = magGradientCoil(coils,mu0,positionArrayGrad,n,gBmax) ; 
% 
% % import average B field data 
% rawOn = 0 ; 
% X = reshape(X,[14,11]) ; 
% Y = reshape(Y,[14,11]) ; 
% avgdata = reshape(avgdata(:,1),[14,11]) ; 
% 
% % calculate the gradient in x and y
% [dx,dy] = gradient( avgdata ) ; 
% 
% % plot trajectory
% directory = '/home/emma/Documents/Research/TP11/Sequence 1 v2' ; 
% [indices, origins] = extractIndOr(directory) ; 
% filename = 's1' ; 
% outputname = 's1-CompositeData' ; 
% [Xq,meanposy,stdposy] = processExpData(indices, origins, ...
%     directory, filename, outputname) ; 
% 
% figure(3)
% quiver(gx2,gy2,gBx,gBy,3,'-b')
% hold on 
% quiver(X,Y,dx,dy,3,'--r') 
% hold on
% xlabel('X [m]','interpreter','latex','FontSize',12) ; 
% ylabel('Y [m]','interpreter','latex','FontSize',12) ;  
% title('Magnetic Field Gradient Produced by Maxwell Coils','Interpreter','latex') ; 
% legend('Exp. Trajectory Data','Model','Exp. Field Data','interpreter','latex') ; 