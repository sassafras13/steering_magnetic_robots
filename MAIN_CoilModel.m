% MAIN_CoilModel.m
% Purpose: Characterizes the magnetic field and magnetic field gradient
% produced by a pair of Helmholtz or Maxwell coils. For comparison with other models to
% prove it is a robust model. 
% Author: Emma Benjaminson

%% 

clc ; clear all ; close all ; 

%% Directories

localDir = fileparts(mfilename('fullpath')) ;
restoredefaultpath ;
addpath(fullfile(localDir, 'active_functions')) ;

%% Setup

mu0 = 4*pi*(10^(-7)) ; % [N/A^2] magnetic permeability

% coil dimensions 
I1 = 1.9 ; % [A] current
I2 = 0 ; % [A] current
a = (0.136 / 2) ; % [m] radius of coils
% d = a/2 ; % [m] Helmholtz separation of coils
d = sqrt(3)*(a/2) ; % [m] Maxwell separation of coils
nturns = 320 ; % number of turns
coils = [I1 ; I2 ; a ; d ; nturns] ; 

% position array for calculations and plotting
% field plotting
dx = 0.001 ; % [m]
xmin = -1.5*d ; % [m]
xmax = 1.5*d ; % [m]
ymin = 0.01 ; % [m]
ymax = 1.5*a ; % [m]
positionArrayField = [dx ; xmin ; xmax ; ymin ; ymax] ; 

% gradient plotting
xmingrad = -0.9*d ; % [m]
xmaxgrad = 0.9*d ; % [m]
ymingrad = 0.01 ; % [m]
ymaxgrad = 0.9*d ; % [m]
positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 

%% Function Calls

%%%%%%%%%%%%%%%%
% magnetic field
%%%%%%%%%%%%%%%%
Bmax = 2E-2 ; % [T] highest value of magnetic field allowed
[gx, gy, Bx, By, B, Bxcheck] = magFieldCoil(coils,mu0,positionArrayField,Bmax) ; % [T]
B = B * 1E3; % convert to mT
Bxcheck = Bxcheck * 1E3; % convert to mT

%%%%%%%%%%%%%%%%%%%%%%%%%
% magnetic field gradient 
%%%%%%%%%%%%%%%%%%%%%%%%%

n = 26 ; % size of position vector
gBmax = 10 ; % [?] highest value of magnetic field gradient allowed
[gx2, gy2, gBx, gBy, gB] = magGradientCoil(coils,mu0,positionArrayGrad,n,gBmax) ; 
% gB = gB*1E6 ; % convert 

%% Experimental Data 

% board-fridge
selpath = '/home/emma/Documents/Research/TP7' ; 
filetype = 'All_Board-Fridge*.csv' ; 
names = dir(sprintf('%s/%s',selpath,filetype)); 
fullfilename1 = fullfile(selpath,names.name); 
A = readtable(fullfilename1) ; 

posx1 = table2array(A(:,1)) ; % [cm]
alphalabs1 = table2array(A(:,3)) ; % [G]

posx1 = posx1 ./ 100 ; % [m]
alphalabs1 = alphalabs1 ./ 10 ; % [mT]

posy1 = zeros(length(posx1),1) ; 
accuracy_alphalabs = 0.02 ; % [%] mfg reported accuracy on sensor
error_alphalabs1 = accuracy_alphalabs * alphalabs1 ; % [mT] 

%% Plots

%%%%%%%%%%%%%%%%
% magnetic field
%%%%%%%%%%%%%%%%
% plot magnetic field magnitude as surface plot
figure(1)
surf(gx,gy,B)
hold on
plot3(posx1',posy1',alphalabs1','or','MarkerFaceColor','r','MarkerSize',4) % experimental data
hold on 
plot3([posx1,posx1]',[posy1,posy1]',[-error_alphalabs1,error_alphalabs1]'+alphalabs1','-g','LineWidth',4) ; % error on experimental data
xlabel('X [m]','interpreter','latex','FontSize',12) ; 
ylabel('Y [m]','interpreter','latex','FontSize',12) ; 
zlabel('B [mT]','interpreter','latex','FontSize',12) ; 
title('Magnitude of Magnetic Field Produced by Maxwell Coils','Interpreter','latex') ; 
legend('Model','Experimental Data') ; 

%%%%%%%%%%%%%%%%%%%%%%%%%
% magnetic field gradient
%%%%%%%%%%%%%%%%%%%%%%%%%

% plot magnetic field as contour lines
figure(2)
% contour(gx,gy,B,30)
% hold on
quiver(gx2,gy2,gBx,gBy,3)
hold on
axis([-(d+0.01),(d+0.01),-(d+0.01),(d+0.01)])
xlabel('X [m]','interpreter','latex','FontSize',12) ; 
ylabel('Y [m]','interpreter','latex','FontSize',12) ;  
title('Magnetic Field Gradient Produced by Maxwell Coils','Interpreter','latex') ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check field against analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3) 
plot(gx(1,:),Bxcheck,'--b') 
hold on 
plot(gx(1,:),B(ceil(size(B,1)/2),:),'-r')
xlabel('Z [m]','interpreter','latex','FontSize',12) ; 
ylabel('B [mT]','interpreter','latex','FontSize',12) ;
% title(sprintf('Maxwell Coils, I1 = %0.1fA, I2 = %0.1fA\n r = %0.2fm, separation = %0.2fm',I1,I2,a,2*d) ,'FontSize',14) ; 
title('Magnitude of Magnetic Field Along Centerline of Maxwell Coils','Interpreter','latex') ; 
legend('Analytical Solution','Model','interpreter','latex','FontSize',12,'location','southeast') ; 

