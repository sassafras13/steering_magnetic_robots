% VARIABLES.m
% Purpose: Stores all fixed variable values used in this paper.
% Author: Emma Benjaminson

global w L h mass mag drag ct cn mu0 a d nturns Bmax gBmax 
%% microswimmer properties
w = 0.001 ; % [m] width of macroswimmer
L = 0.010 ; % [m] length of macroswimmer
h = 0.005 ; % [m] height of macroswimmer
mass = 0.001 ; % [kg] mass of macroswimmer
mag = 2.55E-4 ; % [Am2] magnetic moment of macroswimmer (obtained experimentally)
drag = 0.0471 ; % [-] drag coefficient FIT PARAMETER
ct = drag ; 
cn = 2*ct ; 

%% environment properties
mu0 = 4*pi*(10^(-7)) ; % [N/A^2] magnetic permeability

%% coil dimensions 
a = (0.136 / 2) ; % [m] radius of coils
d = sqrt(3)*(a/2) ; % [m] Maxwell separation of coils
nturns = 320 ; % number of turns

%% max value of B fields
Bmax = 0.01; % [T] highest value of magnetic field allowed
gBmax = 10;  % [] highest value of magnetic field gradient allowed

%% position array for calculations and plotting
% field plotting
dx = 0.005 ; % [m]
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