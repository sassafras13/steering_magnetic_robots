% VAR_SingleLinkModelCoil.m
% Contains variables for use in MAIN_SingleLinkModelCoil.m

%% microswimmer properties
% w = 1 ; % [mm] width of macroswimmer
% L = 10 ; % [mm] length of macroswimmer
% h = 5 ; % [mm] height of macroswimmer
% mass = 1 ; % [g] mass of macroswimmer
% mag = 255 ; % [Amm2] magnetic moment of macroswimmer (obtained experimentally)

w = 0.001 ; % [m] width of macroswimmer
L = 0.010 ; % [m] length of macroswimmer
h = 0.005 ; % [m] height of macroswimmer
mass = 0.001 ; % [kg] mass of macroswimmer
mag = 2.55E-4 ; % [Am2] magnetic moment of macroswimmer (obtained experimentally)
drag = 0.0471 ; % [] drag coefficient FIT PARAMETER
ct = drag ; 
cn = 2*ct ; 


%% environment properties
mu0 = 4*pi*(10^(-7)) ; % [N/A^2] magnetic permeability

%% coil dimensions 
a = (0.136 / 2) ; % [m] radius of coils
d = sqrt(3)*(a/2) ; % [m] Maxwell separation of coils
nturns = 320 ; % number of turns
% I1 = 1.9 ; % [A]
% I2 = 0 ; % [A]
% coils = [I1 ; I2 ; a ; d ; nturns] ; 

%% position array for calculations and plotting
% field plotting
dx = 0.001; % [m]
xmin = -1.5*d ; % [m]
xmax = 1.5*d ; % [m]
ymin = 0.01 ; % [m]
ymax = 1.5*a ; % [m]
positionArrayField = [dx ; xmin ; xmax ; ymin ; ymax] ; 

% % gradient plotting
% xmingrad = -1.5*d ; % [m]
% xmaxgrad = 1.5*d ; % [m]
% ymingrad = 0.010 ; % [m]
% ymaxgrad = 1.5*d ; % [m]
% positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 

% gradient plotting
xmingrad = -0.9*d ; % [m]
xmaxgrad = 0.9*d ; % [m]
ymingrad = 0.01 ; % [m]
ymaxgrad = 0.9*d ; % [m]
positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 

%% max value of B fields
Bmax = 0.01; % [T] highest value of magnetic field allowed
gBmax = 10;  % [?] highest value of magnetic field gradient allowed