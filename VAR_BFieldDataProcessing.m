% 00_VAR_BFieldDataProcessing.m
% Contains variables for use in MAIN_BFieldDataProcessing.m

%% environment properties
mu0 = 4*pi*(10^(-7)) ; % [N/A^2] magnetic permeability

%% coil dimensions 
I1 = 1.9 ; % [A] current
I2 = 0 ; % [A] current
a = (0.136 / 2) ; % [m] radius of coils
d = sqrt(3)*(a/2) ; % [m] Maxwell separation of coils
nturns = 320 ; % number of turns
coils = [I1 ; I2 ; a ; d ; nturns] ; 

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