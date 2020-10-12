% MAIN_03_GenerateFieldData.m
% Purpose: Generate an array of values for the magnetic field and magnetic
% field gradient. 
% Author: Emma Benjaminson

%% Load Variables

%%%%%%%%%% USER DEFINED VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of primitives in trajectory
niter = 3 ; 

% control inputs
I1 = [1.9 ;
      0 ; 
      1.9 ];  % [A] current
I2 = [0 ; 
      1.9 ; 
      0 ]; % [A] current
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
    
    [gx, gy, Bx, By, B, Bxcheck] = magFieldCoil(coils,mu0,positionArrayField,Bmax) ; % [T]
    temp_B = zeros(size(gx,1),size(gx,2),5) ; % create a 3D array 
    
    temp_B(:,:,1) = gx ; 
    temp_B(:,:,2) = gy ; 
    temp_B(:,:,3) = Bx ; 
    temp_B(:,:,4) = By ; 
    temp_B(:,:,5) = B ; 
    Barray{i} = temp_B ; 
    
    %%%%%%%%%%%%%%%%%%%
    % magnetic gradient
    %%%%%%%%%%%%%%%%%%%
    
    n = 26 ; % size of position vector
    [gx2, gy2, gBx, gBy, gB] = magGradientCoil(coils,mu0,positionArrayGrad,n,gBmax) ; 
    temp_gB = zeros(size(gx2,1),size(gx2,2),5) ; % create a 3D array 
    temp_gB(:,:,1) = gx2 ; 
    temp_gB(:,:,2) = gy2 ; 
    temp_gB(:,:,3) = gBx ; 
    temp_gB(:,:,4) = gBy ; 
    temp_gB(:,:,5) = gB ; 
    gBarray{i} = temp_gB ; 

end

%% Save Data

save('MAT_gBarray.mat', 'gBarray') ; 
save('MAT_Barray.mat', 'Barray') ; 

