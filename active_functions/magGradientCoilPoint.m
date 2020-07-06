function [gBx, gBy, gB] = magGradientCoilPoint(coils,mu0,x,y)   
    % magGradientCoil Calculate the magnetic field gradient created by a
    % pair of wire coils with current flowing through them via the
    % Biot-Savart law. 
    %
    % Inputs: 
    %   coils :- data on the coils, including current in each coil, radius,
    %   spacing between coils and number of turns
    %   coils = [I1 ; I2 ; a ; d ; nturns]
    %
    %   mu0 :- magnetic permeability of free space
    %
    %   x, y :- coordinates of the point at which to calculate the field
    %   gradient
    %
    % Outputs: 
    %   gBx :- x-components of magnetic field gradient at each position in
    %   2D space as defined by positionArray
    %
    %   gBy :- y-components of magnetic field gradient at each position in
    %   2D space as defined by positionArray
    %
    %   gB :- magnitude of magnetic field gradient at each position in 2D
    %   space as defined by positionArray
    
    I1 = coils(1,1) ;
    I2 = coils(2,1) ;
    a = coils(3,1) ; % radius 
    d = coils(4,1) ; 
    nturns = coils(5,1) ;  
    
    x1 = x - d ; 
    [gBx1,gBy1] = dBdi(I1,a,mu0,nturns,x1,y) ; 
    
    x2 = x + d ; 
    [gBx2,gBy2] = dBdi(I2,a,mu0,nturns,x2,y) ; 
    
    gBx = gBx1 + gBx2 ; 
    gBy = gBy1 + gBy2 ; 
    gB = sqrt(gBx.^2 + gBy.^2) ; 
    
    gBx = double(gBx) ; 
    gBy = double(gBy) ; 
    gB = double(gB) ; 


end