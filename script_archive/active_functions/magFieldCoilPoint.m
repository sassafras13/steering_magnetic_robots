function [Bx, By, B] = magFieldCoilPoint(coils,mu0,x,y,Bmax)
    % magFieldCoil Calculate the magnetic field created by a pair of wire
    % coils with current flowing through them via the Biot-Savart law. 
    %
    % Inputs: 
    %   coils :- data on the coils, including current in each coil, radius,
    %   spacing between coils and number of turns
    %   coils = [I1 ; I2 ; a ; d ; nturns]
    % 
    %   mu0 :- magnetic permeability of free space
    %
    %   positionArray :- array of values defining the extent of the space in
    %   x and y that the function will use to calculate magnetic field values
    %   positionArray = [dx ; xmin ; xmax ; ymin ; ymax]
    %
    %   Bmax :- highest value of magnetic field (in T) allowed 
    %
    % Outputs: 
    %   gx :- meshgrid of x values, to be used for plotting results
    %
    %   gy :- meshgrid of y values, to be used for plotting results
    %
    %   Bx :- x-components of magnetic field at each position in 2D space
    %   as defined by positionArray
    %
    %   By :- y-components of magnetic field at each position in 2D space
    %   as defined by positionArray
    %
    %   B :- magnitude of magnetic field at each position in 2D space as
    %   defined by positionArray
    %
    %   Bxcheck :- analytical solution of magnetic field along center axis
    %   between coils. Only valid for Helmholtz coil configuration! 
    
    I1 = coils(1,1) ; % [A] current
    I2 = coils(2,1) ; % [A] current
    a = coils(3,1) ; % [m] radius 
    d = coils(4,1) ; % [m] separation distance
    nturns = coils(5,1) ; % number of turns
        
    x1 = x - d ;
    ksq1 = ksquared(x1,y,a) ;
    [K1,E1] = ellipke(ksq1) ; 
    Bx1 = BxHH(x1,y,mu0,nturns,I1,a,E1,K1) ;
    By1 = ByHH(x1,y,mu0,nturns,I1,a,E1,K1) ;
        
    % second coil
    x2 = x + d ;         
    ksq2 = ksquared(x2,y,a) ;
    [K2,E2] = ellipke(ksq2) ; 
    Bx2 = BxHH(x2,y,mu0,nturns,I2,a,E2,K2) ;
    By2 = ByHH(x2,y,mu0,nturns,I2,a,E2,K2) ;
        
    % sum 
    Bx = Bx1 + Bx2 ; 
    By = By1 + By2 ; 
      
    Bx(Bx>Bmax) = Bmax ; 
    By(By>Bmax) = Bmax ; 
    
    B = sqrt(Bx.^2 + By.^2) ;
    B(B>Bmax) = Bmax ; 

end
