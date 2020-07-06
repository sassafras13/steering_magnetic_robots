function [gBx, gBy, gB] = magGradientElectromagnetSinglePoint(mu0,gBmax,z,I,a,nturns,x,y)   
    % magGradientCoil Calculate the magnetic field gradient created by a configuration
    % of electromagnets that are arranged to represent an equivalent coil
    % system.
    %
    % Inputs: 
    %   mu0 :- magnetic permeability of free space
    %
    %   positionArray :- array of values defining the extent of the space
    %   in x and y that the function will use to calculat magnetic field
    %   values
    %   positionArray = [dx ; xmin ; xmax ; ymin ; ymax]
    %
    %   n :- number of points to calculate in linspace formula
    %
    %   gBmax :- highest value of magnetic field gradient allowed
    %
    %   d :- the separation of the equivalent coils sytem
    %
    %   r :- the radius of the coils we are trying to represent with the
    %   electromagnet configuration
    %
    %   I :- the current in the electromagnets
    %
    %   a :- the radius of the electromagnets
    %
    %   nturns :- the number of turns in the electromagnets
    %
    % Outputs: 
    %   gx :- meshgrid of x values, to be used for plotting results
    %
    %   gy :- meshgrid of y values, to be used for plotting results
    %
    %   gBx :- x-components of magnetic field gradient at each position in
    %   2D space as defined by positionArray
    %
    %   gBy :- y-components of magnetic field gradient at each position in
    %   2D space as defined by positionArray
    %
    %   gB :- magnitude of magnetic field gradient at each position in 2D
    %   space as defined by positionArray

    % calculate gradient contribution from each electromagnet separately
    rho1 = sqrt( x.^2 + y.^2 ) ; 
    [~,gBrho1] = dBdi(I,a,mu0,nturns,z,rho1) ;
    gBrho = double(gBrho1) ; 
    theta = atan2( y, x ) ; 
    gBx = gBrho .* cos(theta) ; 
    gBy = gBrho .* sin(theta) ; 
    gB = sqrt(gBx.^2 + gBy.^2) ; 
    gB = fillmissing(gB,'constant',0) ; 
    gB(gB>gBmax) = gBmax ; 

end