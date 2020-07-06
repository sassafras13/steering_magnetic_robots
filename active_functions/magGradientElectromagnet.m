function [gx, gy, gBx, gBy, gB] = magGradientElectromagnet(mu0,positionArray,n,gBmax,d,r,I,a,nturns)   
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
    
    % create the grid of points for which we are calculating the field
    xmin = positionArray(2,1) ; 
    xmax = positionArray(3,1) ; 
    ymin = positionArray(4,1) ; 
    ymax = positionArray(5,1) ; 
    
    xlin = linspace(xmin,xmax,n) ; 
    ylin = linspace(ymin,ymax,n) ; 
    [xgrid,ygrid] = meshgrid(xlin,ylin) ; 
    
    % calculate the center locations of the four electromagnets
    EMlocation = zeros(4,2) ; 
    
    EMlocation(1,1) = -d ; 
    EMlocation(1,2) = r ; 
    
    EMlocation(2,1) = -d ; 
    EMlocation(2,2) = -r ; 
    
    EMlocation(3,1) = d ; 
    EMlocation(3,2) = r ; 
    
    EMlocation(4,1) = d ; 
    EMlocation(4,2) = -r ; 

    % calculate gradient contribution from each electromagnet separately
    % then add together
    % this is actually faster because we can calculate the contribution
    % from each electromagnet all in one go for all points, thus minimizing
    % the number of times we have to do symbolic calculations
    
    % electromagnet 1
    x1 = xgrid - EMlocation(1,1) ; 
    y1 = ygrid - EMlocation(1,2) ; 
    rho1 = sqrt( x1^2 + y1^2 ) ; 
    [~,gBrho1] = dBdi(I,a,mu0,nturns,0,rho1) ; 
    theta1 = atan2( y1, x1 ) ; % can I do atan2 of a meshgrid of values?
    gBrho1x = gBrho1 * cos(theta1) ; 
    gBrho1y = gBrho1 * sin(theta1) ; 
    
    % electromagnet 2
    x2 = xgrid - EMlocation(2,1) ; 
    y2 = ygrid - EMlocation(2,2) ; 
    rho2 = sqrt( x2^2 + y2^2 ) ; 
    [~,gBrho2] = dBdi(I,a,mu0,nturns,0,rho2) ; 
    theta2 = atan2( y2, x2 ) ; % can I do atan2 of a meshgrid of values?
    gBrho2x = gBrho2 * cos(theta2) ; 
    gBrho2y = gBrho2 * sin(theta2) ; 
    
    % electromagnet 3
    x3 = xgrid - EMlocation(3,1) ; 
    y3 = ygrid - EMlocation(3,2) ; 
    rho3 = sqrt( x3^2 + y3^2 ) ; 
    [~,gBrho3] = dBdi(I,a,mu0,nturns,0,rho3) ; 
    theta3 = atan2( y3, x3 ) ; % can I do atan2 of a meshgrid of values?
    gBrho3x = gBrho3 * cos(theta3) ; 
    gBrho3y = gBrho3 * sin(theta3) ; 
    
    % electromagnet 4
    x4 = xgrid - EMlocation(4,1) ; 
    y4 = ygrid - EMlocation(4,2) ; 
    rho4 = sqrt( x4^2 + y4^2 ) ; 
    [~,gBrho4] = dBdi(I,a,mu0,nturns,0,rho4) ; 
    theta4 = atan2( y4, x4 ) ; % can I do atan2 of a meshgrid of values?
    gBrho4x = gBrho4 * cos(theta4) ; 
    gBrho4y = gBrho4 * sin(theta4) ; 
    
    % add all contributions at each point together    
    gBx = gBrho1x + gBrho2x + gBrho3x + gBrho4x ; 
    gBy = gBrho1y + gBrho2y + gBrho3y + gBrho4y ;  
    gB = sqrt(gBx.^2 + gBy.^2) ; 
    
    gBx = double(gBx) ; 
    gBy = double(gBy) ; 
    gB = double(gB) ; 
    
    gBx = reshape(gBx,length(xgrid),length(ygrid)) ; 
    gBx = [gBx(end:-1:2,:) ; gBx] ; 

    gBy = reshape(gBy,length(xgrid),length(ygrid)) ; 
    gBy = [-gBy(end:-1:2,:) ; gBy] ; 

    gB = reshape(gB,length(xgrid),length(ygrid)) ; 
    gB = fillmissing(gB,'constant',0) ; 
    gB = [gB(end:-1:2,:) ; gB] ; 
    gB(gB>gBmax) = gBmax ;
    
    xval = linspace(xmin,xmax,size(gBx,2)) ; 
    yval = linspace(-ymax,ymax,size(gBx,1)) ; 
    [gx,gy] = meshgrid(xval,yval) ; 

%     [gBx,gBy] = signChanges(gBx,gBy,gx,gy,I1,I2) ;  

end