function [gx, gy, gBx, gBy, gB] = magGradientCoil(coils,mu0,positionArray,n,gBmax)   
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
    %   positionArray :- array of values defining the extent of the space
    %   in x and y that the function will use to calculat magnetic field
    %   values
    %   positionArray = [dx ; xmin ; xmax ; ymin ; ymax]
    %
    %   gBmax :- highest value of magnetic field gradient allowed
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
    global a d nturns 
    
    I1 = coils(1,1) ;
    I2 = coils(2,1) ;
%     a = coils(3,1) ; % radius 
%     d = coils(4,1) ; 
%     nturns = coils(5,1) ;  
    
    xmin = positionArray(2,1) ; 
    xmax = positionArray(3,1) ; 
    ymin = positionArray(4,1) ; 
    ymax = positionArray(5,1) ; 
    
    xlin = linspace(xmin,xmax,n) ; 
    ylin = linspace(ymin,ymax,n) ; 
    [x,y] = meshgrid(xlin,ylin) ; 
    
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
    
    gBx = reshape(gBx,length(x),length(y)) ; 
    gBx = [gBx(end:-1:2,:) ; gBx] ; 

    gBy = reshape(gBy,length(x),length(y)) ; 
    gBy = [-gBy(end:-1:2,:) ; gBy] ; 

    gB = reshape(gB,length(x),length(y)) ; 
    gB = fillmissing(gB,'constant',0) ; 
    gB = [gB(end:-1:2,:) ; gB] ; 
    gB(gB>gBmax) = gBmax ;
    
    xval = linspace(xmin,xmax,size(gBx,2)) ; 
    yval = linspace(-ymax,ymax,size(gBx,1)) ; 
    [gx,gy] = meshgrid(xval,yval) ; 

%     [gBx,gBy] = signChanges(gBx,gBy,gx,gy,I1,I2) ;  

end