function [gx, gy, Bx, By, B, Bxcheck] = magFieldCoil(coils,mu0,positionArray,Bmax)
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

    dx = positionArray(1,1) ; 
    xmin = positionArray(2,1) ; 
    xmax = positionArray(3,1) ; 
    ymin = positionArray(4,1) ; 
    ymax = positionArray(5,1) ; 
    
    [gx,gy] = meshgrid(xmin:dx:xmax, ymin:dx:ymax) ; % [m]     
    
    I1 = coils(1,1) ; % [A] current
    I2 = coils(2,1) ; % [A] current
    a = coils(3,1) ; % [m] radius 
    d = coils(4,1) ; % [m] separation distance
    nturns = coils(5,1) ; % number of turns
    
    [coeff1, coeff2] = config(I1, I2) ; 
    
    Bx = zeros(size(gx,1)*size(gx,2),1) ; 
    By = zeros(size(gx,1)*size(gx,2),1) ; 

    for i = 1:(size(gx,1)*size(gx,2))
        
        % first coil
        xi1 = gx(i) - d ; 
        yi1 = gy(i) ; 
        
        ksq1 = ksquared(xi1,yi1,a) ;
        [K1,E1] = ellipke(ksq1) ; 
        Bx1 = BxHH(xi1,yi1,mu0,nturns,I1,a,E1,K1) ;
        By1 = ByHH(xi1,yi1,mu0,nturns,I1,a,E1,K1) ;
        
        % second coil
        xi2 = gx(i) + d ; 
        yi2 = gy(i) ; 
        
        ksq2 = ksquared(xi2,yi2,a) ;
        [K2,E2] = ellipke(ksq2) ; 
        Bx2 = BxHH(xi2,yi2,mu0,nturns,I2,a,E2,K2) ;
        By2 = ByHH(xi2,yi2,mu0,nturns,I2,a,E2,K2) ;
        
        % sum 
        Bx(i) = Bx1 + Bx2 ; 
        By(i) = By1 + By2 ; 
        
    end
    
    Bx = reshape(Bx,size(gx,1),size(gx,2)) ; 
    Bx = [Bx(end:-1:2,:) ; Bx] ; 
    Bx(Bx>Bmax) = Bmax ; 
    
    By = reshape(By,size(gx,1),size(gx,2)) ; 
    By = fillmissing(By,'constant',0) ; 
    By = [By(end:-1:2,:) ; By] ; 
    By(By>Bmax) = Bmax ; 
    
    B = sqrt(Bx.^2 + By.^2) ;
    B(B>Bmax) = Bmax ; 
%     B = sign(Bx).*B ; 
    
    xval = linspace(xmin,xmax,size(Bx,2)) ; 
    yval = linspace(-ymax,ymax,size(Bx,1)) ; 
    [gx,gy] = meshgrid(xval,yval) ; 
    
    % test along the x axis where we know the solution
    Bxcheck = zeros(size(gx,2),1) ; 
    for i = 1:size(gx,2)
        x = gx(1,i) ; 
        Bxcheck(i) = 0.5*nturns*mu0*I1*(a^2)*(coeff1*(1 /(sqrt(a^2 + (x-d)^2))^3) + coeff2*(1 /(sqrt(a^2 + (x+d)^2))^3)) ;
%         Bxcheck(i) = 0.5*nturns*mu0*I1*(a^2)*((1 /(sqrt(a^2 + (x-d)^2))^3) + (1 /(sqrt(a^2 + (x+d)^2))^3)) ;
    end

end
