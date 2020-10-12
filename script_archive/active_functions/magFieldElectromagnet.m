 function [gx, gy, Bx, By, B] = magFieldElectromagnet(mu0,r,d,positionArray,Bmax,nturns,I,a,z)
    % magFieldCoil Calculate the magnetic field created by a configuration
    % of electromagnets that are arranged to represent an equivalent coil
    % system.
    %
    % Inputs: 
    %   mu0 :- magnetic permeability of free space
    %
    %   r :- the radius of the coils we are trying to represent with the
    %   electromagnet configuration
    %
    %   d :- the separation of the equivalent coil system
    %
    %   positionArray :- array of values defining the extent of the space in
    %   x and y that the function will use to calculate magnetic field values
    %   positionArray = [dx ; xmin ; xmax ; ymin ; ymax]
    %
    %   Bmax :- highest value of magnetic field (in T) allowed 
    %
    %   nturns :- the number of turns in the electromagnets
    %
    %   I :- the current in the electromagnets
    %   
    %   a :- the radius of the electromagnets
    %
    %   z :- the height along the electromagnet at which we are modeling
    %   the field
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

    % create the grid of points for which we are calculating the field
    dx = positionArray(1,1) ; 
    xmin = positionArray(2,1) ; 
    xmax = positionArray(3,1) ; 
    ymin = positionArray(4,1) ; 
    ymax = positionArray(5,1) ; 
    
    [gx,gy] = meshgrid(xmin:dx:xmax, ymin:dx:ymax) ; % [m]  
    
    % calculate the center locations of the four electromagnets
    EMlocation = zeros(2,2) ; 
    
%     EMlocation(1,1) = -d ; 
%     EMlocation(1,2) = r ; 
    
    EMlocation(1,1) = -d ; 
    EMlocation(1,2) = -r ; 
    
    EMlocation(2,1) = d ; 
    EMlocation(2,2) = r ; 
    
%     EMlocation(4,1) = d ; 
%     EMlocation(4,2) = -r ; 
        
    % create arrays for magnetic field in x and y directions
    Bx = zeros(size(gx,1)*size(gx,2),1) ; 
    By = zeros(size(gx,1)*size(gx,2),1) ; 

    for i = 1:(size(gx,1)*size(gx,2))
        
        % observation point
        x = gx(i) ; 
        y = gy(i) ; 
        
        % calculate contributions for each electromagnet
        for j = 1:2 
            
            % get the coordinates for the current electromagnet
            xj = EMlocation(j,1) ; 
            yj = EMlocation(j,2) ; 
            
            % calculate the field and its angle at the observation point
            rho = sqrt( (x-xj)^2 + (y-yj)^2 ) ; 
            k = ksquared(z,rho,a) ; 
            [K,E] = ellipke(k) ; 
            Brho = ByHH(z,rho,mu0,nturns,I,a,E,K) ; % y = rho 
            theta = atan2( (y-yj), (x-xj) ) ; 
            
            % add the field from the current electromagnet to the
            % contributions from the other electromagnets at this point
            Bx(i) = Bx(i) + Brho.*cos(theta) ; 
            By(i) = By(i) + Brho.*sin(theta) ; 
        end
        
    end
    
    Bx = reshape(Bx,size(gx,1),size(gx,2)) ; 
%     Bx = [Bx(end:-1:2,:) ; Bx] ; 
    Bx(Bx>Bmax) = Bmax ; 
    
    By = reshape(By,size(gx,1),size(gx,2)) ; 
    By = fillmissing(By,'constant',0) ; 
%     By = [By(end:-1:2,:) ; By] ; 
    By(By>Bmax) = Bmax ; 
    
    B = sqrt(Bx.^2 + By.^2) ;
    B(B>Bmax) = Bmax ; 
%     B = sign(Bx).*B ; 
    
%     xval = linspace(xmin,xmax,size(Bx,2)) ; 
%     yval = linspace(-ymax,ymax,size(Bx,1)) ; 
%     [gx,gy] = meshgrid(xval,yval) ; 

end
