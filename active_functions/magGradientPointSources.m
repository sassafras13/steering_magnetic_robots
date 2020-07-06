function [gx, gy, gBx, gBy, gB] = magGradientPointSources(mu0,nturns,positionArray,pointSources)
    % magGradientPointSources Calculate magnetic field gradient for a
    % collection of 4 point sources arranged to mimic the configuration of
    % two coils. 
    %
    % Inputs: 
    %   mu0 :- magnetic permeability of free space
    %
    %   nturns :- number of turns in the equivalent set of coils that we
    %   are representing. Used to increase the magnetic field produced by
    %   the point source to be equivalent to field produced by a coil of
    %   wire with many turns
    %
    %   positionArray :- array of values defining the extent of the space
    %   in x and y that the function will use to calculate magnetic field
    %   values
    %   positionArray = [dx ; xmin ; xmax ; ymin ; ymax] 
    %   
    %   pointSources :- the array of information on the point sources used
    %   to produce the magnetic field gradient
    %   pointSource matrix: 
    %   col 1 = x coordinate [m]
    %   col 2 = y coordinate [m]
    %   col 3 = orientation (+1 = out of the page, -1 = into the page)
    %   col 4 = current [A]
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
    
    dx = positionArray(1,1) ; 
    xmin = positionArray(2,1) ; 
    xmax = positionArray(3,1) ; 
    ymax = positionArray(5,1) ; 
    ymin = -ymax ; 
    
    [gx,gy] = meshgrid(xmin:dx:xmax, ymin:dx:ymax) ; % [m]     
    gBx = zeros(size(gx)) ; % [T]
    gBy = zeros(size(gy)) ; % [T]
    dthetaarray = [] ; 
    gB = [] ; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE BX, BY AT EVERY POINT %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:(size(gx,1)*size(gx,2))    
            for k = 1:size(pointSources,1)
                dx = gx(i) - pointSources(k,1) ; 
                dy = gy(i) - pointSources(k,2) ; 
                dr = sqrt(dx^2 + dy^2) ; 
                dtheta = atan2(dy,dx) ; 
                dthetaarray = [dthetaarray, dtheta] ;

                % check orientation and set n accordingly
                % if orientation is out of page, vectors should point IN
                % if orientation is into page, vectors should point OUT
                if pointSources(k,3) > 0
                    n = 2 ;
                elseif pointSources(k,3) < 0
                    n = 1 ; 
                end
                
                I = pointSources(k,4) ; 

                gBx(i) = gBx(i) + ((-1)^n)*(-1)*((mu0*nturns*I*cos(dtheta)) / (2*pi*(dr^2))) ; 
                gBy(i) = gBy(i) + ((-1)^n)*(-1)*((mu0*nturns*I*sin(dtheta)) / (2*pi*(dr^2))) ; 

            end 
            gB = [gB ; sqrt(gBx(i)^2 + gBy(i)^2)] ; 
    end

    gB = reshape(gB,size(gx,1),size(gx,2)) ; 

end