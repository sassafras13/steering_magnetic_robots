function [gx, gy, Bx, By, B, gBx, gBy, gB] = magFieldCoil_fit_v2_forMATfile(c1, c2, I1, I2, a, nturns, mu0, Bmax, x, y)
    % magFieldCoil_fit Calculate the magnetic field created by a pair of 
    % wire coils with current flowing through them via the Biot-Savart law.
    % For use with Matlab fit() function to find optimal fit parameters.
    %
    % Inputs: 
    %   I1 :- current to coil 1
    %
    %   I2 :- current to coil 2
    %
    %   a :- coil radius
    %
    %   nturns :- number of turns in one coil
    % 
    %   mu0 :- magnetic permeability of free space
    %
    %   x :- vector of x positions at which B is sampled
    %
    %   y :- vector of y positions at which B is sampled
    %
    %   Bmax :- highest value of magnetic field (in T) allowed 
    %
    % Outputs: 
    %   B :- magnitude of magnetic field at each position in 2D space as
    %   defined by x, y
    
    
    [gx,gy] = meshgrid(linspace(min(x), max(x), 20),...
        linspace(min(y), max(y), 20)) ; % [m]     
    
    d = sqrt(3)*(a/2) ; % [m] Maxwell separation of coils
        
    Bx = zeros(size(gx,1)*size(gx,2),1) ; 
    By = zeros(size(gx,1)*size(gx,2),1) ; 

    for i = 1:(size(gx,1)*size(gx,2))
        
        % first coil
        xi1 = gx(i) - d  ;
        yi1 = abs(gy(i)  );
        
        ksq1 = ksquared(xi1,yi1,a) ;
        [K1,E1] = ellipke(ksq1) ; 
        Bx1 = BxHH(xi1,yi1,mu0,nturns,I1,a,E1,K1) ;
        By1 = ByHH(xi1,yi1,mu0,nturns,I1,a,E1,K1) ;
        
        % second coil
        xi2 = gx(i) + d  ;
        yi2 = abs(gy(i)  );
        
        ksq2 = ksquared(xi2,yi2,a) ;
        [K2,E2] = ellipke(ksq2) ; 
        Bx2 = BxHH(xi2,yi2,mu0,nturns,I2,a,E2,K2) ;
        By2 = ByHH(xi2,yi2,mu0,nturns,I2,a,E2,K2) ;
        
        % sum 
        Bx(i) = Bx1 + Bx2 ; 
        By(i) = By1 + By2 ; 
        
    end
    
    Bx = reshape(Bx,size(gx,1),size(gx,2)) ; 
    Bx = c1 * Bx ; 
    
    By = reshape(By,size(gx,1),size(gx,2)) ; 
    By = fillmissing(By,'constant',0) ; 
    By = c2 * By ; 
    
    B = sqrt(Bx.^2 + By.^2) ;

    B(B>Bmax) = Bmax ; 
    
    % calculate the gradient
    [gBx, gBy] = gradient(B, gx(1, :), gy(:, 1)) ; % [T]
    
    gB = sqrt( gBx.^2 + gBy.^2 ) ; % [T]

end
