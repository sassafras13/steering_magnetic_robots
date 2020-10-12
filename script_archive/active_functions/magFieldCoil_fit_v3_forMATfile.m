function [gx, gy, Bx, By, B, gBx, gBy, gB] = magFieldCoil_fit_v3_forMATfile(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, I1, I2, a, nturns, mu0, Bmax, x, y)
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
%         Bx(i) = (c0 + c1*gx(i) + c2*(gy(i)^2)) * (Bx1 + Bx2) ; 
%         By(i) = (c3*gy(i) + c4*gx(i)*gy(i) + c5*(gy(i)^2)) * (By1 + By2) ; 

%         Bx(i) = (c0 + c1*gx(i)) * (Bx1 + Bx2) ; 
%         By(i) = (c2 + c3*gy(i)) * (By1 + By2) ; 

%         Bx(i) = (c0 + c1*gx(i) + c2*(gx(i)^2)) * (Bx1 + Bx2) ; 
%         By(i) = (c3 + c4*gy(i) + c5*(gy(i)^2)) * (By1 + By2) ; 

%         Bx(i) = (c0 + c1*gx(i) + c2*gy(i)) * (Bx1 + Bx2) ; 
%         By(i) = (c3 + c4*gx(i) + c5*gy(i)) * (By1 + By2) ; 

%         Bx(i) = (c0*gx(i) + c1*(gx(i)^2)) * (Bx1 + Bx2) ; 
%         By(i) = (c2*gy(i) + c3*(gy(i)^2)) * (By1 + By2) ; 

        Bx(i) = (c0 + c1*gx(i) + c2*gy(i) + c3*gx(i)*gy(i) + c4*(gx(i)^2) + c5*(gy(i)^2)) * (Bx1 + Bx2) ; 
        By(i) = (c6 + c7*gx(i) + c8*gy(i) + c9*gx(i)*gy(i) + c10*(gx(i)^2) + c11*(gy(i)^2)) * (By1 + By2) ;
        
    end
    
    Bx = reshape(Bx,size(gx,1),size(gx,2)) ; 
    
    By = reshape(By,size(gx,1),size(gx,2)) ; 
    By = fillmissing(By,'constant',0) ; 
    
    B = sqrt(Bx.^2 + By.^2) ;
    
    B(B>Bmax) = Bmax ; 

    
%     B = interp2(gx, gy, B, x, y) ; 
    
    % calculate the gradient
    [gBx, gBy] = gradient(B, gx(1, :), gy(:, 1)) ; % [T]
    
    gB = sqrt( gBx.^2 + gBy.^2 ) ; % [T]

end
