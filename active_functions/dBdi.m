function [gBx,gBy] = dBdi(I,a,mu0,n,xarray,yarray)
    % dBdi Calculate the gradient of the magnetic field at every point in
    % xarray, yarray. This function is entirely numeric, and runs faster
    % than its symbolic counterpart, dBdi_symbolic.m, in the
    % inactive_functions folder.
    %
    % Inputs: 
    %   I :- current to coil
    %
    %   a :- radius of coil
    %   
    %   mu0 :- magnetic permeability of free space
    %
    %   n :- number of turns in coil
    %
    %   xarray, yarray :- array of points in space for which we want to
    %   calculate gradient
    %
    % Outputs: 
    %   gBx, gBy :- components of the gradient in the x and y directions

    dx = abs(xarray(1,1) - xarray(1,2)) ; 
    dy = abs(yarray(1,1) - yarray(2,1)) ; 
    
    ksq = zeros( (size(xarray,1)*size(xarray,2)), 1 ) ; 
    Bx = zeros( (size(xarray,1)*size(xarray,2)), 1 ) ; 
    By = zeros( (size(xarray,1)*size(xarray,2)), 1 ) ; 
    
    for i = 1:( (size(xarray,1) * size(xarray,2)) ) 
        ksq(i,1) = (4*a*yarray(i)) / ((a+yarray(i))^2 + xarray(i)^2) ; 
        K = ellipseK(ksq(i,1)) ; 
        E = ellipseE(ksq(i,1)) ; % results match MATLAB's ellipke function
        Bx(i,1) = ((mu0*I*n) / (2*pi)) * (1 / sqrt((yarray(i) + a)^2 + xarray(i)^2)) * (((a^2 - yarray(i)^2 - xarray(i)^2) / ((a-yarray(i))^2 + xarray(i)^2))*E + K) ;
        By(i,1) = ((mu0*I*n) / (2*pi)) * (xarray(i) / (yarray(i) * sqrt((yarray(i)+a)^2 + xarray(i)^2))) * (((a^2 + yarray(i)^2 + xarray(i)^2) / ((a-yarray(i))^2 + xarray(i)^2))*E - K) ;
    end

    B = sqrt(Bx.^2 + By.^2) ; % check if a dot operation is different? --> other way produces 676x676 matrix
    B = reshape(B,[size(xarray,1), size(xarray,2)]) ; 
    [gBx,gBy] = gradient(B, dx, dy) ; 
%     gBx = 5E2 * gBx ; 
%     gBy = 5E2 * gBy ; 
    
end