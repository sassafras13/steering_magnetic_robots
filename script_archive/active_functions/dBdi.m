function [gBx,gBy] = dBdi(I,a,mu0,n,xarray,yarray)
    % dBdi Calculate the gradient of the magnetic field at every point in
    % xarray, yarray. This function uses an embedded series of symbolic
    % equations to solve the gradient of the magnetic field and then
    % converts the equation to a series of values which it outputs for
    % plotting in a quiver plot or equivalent.
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
    
    syms x y 'real'
    ksq = (4*a*y) / ((a+y)^2 + x^2) ; 
    [K,E] = ellipke(ksq) ; 
    Bx = ((mu0*I*n) / (2*pi)) * (1 / sqrt((y + a)^2 + x^2)) * (((a^2 - y^2 - x^2) / ((a-y)^2 + x^2))*E + K) ;
    By = ((mu0*I*n) / (2*pi)) * (x / (y * sqrt((y+a)^2 + x^2))) * (((a^2 + y^2 + x^2) / ((a-y)^2 + x^2))*E - K) ;



%     ksq = zeros( (size(xarray,1)*size(xarray,2)), 1 ) ; 
%     Bx = zeros( (size(xarray,1)*size(xarray,2)), 1 ) ; 
%     By = zeros( (size(xarray,1)*size(xarray,2)), 1 ) ; 
%     
%     for i = 1:( (size(xarray,1) * size(xarray,2)) ) 
%         ksq(i,1) = (4*a*yarray(i)) / ((a+yarray(i))^2 + xarray(i)^2) ; 
%         K = ellipseK(ksq(i,1)) ; 
%         E = ellipseE(ksq(i,1)) ; % results match MATLAB's ellipke function
%         Bx(i,1) = ((mu0*I*n) / (2*pi)) * (1 / sqrt((yarray(i) + a)^2 + xarray(i)^2)) * (((a^2 - yarray(i)^2 - xarray(i)^2) / ((a-yarray(i))^2 + xarray(i)^2))*E + K) ;
%         By(i,1) = ((mu0*I*n) / (2*pi)) * (xarray(i) / (yarray(i) * sqrt((yarray(i)+a)^2 + xarray(i)^2))) * (((a^2 + yarray(i)^2 + xarray(i)^2) / ((a-yarray(i))^2 + xarray(i)^2))*E - K) ;
%     end
    B = sqrt(Bx.^2 + By.^2) ; % check if a dot operation is different? --> other way produces 676x676 matrix
%     B = reshape(B,[size(xarray,1), size(xarray,2)]) ; 
%     [gBx,gBy] = gradient(B) ; 
    
    
%     B = subs(B,{x,y},{xarray,yarray}) ; 
%     dBdx = simplify(diff(B,x)) ; 
%     dBdy = simplify(diff(B,y)) ;
    dBdx = diff(B,x) ; 
    dBdy = diff(B,y) ;
    gBx = subs(dBdx,{x,y},{xarray,yarray}) ; 
    gBy = subs(dBdy,{x,y},{xarray,yarray}) ; 
end