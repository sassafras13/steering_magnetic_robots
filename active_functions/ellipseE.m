function E = ellipseE(ksq) 
    %ellipseE calculates the complete elliptic integral of the second kind.
    %It will fail if you pass ksq = 1. 
    
    eE = @(theta) sqrt(1 - ksq.*(sin(theta).^2)) ; 
    
    E = integral(eE,0,(pi/2)) ; 
end