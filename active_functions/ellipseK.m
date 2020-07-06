function K = ellipseK(ksq)
    %ellipseK calculates the complete elliptic integral of the first kind.
    %It will fail if you pass ksq = 1. 
    
    eK = @(theta) 1 ./ sqrt(1 - ksq.*(sin(theta).^2)) ; 
    
    K = integral(eK,0,(pi/2)) ; 
end